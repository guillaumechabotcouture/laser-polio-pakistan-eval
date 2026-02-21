#!/usr/bin/env python3
"""
Chad Guinea Worm Dual-Host SEIS Spatial Model

Models Dracunculus medinensis transmission across 8 endemic districts in
southern Chad with coupled human and dog populations sharing water sources.

Key biology:
  - SEIS dynamics: no lasting immunity (recovered → susceptible)
  - Water-mediated transmission via copepods in stagnant water
  - 365-day pre-patent period (E), 21-day worm emergence (I)
  - Dogs are dominant reservoir (~95% of detected infections in Chad)
  - Dry-season peak (Jan-Mar) when stagnant pools concentrate copepods

Architecture:
  Two parallel LASER Models (human, dog) coupled via DualHostWaterTransmission.
  Each model has its own within-species gravity-coupled TransmissionSE.
  Cross-species FOI flows through shared water contamination.

Usage:
    /opt/anaconda3/bin/python3 scripts/guinea_worm_chad.py
"""

import sys
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance as haversine_distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from guinea_worm_components import (
    SEISRecovery, DualHostWaterTransmission,
    PatchImportation, build_chad_dry_season,
)

# ============================================================================
# District Configuration — 8 endemic districts in southern Chad
# ============================================================================

DISTRICTS = pd.DataFrame([
    # name,              human_pop, dog_pop, lat,     lon
    ("Chari-Baguirmi",   450000,    40000,   12.11,   15.05),
    ("Moyen-Chari",      400000,    35000,   9.15,    18.45),
    ("Salamat",          350000,    30000,   10.97,   20.67),
    ("Mandoul",          300000,    25000,   8.63,    17.47),
    ("Logone_Oriental",  250000,    20000,   8.62,    16.08),
    ("Tandjile",         200000,    15000,   9.65,    16.72),
    ("Mayo-Kebbi_Est",   150000,    10000,   9.37,    15.35),
    ("Lac",              100000,     5000,   13.30,   14.10),
], columns=["name", "human_pop", "dog_pop", "lat", "lon"])

# ============================================================================
# Disease Parameters
# ============================================================================

# Guinea worm epidemiological parameters
GW_PARAMS_COMMON = {
    "prng_seed": 42,
    "nticks": 5 * 365,           # 5-year simulation
    # Pre-patent period: 10-14 months → ~365 days (gamma distribution)
    "exp_shape": 12,             # gamma shape (lower variance at shape=12)
    "exp_scale": 30.0,           # gamma scale → mean = 12 × 30 = 360 days
    # Worm emergence (infectious): 2-4 weeks → ~21 days
    "inf_mean": 21,
    "inf_sigma": 5,
    # Spatial coupling via shared water sources
    "gravity_b": 0.3,           # destination population exponent (water access)
    "gravity_c": 2.0,           # distance decay (water sources are local)
    # Vital dynamics
    "capacity_safety_factor": 1.5,  # low growth expected in 5yr
}

HUMAN_PARAMS = PropertySet({
    **GW_PARAMS_COMMON,
    "beta": 0.02,               # R0 ≈ 0.4 within-species (subcritical alone)
    "cbr": 44.0,                # Chad CBR per 1000/year (high fertility)
    "gravity_k": 0.005,         # human spatial coupling strength
})

DOG_PARAMS = PropertySet({
    **GW_PARAMS_COMMON,
    "beta": 0.10,               # R0 ≈ 2.1 within-species (dogs sustain epidemic)
    "cbr": 80.0,                # dog crude birth rate per 1000/year
    "gravity_k": 0.003,         # dogs less mobile than humans
})

# Cross-species coupling via shared water sources.
# Dogs are the dominant reservoir (~95% of cases in Chad).
# dog_to_human: dog prevalence spills into human FOI (moderate — most
#   humans use filtered water, so coupling is attenuated)
# human_to_dog: human prevalence spills into dog FOI (very low — humans
#   contribute little water contamination relative to dogs)
DOG_TO_HUMAN_COUPLING = 0.15   # attenuated by filter use / behavior
HUMAN_TO_DOG_COUPLING = 0.05   # humans contribute very little to dog FOI


# ============================================================================
# Build Scenarios
# ============================================================================

def build_scenario(pop_col, seed_patches, seed_counts):
    """Build GeoDataFrame scenario for one species.

    Args:
        pop_col: column name in DISTRICTS ('human_pop' or 'dog_pop')
        seed_patches: list of patch indices for initial infections
        seed_counts: list of initial I counts per seeded patch
    """
    geometry = [Point(lon, lat) for lat, lon in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(DISTRICTS.copy(), geometry=geometry, crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))
    scenario["population"] = scenario[pop_col].astype(np.uint32)

    # Initialize: all susceptible except seeded infections
    scenario["E"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["I"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["R"] = np.zeros(len(scenario), dtype=np.uint32)

    for patch_idx, n_infected in zip(seed_patches, seed_counts):
        scenario.loc[patch_idx, "I"] = n_infected

    scenario["S"] = (scenario["population"] - scenario["E"]
                     - scenario["I"] - scenario["R"]).astype(np.uint32)

    # Validate initial conditions
    total = scenario.S + scenario.E + scenario.I + scenario.R
    assert (total == scenario.population).all(), \
        "Initial S+E+I+R must equal population in every patch"

    return scenario


def build_distance_matrix(scenario):
    """Compute pairwise haversine distance matrix (km)."""
    nnodes = len(scenario)
    lats = np.array(scenario.lat)
    lons = np.array(scenario.lon)
    dist_matrix = np.zeros((nnodes, nnodes))
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                dist_matrix[i, j] = haversine_distance(
                    lats[i], lons[i], lats[j], lons[j]
                )
    return dist_matrix


def build_gravity_network(scenario, dist_matrix, params):
    """Build gravity-model migration network for water-source coupling."""
    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(pops, dist_matrix, 1, 0, params.gravity_b, params.gravity_c)

    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * params.gravity_k
    network = row_normalizer(network, 0.15)  # cap at 15% export per node

    # Validate
    assert network.sum() > 0, "Network is all zeros — check gravity_k"
    assert network.sum(axis=1).max() < 0.3, \
        f"Max row sum = {network.sum(axis=1).max():.3f} — too high"

    return network


def build_age_pyramid(mean_life):
    """Build simple exponential age pyramid."""
    ages = np.arange(100) if mean_life > 20 else np.arange(20)
    max_age = len(ages)
    stable_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_dist = np.maximum(stable_dist, 1)
    return AliasedDistribution(stable_dist), stable_dist


def initialize_ages(model, pyramid):
    """Set initial agent DOBs from age pyramid before model.run()."""
    count = model.people.count
    max_years = len(pyramid[1])  # length of the distribution array
    ages_years = pyramid[0].sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, max_years - 1)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


# ============================================================================
# Build and Run Model
# ============================================================================

def run_model():
    """Build and run the dual-host guinea worm SEIS model."""

    nnodes = len(DISTRICTS)
    nticks = HUMAN_PARAMS.nticks

    # --- Seed initial infections ---
    # Near-elimination: 2 infected humans, 50 infected dogs
    # Concentrated in the 3 most endemic districts (indices 0, 1, 2)
    human_seed_patches = [0, 1, 2]
    human_seed_counts = [1, 1, 0]      # 2 humans total

    dog_seed_patches = [0, 1, 2]
    dog_seed_counts = [20, 20, 10]     # 50 dogs total

    # --- Build scenarios ---
    human_scenario = build_scenario("human_pop", human_seed_patches, human_seed_counts)
    dog_scenario = build_scenario("dog_pop", dog_seed_patches, dog_seed_counts)

    # --- Distance matrix (shared — same geography) ---
    dist_matrix = build_distance_matrix(human_scenario)

    # --- Duration distributions (shared across species) ---
    expdurdist = dists.gamma(shape=GW_PARAMS_COMMON["exp_shape"],
                             scale=GW_PARAMS_COMMON["exp_scale"])
    infdurdist = dists.normal(loc=GW_PARAMS_COMMON["inf_mean"],
                              scale=GW_PARAMS_COMMON["inf_sigma"])

    # --- Seasonality (shared — same climate) ---
    seasonality_map, season_365 = build_chad_dry_season(nticks, nnodes)

    # --- Birth/death rates ---
    human_cbr = HUMAN_PARAMS.cbr  # per 1000/year
    assert 1 <= human_cbr <= 60, f"Human CBR={human_cbr} out of range for per-1000/year"
    human_birthrates = np.full((nticks, nnodes), human_cbr, dtype=np.float32)
    human_cdr = 13.0  # Chad CDR ~13/1000/year
    human_deathrates = np.full((nticks, nnodes), human_cdr, dtype=np.float32)

    dog_cbr = DOG_PARAMS.cbr
    assert 1 <= dog_cbr <= 200, f"Dog CBR={dog_cbr} out of expected range"
    dog_birthrates = np.full((nticks, nnodes), dog_cbr, dtype=np.float32)
    dog_cdr = 100.0  # dogs: ~10-year average lifespan
    dog_deathrates = np.full((nticks, nnodes), dog_cdr, dtype=np.float32)

    # --- Age pyramids ---
    human_pyramid = build_age_pyramid(mean_life=55.0)  # Chad life expectancy
    dog_pyramid = build_age_pyramid(mean_life=8.0)     # dog life expectancy

    # ==== BUILD HUMAN MODEL ====
    print("Building human model...")
    human_model = Model(human_scenario, HUMAN_PARAMS, birthrates=human_birthrates)
    human_model.network = build_gravity_network(human_scenario, dist_matrix, HUMAN_PARAMS)

    # Human components (SEIS: Susceptible → Exposed → Infectious → Susceptible)
    h_susceptible = SEIR.Susceptible(human_model)
    h_exposed = SEIR.Exposed(human_model, expdurdist, infdurdist)
    h_infectious = SEIR.Infectious(human_model, infdurdist)
    h_recovered = SEIR.Recovered(human_model)  # counts R (will be emptied by SEISRecovery)
    h_seis_recovery = SEISRecovery(human_model)  # I→R→S (no lasting immunity)

    human_model.components = [
        h_susceptible,
        h_exposed,
        h_infectious,
        h_recovered,
        h_seis_recovery,           # SEIS: move R back to S
        BirthsByCBR(human_model, birthrates=human_birthrates, pyramid=human_pyramid[0]),
        MortalityByCDR(human_model, mortalityrates=human_deathrates),
        SEIR.Transmission(human_model, expdurdist, seasonality=seasonality_map),
    ]

    initialize_ages(human_model, human_pyramid)

    # ==== BUILD DOG MODEL ====
    print("Building dog model...")
    dog_model = Model(dog_scenario, DOG_PARAMS, birthrates=dog_birthrates)
    dog_model.network = build_gravity_network(dog_scenario, dist_matrix, DOG_PARAMS)

    d_susceptible = SEIR.Susceptible(dog_model)
    d_exposed = SEIR.Exposed(dog_model, expdurdist, infdurdist)
    d_infectious = SEIR.Infectious(dog_model, infdurdist)
    d_recovered = SEIR.Recovered(dog_model)
    d_seis_recovery = SEISRecovery(dog_model)

    # Dog seasonality (same profile)
    dog_seasonality_map, _ = build_chad_dry_season(nticks, nnodes)

    dog_model.components = [
        d_susceptible,
        d_exposed,
        d_infectious,
        d_recovered,
        d_seis_recovery,
        BirthsByCBR(dog_model, birthrates=dog_birthrates, pyramid=dog_pyramid[0]),
        MortalityByCDR(dog_model, mortalityrates=dog_deathrates),
        SEIR.Transmission(dog_model, expdurdist, seasonality=dog_seasonality_map),
    ]

    initialize_ages(dog_model, dog_pyramid)

    # ==== CROSS-SPECIES COUPLING ====
    cross_host = DualHostWaterTransmission(
        human_model, dog_model, expdurdist,
        dog_to_human=DOG_TO_HUMAN_COUPLING,
        human_to_dog=HUMAN_TO_DOG_COUPLING,
        seasonality_values=np.tile(season_365, nticks // 365 + 1)[:nticks],
    )

    # ==== IMPORTATION (near-elimination seeding) ====
    endemic_patches = [0, 1, 2]  # Chari-Baguirmi, Moyen-Chari, Salamat
    human_importation = PatchImportation(
        human_model, expdurdist, endemic_patches,
        period=180, count=1, end_tick=3 * 365,  # 1 case every 6 months for 3 years
    )
    dog_importation = PatchImportation(
        dog_model, expdurdist, endemic_patches,
        period=90, count=3, end_tick=3 * 365,   # 3 dogs every 90 days for 3 years
    )

    # ==== RUN COUPLED SIMULATION ====
    print(f"\nRunning {nticks // 365}-year dual-host guinea worm SEIS simulation...")
    print(f"  Districts: {nnodes}")
    print(f"  Human pop: {human_scenario.population.sum():,}")
    print(f"  Dog pop: {dog_scenario.population.sum():,}")
    print(f"  Human R0≈{HUMAN_PARAMS.beta * GW_PARAMS_COMMON['inf_mean']:.1f}, "
          f"Dog R0≈{DOG_PARAMS.beta * GW_PARAMS_COMMON['inf_mean']:.1f}")
    print(f"  Cross-species: dog→human={DOG_TO_HUMAN_COUPLING}, "
          f"human→dog={HUMAN_TO_DOG_COUPLING}")
    print(f"  Seasonality: dry-season peak (Jan-Mar 1.8x, rainy Jun-Sep 0.4x)")

    # Manual tick loop: run both models + cross-host coupling each tick.
    # Must call _initialize_flows() to propagate compartment counts (S/E/I/R)
    # from tick t to tick t+1 before components modify them — this is what
    # Model.run() does internally.
    import laser.core.random as lcr
    lcr.seed(HUMAN_PARAMS.prng_seed)

    for tick in range(nticks):
        # Initialize flows: copy S[t]→S[t+1], E[t]→E[t+1], etc.
        human_model._initialize_flows(tick)
        dog_model._initialize_flows(tick)

        # Step each model's components
        for component in human_model.components:
            component.step(tick)
        for component in dog_model.components:
            component.step(tick)

        # Cross-species transmission via shared water sources
        cross_host.step(tick)

        # Importation (near-elimination seeding)
        human_importation.step(tick)
        dog_importation.step(tick)

        if tick % 365 == 0:
            year = tick // 365
            h_I = human_model.nodes.I[tick].sum()
            d_I = dog_model.nodes.I[tick].sum()
            print(f"  Year {year}: Human I={h_I}, Dog I={d_I}")

    print("Simulation complete.")
    return human_model, dog_model, human_scenario, dog_scenario, season_365


# ============================================================================
# Analysis and Output
# ============================================================================

def extract_annual_cases(model, nticks):
    """Extract annual new infection counts per patch."""
    nnodes = model.nodes.count
    nyears = nticks // 365

    # newly_infected is computed by TransmissionSE
    newly_inf = model.nodes.newly_infected[:nticks, :]

    annual = np.zeros((nyears, nnodes), dtype=np.int64)
    for yr in range(nyears):
        t0 = yr * 365
        t1 = (yr + 1) * 365
        annual[yr] = newly_inf[t0:t1].sum(axis=0)

    return annual


def print_summary(human_model, dog_model, human_scenario, dog_scenario):
    """Print annual case summary for both species."""
    nticks = human_model.params.nticks
    nyears = nticks // 365

    h_annual = extract_annual_cases(human_model, nticks)
    d_annual = extract_annual_cases(dog_model, nticks)

    print("\n" + "=" * 70)
    print("ANNUAL GUINEA WORM CASES — CHAD DUAL-HOST MODEL")
    print("=" * 70)

    # National totals
    print(f"\n{'Year':<6} {'Human Cases':>12} {'Dog Cases':>12} {'Dog:Human':>12}")
    print("-" * 44)
    for yr in range(nyears):
        h_total = h_annual[yr].sum()
        d_total = d_annual[yr].sum()
        ratio = f"{d_total / max(h_total, 1):.0f}:1" if h_total > 0 else "N/A"
        print(f"{yr + 1:<6} {h_total:>12,} {d_total:>12,} {ratio:>12}")

    # Per-district breakdown
    print(f"\n{'District':<20}", end="")
    for yr in range(nyears):
        print(f"  {'Yr' + str(yr + 1) + ' H':>7}  {'Yr' + str(yr + 1) + ' D':>7}", end="")
    print()
    print("-" * (20 + nyears * 18))

    for i in range(len(human_scenario)):
        name = human_scenario.iloc[i]["name"]
        print(f"{name:<20}", end="")
        for yr in range(nyears):
            print(f"  {h_annual[yr, i]:>7}  {d_annual[yr, i]:>7}", end="")
        print()


def plot_results(human_model, dog_model, season_365):
    """Generate diagnostic plots for the dual-host model."""
    nticks = human_model.params.nticks
    nnodes = human_model.nodes.count

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Chad Guinea Worm Dual-Host SEIS Model", fontsize=14, fontweight="bold")

    # 1. Seasonal forcing profile
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, season_365, "b-", linewidth=2)
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    ax.fill_between(days, season_365, alpha=0.2)
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Transmission Multiplier")
    ax.set_title("Dry-Season Transmission Forcing")
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    month_ticks = [i * 365 // 12 for i in range(12)]
    ax.set_xticks(month_ticks)
    ax.set_xticklabels(months, rotation=45, fontsize=8)

    # 2. Daily incidence (both species)
    ax = axes[0, 1]
    h_daily = human_model.nodes.newly_infected[:nticks, :].sum(axis=1)
    d_daily = dog_model.nodes.newly_infected[:nticks, :].sum(axis=1)
    weeks = np.arange(0, nticks - 6, 7)
    h_weekly = np.array([h_daily[w:w + 7].sum() for w in weeks])
    d_weekly = np.array([d_daily[w:w + 7].sum() for w in weeks])
    ax.plot(weeks / 365, h_weekly, "b-", label="Human", alpha=0.8)
    ax.plot(weeks / 365, d_weekly, "r-", label="Dog", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly New Infections")
    ax.set_title("Epidemic Curves (Weekly)")
    ax.legend()

    # 3. Compartment dynamics — humans
    ax = axes[1, 0]
    t = np.arange(nticks)
    h_S = human_model.nodes.S[:nticks, :].sum(axis=1)
    h_E = human_model.nodes.E[:nticks, :].sum(axis=1)
    h_I = human_model.nodes.I[:nticks, :].sum(axis=1)
    h_N = h_S + h_E + h_I
    ax.plot(t / 365, h_S / np.maximum(h_N, 1), label="S", alpha=0.8)
    ax.plot(t / 365, h_E / np.maximum(h_N, 1), label="E", alpha=0.8)
    ax.plot(t / 365, h_I / np.maximum(h_N, 1), label="I", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Fraction")
    ax.set_title("Human Compartments (SEIS)")
    ax.legend()

    # 4. Compartment dynamics — dogs
    ax = axes[1, 1]
    d_S = dog_model.nodes.S[:nticks, :].sum(axis=1)
    d_E = dog_model.nodes.E[:nticks, :].sum(axis=1)
    d_I = dog_model.nodes.I[:nticks, :].sum(axis=1)
    d_N = d_S + d_E + d_I
    ax.plot(t / 365, d_S / np.maximum(d_N, 1), label="S", alpha=0.8)
    ax.plot(t / 365, d_E / np.maximum(d_N, 1), label="E", alpha=0.8)
    ax.plot(t / 365, d_I / np.maximum(d_N, 1), label="I", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Fraction")
    ax.set_title("Dog Compartments (SEIS)")
    ax.legend()

    plt.tight_layout()
    outpath = Path(__file__).parent.parent / "eval" / "guinea-worm" / "outputs"
    outpath.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath / "guinea_worm_dual_host.png", dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to {outpath / 'guinea_worm_dual_host.png'}")
    plt.close(fig)


# ============================================================================
# Verification
# ============================================================================

def verify_seis_model(human_model, dog_model):
    """Run basic verification checks for the SEIS dual-host model."""
    nticks = human_model.params.nticks
    checks = []

    for label, model in [("Human", human_model), ("Dog", dog_model)]:
        # 1. R compartment should be near zero (SEIS moves R→S each tick)
        R_final = model.nodes.R[nticks - 1, :].sum()
        checks.append((f"{label} SEIS (R≈0)", R_final < 10,
                        f"R={R_final} (should be ~0 for SEIS)"))

        # 2. Population should be non-negative
        S = model.nodes.S[:nticks, :]
        E = model.nodes.E[:nticks, :]
        I = model.nodes.I[:nticks, :]
        neg_S = (S < 0).any()
        neg_E = (E < 0).any()
        neg_I = (I < 0).any()
        checks.append((f"{label} non-negativity", not (neg_S or neg_E or neg_I),
                        f"S<0:{neg_S}, E<0:{neg_E}, I<0:{neg_I}"))

        # 3. Some infections should have occurred
        total_inf = model.nodes.newly_infected[:nticks, :].sum()
        checks.append((f"{label} infections>0", total_inf > 0,
                        f"Total infections: {total_inf:,}"))

        # 4. Spatial coupling — check timing varies across patches
        peak_ticks = []
        for p in range(model.nodes.count):
            ts = model.nodes.newly_infected[:nticks, p]
            if ts.sum() > 0:
                peak_ticks.append(np.argmax(ts))
        if len(peak_ticks) > 1:
            spread = max(peak_ticks) - min(peak_ticks)
            checks.append((f"{label} spatial spread", spread > 0,
                            f"Peak tick spread: {spread} days"))

    print("\n" + "=" * 60)
    print("VERIFICATION REPORT")
    print("=" * 60)
    all_pass = True
    for name, passed, detail in checks:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  [{status}] {name}: {detail}")

    if all_pass:
        print("\nAll verification checks passed.")
    else:
        print("\nSome checks failed — review model configuration.")

    return all_pass


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    human_model, dog_model, human_scenario, dog_scenario, season_365 = run_model()
    print_summary(human_model, dog_model, human_scenario, dog_scenario)
    verify_seis_model(human_model, dog_model)
    plot_results(human_model, dog_model, season_365)
