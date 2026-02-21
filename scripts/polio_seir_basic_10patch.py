#!/usr/bin/env python3
"""
Basic 10-Patch Polio SEIR Spatial Model — LASER Framework

Spatial SEIR simulation for WPV1 transmission dynamics across 10 patches
(districts), each with 100,000 population.

Parameters:
    R0 ~ 6 (beta = R0/D_inf = 6/28 ~ 0.2143)
    Latent period: 3 days (gamma, shape=3, scale=1)
    Infectious period: 28 days (normal, mean=28, sigma=3)
    Initial: 95% recovered (immune), 5 infectious per patch, rest susceptible
    Crude birth rate: 29/1000/year (Pakistan average)
    Crude death rate: 7/1000/year
    Monsoon seasonal forcing (Jul-Oct peak, +/-30%)
    Gravity-model spatial coupling between patches
    Endemic importation in 3 patches to prevent stochastic fadeout

Usage:
    /opt/anaconda3/bin/python3 scripts/polio_seir_basic_10patch.py
"""

import sys
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance as haversine_distance

# Add scripts dir to path for custom_components
sys.path.insert(0, str(Path(__file__).parent))
from custom_components import PatchImportation

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# Configuration
# ============================================================================

N_PATCHES = 10
POP_PER_PATCH = 100_000
N_SEED_INFECTIOUS = 5  # per patch

PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": 10 * 365,            # 10-year simulation
    "beta": 6.0 / 28.0,            # R0=6, D_inf=28d → beta = R0/D ≈ 0.2143
    "exp_shape": 3,                 # Latent period ~3 days (gamma shape=3, scale=1)
    "exp_scale": 1.0,
    "inf_mean": 28,                 # Infectious period ~28 days
    "inf_sigma": 3,
    "cbr": 29.0,                    # Crude birth rate per 1000/year (Pakistan avg)
    "gravity_k": 0.005,             # Spatial coupling strength
    "gravity_b": 0.5,               # Destination population exponent
    "gravity_c": 1.5,               # Distance decay exponent
    "cdr": 7.0,                     # Crude death rate per 1000/year
    "capacity_safety_factor": 2.5,  # Pre-allocate 2.5x initial pop for growth
})

# Patch layout: 10 districts spread across Pakistan (approximate coordinates)
PATCH_NAMES = [f"District_{i+1:02d}" for i in range(N_PATCHES)]
LATS = np.array([34.0, 33.7, 31.5, 30.2, 31.0, 30.5, 31.5, 27.5, 25.4, 24.9])
LONS = np.array([71.5, 73.0, 74.4, 67.0, 70.5, 72.0, 73.1, 65.5, 68.4, 67.0])
INITIAL_IMMUNE_FRAC = 0.95


# ============================================================================
# Seasonal Forcing — Monsoon peak (Jul-Oct)
# ============================================================================

def build_monsoon_seasonality(nticks, nnodes):
    """365-day cosine seasonal profile with monsoon peak (early Sep).

    Pakistan polio cases peak during/after monsoon season (Jul-Oct),
    driven by poor sanitation and population displacement.
    """
    days = np.arange(365)
    peak_day = 245       # Early September
    amplitude = 0.30     # +/-30% modulation around baseline
    beta_season = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
    season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
    return ValuesMap.from_timeseries(season_tiled, nnodes), beta_season


# ============================================================================
# Build and Run
# ============================================================================

def build_scenario():
    """Create GeoDataFrame with 10 equal-population patches."""
    data = {
        "nodeid": range(N_PATCHES),
        "name": PATCH_NAMES,
        "population": [POP_PER_PATCH] * N_PATCHES,
        "lat": LATS,
        "lon": LONS,
    }
    geometry = [Point(lon, lat) for lat, lon in zip(LATS, LONS)]
    scenario = gpd.GeoDataFrame(data, geometry=geometry, crs="EPSG:4326")

    # Initial compartments: 95% recovered, 5 infectious, rest susceptible
    scenario["I"] = np.full(N_PATCHES, N_SEED_INFECTIOUS, dtype=np.uint32)
    scenario["E"] = np.zeros(N_PATCHES, dtype=np.uint32)
    scenario["R"] = np.full(N_PATCHES, int(INITIAL_IMMUNE_FRAC * POP_PER_PATCH), dtype=np.uint32)
    scenario["S"] = (scenario["population"] - scenario["E"]
                     - scenario["I"] - scenario["R"]).astype(np.uint32)

    return scenario


def build_age_pyramid():
    """Approximate Pakistan age distribution (young population, ~65yr life expectancy)."""
    mean_life = 65.0
    ages = np.arange(100)
    stable_age_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_age_dist = np.maximum(stable_age_dist, 1)
    return AliasedDistribution(stable_age_dist)


def initialize_ages(model, pyramid):
    """Assign realistic ages from pyramid so agents aren't all born at tick 0."""
    count = model.people.count
    ages_years = pyramid.sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, 99)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def run_model():
    """Build, configure, and run the 10-patch polio SEIR model."""
    scenario = build_scenario()
    nnodes = len(scenario)
    nticks = PARAMS.nticks

    # Verify initial compartments sum correctly
    for _, row in scenario.iterrows():
        assert row["S"] + row["E"] + row["I"] + row["R"] == row["population"], \
            f"Compartment sum mismatch for {row['name']}"

    # Duration distributions
    expdurdist = dists.gamma(shape=PARAMS.exp_shape, scale=PARAMS.exp_scale)
    infdurdist = dists.normal(loc=PARAMS.inf_mean, scale=PARAMS.inf_sigma)

    # Monsoon seasonal forcing
    seasonality, beta_season_365 = build_monsoon_seasonality(nticks, nnodes)

    # Birth rates — LASER expects per-1000/year (BirthsByCBR divides by 1000 internally)
    assert PARAMS.cbr >= 1 and PARAMS.cbr <= 60, \
        f"CBR must be per-1000/year (got {PARAMS.cbr})"
    birthrate_array = np.full((nticks, nnodes), PARAMS.cbr, dtype=np.float32)

    # Death rates — LASER expects per-1000/year (MortalityByCDR divides by 1000 internally)
    deathrate_array = np.full((nticks, nnodes), PARAMS.cdr, dtype=np.float32)

    # Age pyramid for births
    pyramid = build_age_pyramid()

    # Build model
    model = Model(scenario, PARAMS, birthrates=birthrate_array)

    # Gravity network (manual setup for custom normalization)
    lats = np.array(scenario.lat)
    lons = np.array(scenario.lon)
    dist_matrix = np.zeros((nnodes, nnodes))
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                dist_matrix[i, j] = haversine_distance(
                    lats[i], lons[i], lats[j], lons[j]
                )

    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(pops, dist_matrix, 1, 0, PARAMS.gravity_b, PARAMS.gravity_c)
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * PARAMS.gravity_k
    network = row_normalizer(network, 0.2)  # Cap at 20% export per node
    model.network = network

    # Endemic importation in 3 patches to prevent stochastic fadeout.
    # Initial R_eff = R0 x 0.05 = 0.3 < 1, so disease would die out without
    # importation; births rebuild the susceptible pool over time.
    endemic_patches = [0, 4, 9]  # Districts 1, 5, 10

    # Assemble components (order matters for S+E+I+R=N invariant)
    model.components = [
        SEIR.Susceptible(model),                                          # Count S
        SEIR.Exposed(model, expdurdist, infdurdist),                     # E->I transitions
        SEIR.Infectious(model, infdurdist),                               # I->R transitions
        SEIR.Recovered(model),                                            # Count R
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid), # Births
        MortalityByCDR(model, mortalityrates=deathrate_array),           # Deaths
        PatchImportation(model, infdurdist, endemic_patches,             # Endemic seeding
                         period=30, count=2),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),   # Transmission (S->E)
    ]

    # Initialize ages from pyramid BEFORE model.run()
    initialize_ages(model, pyramid)

    total_pop = scenario.population.sum()
    R0 = PARAMS.beta * PARAMS.inf_mean
    print(f"Running 10-year spatial SEIR simulation...")
    print(f"  Patches: {nnodes}, Population/patch: {POP_PER_PATCH:,}")
    print(f"  Total population: {total_pop:,}")
    print(f"  R0 = beta x D_inf = {PARAMS.beta:.4f} x {PARAMS.inf_mean} = {R0:.1f}")
    print(f"  Latent period: ~{PARAMS.exp_shape * PARAMS.exp_scale:.0f}d (gamma)")
    print(f"  Infectious period: ~{PARAMS.inf_mean}d (normal)")
    print(f"  Initial per patch: S={int(POP_PER_PATCH - INITIAL_IMMUNE_FRAC*POP_PER_PATCH) - N_SEED_INFECTIOUS:,}"
          f", I={N_SEED_INFECTIOUS}, R={int(INITIAL_IMMUNE_FRAC*POP_PER_PATCH):,}")
    print(f"  Seasonal forcing: monsoon peak +/-30%")
    print(f"  Gravity coupling: k={PARAMS.gravity_k}, b={PARAMS.gravity_b}, c={PARAMS.gravity_c}")
    print(f"  Importation: 2 infections in 3 patches every 30d")
    print(f"  Ticks: {nticks}")

    model.run("Polio SEIR 10-patch")

    return model, scenario, beta_season_365


# ============================================================================
# Summary Statistics
# ============================================================================

def print_summary(model, scenario):
    """Print key summary statistics."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    R0 = PARAMS.beta * PARAMS.inf_mean
    incidence = model.nodes.newly_infected
    nyears = nticks // 365

    print("\n" + "=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)

    print(f"\n{'Patch':<14} {'Avg Inf/yr':<14} {'Total Inf':<14} {'S/N(final)':<10}")
    print("-" * 52)

    total_infections = 0
    for i in range(nnodes):
        annual = np.array([incidence[y*365:(y+1)*365, i].sum() for y in range(nyears)])
        mean_ann = annual.mean()
        total = annual.sum()
        total_infections += total

        final_tick = nticks - 1
        S_f = float(model.nodes.S[final_tick, i])
        N_f = float(model.nodes.S[final_tick, i] + model.nodes.E[final_tick, i] +
                     model.nodes.I[final_tick, i] + model.nodes.R[final_tick, i])
        sn = S_f / max(N_f, 1)
        print(f"{names[i]:<14} {mean_ann:<14.0f} {total:<14.0f} {sn:<10.3f}")

    total_pop = scenario.population.sum()
    print(f"\n  Total infections over {nyears} years: {total_infections:,.0f}")
    print(f"  Avg annual infections: {total_infections / nyears:,.0f}")
    print(f"  Annual rate per million: {total_infections / nyears / total_pop * 1e6:.0f}")

    # Population trajectory
    print(f"\nPopulation trajectory:")
    for y in range(0, nyears + 1, 2):
        t = min(y * 365, nticks - 1)
        pop_t = int(model.nodes.S[t].sum() + model.nodes.E[t].sum() +
                     model.nodes.I[t].sum() + model.nodes.R[t].sum())
        print(f"  Year {y:>2}: {pop_t:>10,}")
    print(f"  Agent capacity: {model.people.capacity:,}, Active: {model.people.count:,}")

    # Compartment integrity check
    print(f"\nCompartment integrity:")
    for tick in [0, nticks // 2, nticks - 1]:
        S = model.nodes.S[tick].sum()
        E = model.nodes.E[tick].sum()
        I = model.nodes.I[tick].sum()
        R = model.nodes.R[tick].sum()
        print(f"  Tick {tick:>5}: S={S:>8} E={E:>5} I={I:>5} R={R:>8} N={S+E+I+R:>8}")

    # R_eff in second half of simulation
    sample = np.arange(nticks // 2, nticks, 30)
    S_all = model.nodes.S[sample].sum(axis=1).astype(np.float64)
    N_all = (model.nodes.S[sample].sum(axis=1) + model.nodes.E[sample].sum(axis=1) +
             model.nodes.I[sample].sum(axis=1) + model.nodes.R[sample].sum(axis=1)).astype(np.float64)
    reff = R0 * S_all / np.maximum(N_all, 1)
    print(f"\n  R_eff (years 5-10): mean={reff.mean():.2f}, "
          f"range=[{reff.min():.2f}, {reff.max():.2f}]")


# ============================================================================
# Diagnostic Plots
# ============================================================================

def plot_results(model, scenario, beta_season_365):
    """Generate 4-panel diagnostic plot."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    R0 = PARAMS.beta * PARAMS.inf_mean
    incidence = model.nodes.newly_infected
    nyears = nticks // 365

    cmap = plt.cm.tab10
    colors = [cmap(i) for i in range(nnodes)]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    # --- (A) Seasonal forcing ---
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, beta_season_365, "b-", lw=2)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
    ax.axvspan(182, 304, alpha=0.12, color="blue", label="Monsoon (Jul-Oct)")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Seasonal Multiplier")
    ax.set_title("(A) Monsoon Seasonal Forcing")
    ax.legend()

    # --- (B) Weekly incidence ---
    ax = axes[0, 1]
    nweeks = nticks // 7
    for i in range(nnodes):
        weekly = incidence[:nweeks*7, i].reshape(nweeks, 7).sum(axis=1)
        weeks_x = np.arange(nweeks) / 52.0
        ax.plot(weeks_x, weekly, color=colors[i], alpha=0.7, lw=0.8, label=names[i])
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly Infections")
    ax.set_title("(B) Weekly Incidence by Patch")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (C) Susceptible fraction ---
    ax = axes[1, 0]
    sample30 = np.arange(0, nticks, 30)
    time30 = sample30 / 365.0
    for i in range(nnodes):
        S = model.nodes.S[sample30, i].astype(np.float64)
        E = model.nodes.E[sample30, i].astype(np.float64)
        I = model.nodes.I[sample30, i].astype(np.float64)
        R = model.nodes.R[sample30, i].astype(np.float64)
        N = np.maximum(S + E + I + R, 1.0)
        ax.plot(time30, S / N, color=colors[i], alpha=0.7, lw=0.8, label=names[i])
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.5, label=f"S* = 1/R0 = {1/R0:.3f}")
    ax.set_xlabel("Year")
    ax.set_ylabel("Susceptible Fraction (S/N)")
    ax.set_title("(C) Susceptible Fraction")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (D) R_eff over time ---
    ax = axes[1, 1]
    for i in range(nnodes):
        S = model.nodes.S[sample30, i].astype(np.float64)
        E = model.nodes.E[sample30, i].astype(np.float64)
        I = model.nodes.I[sample30, i].astype(np.float64)
        R = model.nodes.R[sample30, i].astype(np.float64)
        N = np.maximum(S + E + I + R, 1.0)
        ax.plot(time30, R0 * S / N, color=colors[i], alpha=0.7, lw=0.8)
    ax.axhline(1.0, color="black", ls="-", lw=2, alpha=0.4, label="R_eff = 1")
    ax.set_xlabel("Year")
    ax.set_ylabel("R_eff = R0 x S/N")
    ax.set_title("(D) Effective Reproduction Number")
    ax.set_ylim(0, 4)
    ax.legend(fontsize=8)

    plt.suptitle(
        f"10-Patch Polio SEIR Model (R0={R0:.0f}, {POP_PER_PATCH/1000:.0f}K/patch, "
        f"{nyears}yr, monsoon forcing)",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()

    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "polio_seir_basic_10patch.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nDiagnostic plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    model, scenario, beta_season_365 = run_model()
    print_summary(model, scenario)
    plot_results(model, scenario, beta_season_365)
