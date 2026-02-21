#!/usr/bin/env python3
"""
Chad Guinea Worm Dual-Host SEIS Spatial Model with Interventions

Complete 10-year simulation with time-varying prevention interventions
and elimination tracking. Models Dracunculus medinensis transmission
across 8 endemic districts in southern Chad.

Architecture:
  Two parallel LASER Models (human, dog) coupled via shared water sources.
  Each model uses IntervenedWaterTransmission (gravity-coupled FOI with
  containment, tethering, ABATE, and filter modifiers). Cross-species
  coupling via IntervenedDualHostTransmission.

Disease dynamics:
  - SEIS: no lasting immunity (recovered → susceptible)
  - Humans: R0≈1.5, 365-day pre-patent period, 21-day worm emergence
  - Dogs:   R0≈2.5, 365-day pre-patent period, 21-day worm emergence
  - Cross-species: dogs dominant reservoir (~95%), asymmetric coupling

Interventions (all prevention-based, no vaccine exists):
  - ABATE larvicide:      30% → 90% over years (kills copepods in water)
  - Filter distribution:  40% → 90% (removes copepods from drinking water)
  - Case containment:     70% → 97% (prevents infectious humans from water)
  - Dog tethering:        40% → 85% (prevents infectious dogs from water)

Usage:
    /opt/anaconda3/bin/python3 scripts/guinea_worm_chad_interventions.py
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
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance as haversine_distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from guinea_worm_components import (
    SEISRecovery, build_chad_dry_season,
    ABATELarvicide, WaterFilterDistribution, CaseContainment,
    DogTethering, IntervenedWaterTransmission, IntervenedDualHostTransmission,
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
    ("Logone Oriental",  250000,    20000,   8.62,    16.08),
    ("Tandjile",         200000,    15000,   9.65,    16.72),
    ("Mayo-Kebbi Est",   150000,    10000,   9.37,    15.35),
    ("Lac",              100000,     5000,   13.30,   14.10),
], columns=["name", "human_pop", "dog_pop", "lat", "lon"])

# ============================================================================
# Disease Parameters
# ============================================================================

NTICKS = 10 * 365  # 10-year simulation
NNODES = len(DISTRICTS)

# Guinea worm natural history (shared across species)
EXP_SHAPE = 12        # gamma shape for pre-patent period
EXP_SCALE = 30.0      # gamma scale → mean = 12 × 30 = 360 days ≈ 365
INF_MEAN = 21          # worm emergence period (days)
INF_SIGMA = 5          # emergence variability

HUMAN_PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": NTICKS,
    "exp_shape": EXP_SHAPE,
    "exp_scale": EXP_SCALE,
    "inf_mean": INF_MEAN,
    "inf_sigma": INF_SIGMA,
    "beta": 1.5 / INF_MEAN,       # R0 ≈ 1.5 → beta ≈ 0.0714
    "cbr": 44.0,                   # Chad CBR (per 1000/yr)
    "gravity_k": 0.005,
    "gravity_b": 0.3,              # destination population exponent
    "gravity_c": 2.0,              # distance decay (water sources are local)
    "capacity_safety_factor": 2.0,
})

DOG_PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": NTICKS,
    "exp_shape": EXP_SHAPE,
    "exp_scale": EXP_SCALE,
    "inf_mean": INF_MEAN,
    "inf_sigma": INF_SIGMA,
    "beta": 2.5 / INF_MEAN,       # R0 ≈ 2.5 → beta ≈ 0.119
    "cbr": 80.0,                   # dog CBR (per 1000/yr)
    "gravity_k": 0.003,            # dogs less mobile than humans
    "gravity_b": 0.3,
    "gravity_c": 2.0,
    "capacity_safety_factor": 1.5,
})

# Cross-species coupling via shared water sources
# Dogs are dominant reservoir (~95% of infections in Chad)
DOG_TO_HUMAN_COUPLING = 0.15   # dog prevalence → human FOI (attenuated by filters)
HUMAN_TO_DOG_COUPLING = 0.05   # human prevalence → dog FOI (very weak)

# Vital dynamics
HUMAN_CDR = 13.0   # per 1000/yr (Chad national)
DOG_CDR = 100.0    # per 1000/yr (~10yr average lifespan)

# Intervention efficacies (fixed — coverage ramps up, efficacy stays constant)
CONTAINMENT_EFFICACY = 0.90   # water contamination reduction for contained cases
TETHER_EFFICACY = 0.85        # water contamination reduction for tethered dogs
FILTER_EFFICACY = 0.95        # copepod ingestion reduction per filter

# Endemic patches (indices 0, 1, 2 = Chari-Baguirmi, Moyen-Chari, Salamat)
ENDEMIC_PATCHES = np.array([0, 1, 2])


# ============================================================================
# Intervention Schedule
# ============================================================================

def get_intervention_coverages(year):
    """Return intervention coverage for simulation year (0-indexed).

    ABATE larvicide:     30% year 0, +10%/year to max 90%
    Filter distribution: 40% year 0, +10%/year to max 90%
    Case containment:    70% year 0, +3%/year (detection rate)
    Dog tethering:       40% year 0, +5%/year
    """
    return {
        "abate": min(0.30 + 0.10 * year, 0.90),
        "filter": min(0.40 + 0.10 * year, 0.90),
        "containment": min(0.70 + 0.03 * year, 1.0),
        "tethering": min(0.40 + 0.05 * year, 1.0),
    }


# ============================================================================
# Helper Functions
# ============================================================================

def build_scenario(pop_col, seed_patches, seed_counts):
    """Build GeoDataFrame scenario for one species."""
    geometry = [Point(lon, lat) for lat, lon in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(DISTRICTS.copy(), geometry=geometry, crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))
    scenario["population"] = scenario[pop_col].astype(np.uint32)

    scenario["E"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["I"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["R"] = np.zeros(len(scenario), dtype=np.uint32)

    for patch_idx, n_infected in zip(seed_patches, seed_counts):
        scenario.loc[patch_idx, "I"] = n_infected

    scenario["S"] = (scenario["population"] - scenario["E"]
                     - scenario["I"] - scenario["R"]).astype(np.uint32)

    assert (scenario.S + scenario.E + scenario.I + scenario.R
            == scenario.population).all(), \
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
    network = row_normalizer(network, 0.15)

    assert network.sum() > 0, "Network is all zeros — check gravity_k"
    assert network.sum(axis=1).max() < 0.3, \
        f"Max row sum = {network.sum(axis=1).max():.3f} — too high"

    return network


def build_age_pyramid(mean_life):
    """Build simple exponential age pyramid."""
    max_age = 100 if mean_life > 20 else 20
    ages = np.arange(max_age)
    stable_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_dist = np.maximum(stable_dist, 1)
    return AliasedDistribution(stable_dist), stable_dist


def initialize_ages(model, pyramid_dist, pyramid_arr):
    """Set initial agent DOBs from age pyramid."""
    count = model.people.count
    max_years = len(pyramid_arr)
    ages_years = pyramid_dist.sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, max_years - 1)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def update_filter_coverage(model, filter_component, new_rate):
    """Distribute additional filters to reach new target adoption rate.

    Finds agents without filters and distributes enough new filters
    so the overall population rate reaches new_rate.
    """
    count = model.people.count
    if count == 0:
        return

    current_filters = model.people.has_filter[:count].astype(bool)
    n_have = current_filters.sum()
    target_with = int(new_rate * count)
    n_needed = target_with - n_have

    if n_needed <= 0:
        return

    no_filter_idx = np.nonzero(~current_filters)[0]
    if len(no_filter_idx) == 0:
        return

    n_give = min(n_needed, len(no_filter_idx))
    chosen = np.random.choice(no_filter_idx, size=n_give, replace=False)
    model.people.has_filter[chosen] = 1

    # Update component adoption rate for future newborns
    filter_component.adoption_rate[:] = new_rate


# ============================================================================
# Main Model
# ============================================================================

def run_guinea_worm_model():
    """Build and run the dual-host guinea worm SEIS model with interventions."""

    # --- Seed initial infections ---
    # 5 infected humans in 3 most endemic districts
    human_seed_patches = [0, 1, 2]
    human_seed_counts = [2, 2, 1]

    # 200 infected dogs in 3 most endemic districts
    dog_seed_patches = [0, 1, 2]
    dog_seed_counts = [80, 70, 50]

    # --- Build scenarios ---
    human_scenario = build_scenario("human_pop", human_seed_patches,
                                    human_seed_counts)
    dog_scenario = build_scenario("dog_pop", dog_seed_patches, dog_seed_counts)

    # --- Distance matrix (shared geography) ---
    dist_matrix = build_distance_matrix(human_scenario)

    # --- Duration distributions (shared across species) ---
    expdurdist = dists.gamma(shape=EXP_SHAPE, scale=EXP_SCALE)
    infdurdist = dists.normal(loc=INF_MEAN, scale=INF_SIGMA)

    # --- Seasonality (shared climate) ---
    seasonality_map_h, season_365 = build_chad_dry_season(NTICKS, NNODES)
    seasonality_map_d, _ = build_chad_dry_season(NTICKS, NNODES)
    season_tiled = np.tile(season_365, NTICKS // 365 + 1)[:NTICKS]

    # --- Birth/death rates ---
    human_birthrates = np.full((NTICKS, NNODES), HUMAN_PARAMS.cbr,
                               dtype=np.float32)
    human_deathrates = np.full((NTICKS, NNODES), HUMAN_CDR, dtype=np.float32)
    dog_birthrates = np.full((NTICKS, NNODES), DOG_PARAMS.cbr,
                              dtype=np.float32)
    dog_deathrates = np.full((NTICKS, NNODES), DOG_CDR, dtype=np.float32)

    # --- Age pyramids ---
    h_pyr_dist, h_pyr_arr = build_age_pyramid(55.0)
    d_pyr_dist, d_pyr_arr = build_age_pyramid(8.0)

    # --- Initial intervention coverages (year 0) ---
    y0 = get_intervention_coverages(0)

    # ==================================================================
    # BUILD HUMAN MODEL
    # ==================================================================
    print("Building human model...")
    human_model = Model(human_scenario, HUMAN_PARAMS,
                        birthrates=human_birthrates)
    human_model.network = build_gravity_network(human_scenario, dist_matrix,
                                                HUMAN_PARAMS)

    # SEIR state components (must create Susceptible first — it initializes
    # people.state and people.nodeid, which intervention components need)
    h_susceptible = SEIR.Susceptible(human_model)
    h_exposed = SEIR.Exposed(human_model, expdurdist, infdurdist)
    h_infectious = SEIR.Infectious(human_model, infdurdist)
    h_recovered = SEIR.Recovered(human_model)
    h_seis = SEISRecovery(human_model)

    # Intervention components (created after Susceptible so nodeid exists)
    h_abate = ABATELarvicide(human_model, coverage_by_patch=y0["abate"])
    h_filters = WaterFilterDistribution(human_model,
                                        adoption_rate=y0["filter"],
                                        filter_efficacy=FILTER_EFFICACY)
    h_containment = CaseContainment(human_model,
                                    detection_rate=y0["containment"],
                                    containment_efficacy=CONTAINMENT_EFFICACY)

    # Within-species transmission with all intervention modifiers
    h_transmission = IntervenedWaterTransmission(
        human_model, expdurdist, seasonality_values=season_tiled,
        containment_efficacy=CONTAINMENT_EFFICACY,
        filter_efficacy=FILTER_EFFICACY,
    )

    human_model.components = [
        h_susceptible,
        h_exposed,
        h_infectious,
        h_recovered,
        h_seis,                              # R→S (no lasting immunity)
        BirthsByCBR(human_model, birthrates=human_birthrates,
                    pyramid=h_pyr_dist),
        MortalityByCDR(human_model, mortalityrates=human_deathrates),
        h_abate,                             # set abate_factor
        h_filters,                           # passive (on_birth assigns)
        h_containment,                       # containment decisions
        h_transmission,                      # gravity-coupled FOI
    ]

    initialize_ages(human_model, h_pyr_dist, h_pyr_arr)

    # ==================================================================
    # BUILD DOG MODEL
    # ==================================================================
    print("Building dog model...")
    dog_model = Model(dog_scenario, DOG_PARAMS, birthrates=dog_birthrates)
    dog_model.network = build_gravity_network(dog_scenario, dist_matrix,
                                              DOG_PARAMS)

    # SEIR state components (Susceptible first for nodeid)
    d_susceptible = SEIR.Susceptible(dog_model)
    d_exposed = SEIR.Exposed(dog_model, expdurdist, infdurdist)
    d_infectious = SEIR.Infectious(dog_model, infdurdist)
    d_recovered = SEIR.Recovered(dog_model)
    d_seis = SEISRecovery(dog_model)

    # Dog intervention components
    d_abate = ABATELarvicide(dog_model, coverage_by_patch=y0["abate"])
    d_tethering = DogTethering(dog_model,
                               village_coverage=y0["tethering"],
                               tether_efficacy=TETHER_EFFICACY)

    d_transmission = IntervenedWaterTransmission(
        dog_model, expdurdist, seasonality_values=season_tiled,
        tether_efficacy=TETHER_EFFICACY,
    )

    dog_model.components = [
        d_susceptible,
        d_exposed,
        d_infectious,
        d_recovered,
        d_seis,
        BirthsByCBR(dog_model, birthrates=dog_birthrates,
                    pyramid=d_pyr_dist),
        MortalityByCDR(dog_model, mortalityrates=dog_deathrates),
        d_abate,
        d_tethering,
        d_transmission,
    ]

    initialize_ages(dog_model, d_pyr_dist, d_pyr_arr)

    # ==================================================================
    # CROSS-SPECIES COUPLING
    # ==================================================================
    cross_host = IntervenedDualHostTransmission(
        human_model, dog_model, expdurdist,
        dog_to_human=DOG_TO_HUMAN_COUPLING,
        human_to_dog=HUMAN_TO_DOG_COUPLING,
        seasonality_values=season_tiled,
        containment_efficacy=CONTAINMENT_EFFICACY,
        tether_efficacy=TETHER_EFFICACY,
        filter_efficacy=FILTER_EFFICACY,
    )

    # ==================================================================
    # RUN 10-YEAR SIMULATION
    # ==================================================================
    print(f"\nRunning 10-year dual-host guinea worm SEIS simulation...")
    print(f"  Districts:    {NNODES}")
    print(f"  Human pop:    {human_scenario.population.sum():,}")
    print(f"  Dog pop:      {dog_scenario.population.sum():,}")
    print(f"  Human R0:     {HUMAN_PARAMS.beta * INF_MEAN:.1f}")
    print(f"  Dog R0:       {DOG_PARAMS.beta * INF_MEAN:.1f}")
    print(f"  Coupling:     dog→human={DOG_TO_HUMAN_COUPLING}, "
          f"human→dog={HUMAN_TO_DOG_COUPLING}")
    print(f"  Seasonality:  dry-season peak (Jan-Mar 1.8x, Jun-Sep 0.4x)")
    print(f"  Interventions: ABATE + filters + containment + tethering")
    print()

    import laser.core.random as lcr
    lcr.seed(HUMAN_PARAMS.prng_seed)

    for tick in range(NTICKS):
        # --- Yearly intervention coverage updates ---
        if tick > 0 and tick % 365 == 0:
            year = tick // 365
            cov = get_intervention_coverages(year)

            # ABATE (both models — same water sources)
            h_abate.coverage[:] = cov["abate"]
            d_abate.coverage[:] = cov["abate"]

            # Filter distribution (humans only)
            update_filter_coverage(human_model, h_filters, cov["filter"])

            # Case containment detection rate
            h_containment.detection_rate = cov["containment"]

            # Dog tethering coverage
            d_tethering.village_coverage = cov["tethering"]

            # Annual progress report
            h_new = human_model.nodes.newly_infected[
                (year - 1) * 365:tick, :].sum()
            d_new = dog_model.nodes.newly_infected[
                (year - 1) * 365:tick, :].sum()
            print(f"  Year {year:>2}: "
                  f"H_cases={h_new:>6,}  D_cases={d_new:>6,}  |  "
                  f"ABATE={cov['abate']:.0%}  Filt={cov['filter']:.0%}  "
                  f"Cont={cov['containment']:.0%}  Teth={cov['tethering']:.0%}")

        # --- Initialize flows: copy compartment counts tick → tick+1 ---
        human_model._initialize_flows(tick)
        dog_model._initialize_flows(tick)

        # --- Step all model components ---
        for component in human_model.components:
            component.step(tick)
        for component in dog_model.components:
            component.step(tick)

        # --- Cross-species transmission via shared water ---
        cross_host.step(tick)

        # --- Human importation: 1 traveler every 90 days for 3 years ---
        if tick > 0 and tick % 90 == 0 and tick < 3 * 365:
            patch_id = int(np.random.choice(ENDEMIC_PATCHES))
            count = human_model.people.count
            susc_in_patch = np.nonzero(
                (human_model.people.state[:count]
                 == SEIR.State.SUSCEPTIBLE.value) &
                (human_model.people.nodeid[:count] == patch_id)
            )[0]
            if len(susc_in_patch) > 0:
                agent = np.random.choice(susc_in_patch)
                human_model.people.state[agent] = SEIR.State.EXPOSED.value
                etimer_val = max(1, int(round(dists.sample_floats(
                    expdurdist, np.zeros(1, dtype=np.float32)
                )[0])))
                human_model.people.etimer[agent] = etimer_val
                human_model.nodes.S[tick + 1, patch_id] -= 1
                human_model.nodes.S[tick + 1, patch_id] = max(
                    0, human_model.nodes.S[tick + 1, patch_id])
                human_model.nodes.E[tick + 1, patch_id] += 1

    # Final year report
    final_year = NTICKS // 365
    h_new = human_model.nodes.newly_infected[
        (final_year - 1) * 365:NTICKS, :].sum()
    d_new = dog_model.nodes.newly_infected[
        (final_year - 1) * 365:NTICKS, :].sum()
    cov = get_intervention_coverages(final_year - 1)
    print(f"  Year {final_year:>2}: "
          f"H_cases={h_new:>6,}  D_cases={d_new:>6,}  |  "
          f"ABATE={cov['abate']:.0%}  Filt={cov['filter']:.0%}  "
          f"Cont={cov['containment']:.0%}  Teth={cov['tethering']:.0%}")
    print("\nSimulation complete.")

    return {
        "human_model": human_model,
        "dog_model": dog_model,
        "human_scenario": human_scenario,
        "dog_scenario": dog_scenario,
        "season_365": season_365,
    }


# ============================================================================
# Results Extraction
# ============================================================================

def extract_annual_results(model):
    """Extract annual new infection counts per patch."""
    nnodes = model.nodes.count
    nyears = NTICKS // 365
    newly_inf = model.nodes.newly_infected[:NTICKS, :]

    annual = np.zeros((nyears, nnodes), dtype=np.int64)
    for yr in range(nyears):
        t0 = yr * 365
        t1 = (yr + 1) * 365
        annual[yr] = newly_inf[t0:t1].sum(axis=0)

    return annual


def check_elimination(h_annual, d_annual):
    """Check if elimination (zero cases for 3 consecutive years) is achieved.

    Returns:
        eliminated: bool
        elimination_year: int or None (1-indexed year starting the window)
    """
    nyears = h_annual.shape[0]

    for yr in range(nyears - 2):
        h_window = h_annual[yr:yr + 3].sum()
        d_window = d_annual[yr:yr + 3].sum()
        if h_window == 0 and d_window == 0:
            return True, yr + 1  # 1-indexed

    return False, None


def print_results(results):
    """Print comprehensive annual results and elimination status."""
    h_annual = extract_annual_results(results["human_model"])
    d_annual = extract_annual_results(results["dog_model"])
    human_scenario = results["human_scenario"]
    nyears = NTICKS // 365

    print("\n" + "=" * 78)
    print("GUINEA WORM ELIMINATION MODEL — CHAD (10-YEAR PROJECTION)")
    print("=" * 78)

    # Intervention schedule
    print("\nIntervention Coverage Schedule:")
    print(f"{'Year':<6} {'ABATE':>8} {'Filters':>8} {'Contain':>8} "
          f"{'Tether':>8}")
    print("-" * 40)
    for yr in range(nyears):
        cov = get_intervention_coverages(yr)
        print(f"{yr + 1:<6} {cov['abate']:>7.0%} {cov['filter']:>7.0%} "
              f"{cov['containment']:>7.0%} {cov['tethering']:>7.0%}")

    # National totals
    print(f"\n{'Year':<6} {'Human Cases':>12} {'Dog Infections':>15} "
          f"{'Dog:Human':>12}")
    print("-" * 47)
    for yr in range(nyears):
        h_total = h_annual[yr].sum()
        d_total = d_annual[yr].sum()
        if h_total > 0:
            ratio = f"{d_total / h_total:.0f}:1"
        elif d_total > 0:
            ratio = "inf:1"
        else:
            ratio = "0:0"
        print(f"{yr + 1:<6} {h_total:>12,} {d_total:>15,} {ratio:>12}")

    # Per-district breakdown
    print(f"\n{'District':<20}", end="")
    for yr in range(nyears):
        print(f" {'Y' + str(yr + 1) + ' H':>6} {'Y' + str(yr + 1) + ' D':>6}",
              end="")
    print()
    print("-" * (20 + nyears * 13))

    for i in range(NNODES):
        name = DISTRICTS.iloc[i]["name"][:18]
        print(f"{name:<20}", end="")
        for yr in range(nyears):
            print(f" {h_annual[yr, i]:>6} {d_annual[yr, i]:>6}", end="")
        print()

    # Elimination check
    eliminated, elim_year = check_elimination(h_annual, d_annual)
    print("\n" + "=" * 78)
    if eliminated:
        print(f"ELIMINATION ACHIEVED: Zero cases in both species for "
              f"3 consecutive years starting year {elim_year}.")
    else:
        print("ELIMINATION NOT ACHIEVED within the 10-year simulation "
              "window.")
        total_by_year = h_annual.sum(axis=1) + d_annual.sum(axis=1)
        min_year = np.argmin(total_by_year)
        print(f"  Lowest transmission year: Year {min_year + 1} "
              f"(H={h_annual[min_year].sum():,}, "
              f"D={d_annual[min_year].sum():,})")
    print("=" * 78)

    return h_annual, d_annual


# ============================================================================
# Verification
# ============================================================================

def verify_model(results):
    """Run verification checks on the dual-host SEIS model."""
    human_model = results["human_model"]
    dog_model = results["dog_model"]
    checks = []

    for label, model in [("Human", human_model), ("Dog", dog_model)]:
        # SEIS: R should be near zero (SEISRecovery moves R→S each tick)
        R_final = model.nodes.R[NTICKS - 1, :].sum()
        checks.append((f"{label} SEIS (R≈0)", R_final < 50,
                        f"R={R_final}"))

        # Non-negativity
        S = model.nodes.S[:NTICKS, :]
        E = model.nodes.E[:NTICKS, :]
        I = model.nodes.I[:NTICKS, :]
        checks.append((f"{label} non-negativity",
                        not ((S < 0).any() or (E < 0).any()
                             or (I < 0).any()),
                        f"S<0:{(S < 0).any()}, E<0:{(E < 0).any()}, "
                        f"I<0:{(I < 0).any()}"))

        # Some infections occurred
        total_inf = model.nodes.newly_infected[:NTICKS, :].sum()
        checks.append((f"{label} infections>0", total_inf > 0,
                        f"Total: {total_inf:,}"))

        # Population growth (humans should grow with CBR>CDR)
        pop_start = (S[0] + E[0] + I[0] + model.nodes.R[0]).sum()
        pop_end = (S[-1] + E[-1] + I[-1] + model.nodes.R[-1]).sum()
        if label == "Human":
            checks.append((f"{label} pop growth", pop_end > pop_start,
                            f"{pop_start:,} → {pop_end:,}"))

        # Network connectivity
        checks.append((f"{label} network active", model.network.sum() > 0,
                        f"Sum: {model.network.sum():.4f}"))

    print("\n" + "=" * 60)
    print("VERIFICATION REPORT")
    print("=" * 60)
    all_pass = True
    for name, passed, detail in checks:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  [{status}] {name}: {detail}")

    print(f"\n{'All checks passed.' if all_pass else 'Some checks failed.'}")
    return all_pass


# ============================================================================
# Diagnostic Plots
# ============================================================================

def plot_results(results, h_annual, d_annual):
    """Generate diagnostic plots."""
    human_model = results["human_model"]
    dog_model = results["dog_model"]
    season_365 = results["season_365"]
    nyears = NTICKS // 365

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Chad Guinea Worm SEIS Model — 10-Year Intervention Scenario",
                 fontsize=14, fontweight="bold")

    # 1. Seasonal forcing profile
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, season_365, "b-", linewidth=2)
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    ax.fill_between(days, season_365, alpha=0.2)
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Transmission Multiplier")
    ax.set_title("Dry-Season Forcing")
    months = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
    ax.set_xticks([i * 365 // 12 for i in range(12)])
    ax.set_xticklabels(months)

    # 2. Weekly incidence (both species)
    ax = axes[0, 1]
    h_daily = human_model.nodes.newly_infected[:NTICKS, :].sum(axis=1)
    d_daily = dog_model.nodes.newly_infected[:NTICKS, :].sum(axis=1)
    weeks = np.arange(0, NTICKS - 6, 7)
    h_weekly = np.array([h_daily[w:w + 7].sum() for w in weeks])
    d_weekly = np.array([d_daily[w:w + 7].sum() for w in weeks])
    ax.plot(weeks / 365, h_weekly, "b-", label="Human", alpha=0.8)
    ax.plot(weeks / 365, d_weekly, "r-", label="Dog", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly New Infections")
    ax.set_title("Epidemic Curves")
    ax.legend()

    # 3. Annual cases bar chart
    ax = axes[0, 2]
    x = np.arange(1, nyears + 1)
    w = 0.35
    ax.bar(x - w / 2, h_annual.sum(axis=1), w, label="Human",
           color="steelblue")
    ax.bar(x + w / 2, d_annual.sum(axis=1), w, label="Dog",
           color="indianred")
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual Infections")
    ax.set_title("Annual Infections by Species")
    ax.legend()
    ax.set_xticks(x)

    # 4. Intervention coverage timeline
    ax = axes[1, 0]
    years = np.arange(nyears)
    covs = [get_intervention_coverages(yr) for yr in years]
    ax.plot(years + 1, [c["abate"] for c in covs], "g-o",
            label="ABATE", markersize=4)
    ax.plot(years + 1, [c["filter"] for c in covs], "b-s",
            label="Filters", markersize=4)
    ax.plot(years + 1, [c["containment"] for c in covs], "r-^",
            label="Containment", markersize=4)
    ax.plot(years + 1, [c["tethering"] for c in covs], "m-d",
            label="Tethering", markersize=4)
    ax.set_xlabel("Year")
    ax.set_ylabel("Coverage / Efficacy")
    ax.set_title("Intervention Scale-Up")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    # 5. Spatial heatmap (dog infections by district)
    ax = axes[1, 1]
    nweeks = NTICKS // 7
    d_by_patch = dog_model.nodes.newly_infected[:nweeks * 7, :].reshape(
        nweeks, 7, NNODES).sum(axis=1)
    im = ax.imshow(d_by_patch.T, aspect="auto", cmap="YlOrRd",
                   extent=[0, NTICKS / 365, NNODES - 0.5, -0.5])
    ax.set_xlabel("Year")
    ax.set_ylabel("District")
    ax.set_title("Dog Infections by District (Weekly)")
    ax.set_yticks(range(NNODES))
    names = [n[:12] for n in DISTRICTS.name]
    ax.set_yticklabels(names, fontsize=7)
    plt.colorbar(im, ax=ax, label="Cases/week")

    # 6. Cumulative infections
    ax = axes[1, 2]
    h_cum = np.cumsum(h_annual.sum(axis=1))
    d_cum = np.cumsum(d_annual.sum(axis=1))
    ax.plot(x, h_cum, "b-o", label="Human (cumul.)", markersize=4)
    ax.plot(x, d_cum, "r-o", label="Dog (cumul.)", markersize=4)
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative Infections")
    ax.set_title("Cumulative Infections")
    ax.legend()

    plt.tight_layout()
    outpath = Path(__file__).parent.parent / "eval" / "guinea-worm" / "outputs"
    outpath.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath / "guinea_worm_interventions.png", dpi=150,
                bbox_inches="tight")
    print(f"\nPlot saved to {outpath / 'guinea_worm_interventions.png'}")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    results = run_guinea_worm_model()
    h_annual, d_annual = print_results(results)
    verify_model(results)
    plot_results(results, h_annual, d_annual)
