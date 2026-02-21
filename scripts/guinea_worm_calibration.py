#!/usr/bin/env python3
"""
Guinea Worm Calibration Framework — Chad Dual-Host SEIS Model

Calibrates the LASER dual-host SEIS model against Carter Center
Guinea Worm Eradication Program surveillance data (2018-2023).

Calibration targets:
  - Annual human cases:    17, 48, 12, 14, 13, 6
  - Annual dog infections: 1040, 1935, 1507, 830, 688, 526
  - Dog-to-human ratio:    ~50:1 to 100:1
  - Declining trend reflecting intervention scale-up
  - Spatial concentration: cases cluster in 3-4 of 8 districts

Parameter space (5 dimensions):
  - beta_human:             U(0.0005, 0.005)
  - beta_dog:               U(0.005,  0.05)
  - cross_species_coupling: U(0.05,   0.5)
  - seasonal_amplitude:     U(0.3,    0.8)
  - containment_efficacy:   U(0.5,    0.9)

Multi-objective fitness:
  1. Log-scale MSE of annual human cases        (weight 2.0)
  2. Log-scale MSE of annual dog infections     (weight 2.0)
  3. Dog-to-human ratio penalty                 (weight 1.0)
  4. Declining trend penalty (human + dog)      (weight 1.5 each)
  5. Spatial concentration penalty              (weight 1.0 each)

Usage:
    /opt/anaconda3/bin/python3 scripts/guinea_worm_calibration.py
    /opt/anaconda3/bin/python3 scripts/guinea_worm_calibration.py --n-trials 128
"""

import sys
import argparse
import json
import time
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
import laser.core.random as lcr
from laser.generic import SEIR, Model
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance as haversine_distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from guinea_worm_components import (
    SEISRecovery, PatchImportation, build_chad_dry_season,
    CaseContainment, DogTethering,
    IntervenedWaterTransmission, IntervenedDualHostTransmission,
)


# ============================================================================
# Carter Center Surveillance Data (Chad, 2018-2023)
# ============================================================================

OBSERVED = pd.DataFrame({
    "year": [2018, 2019, 2020, 2021, 2022, 2023],
    "human_cases": [17, 48, 12, 14, 13, 6],
    "dog_infections": [1040, 1935, 1507, 830, 688, 526],
})

# Patches where cases concentrate (indices into DISTRICTS)
ENDEMIC_PATCHES = [0, 1, 2, 3]  # Chari-Baguirmi, Moyen-Chari, Salamat, Mandoul
N_DATA_YEARS = len(OBSERVED)


# ============================================================================
# District Configuration — 8 endemic districts in southern Chad
# ============================================================================

DISTRICTS = pd.DataFrame([
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
# Simulation Constants (fixed, not calibrated)
# ============================================================================

BURNIN_YEARS = 4
DATA_YEARS = 6
BURNIN_TICKS = BURNIN_YEARS * 365
DATA_TICKS = DATA_YEARS * 365
TOTAL_TICKS = BURNIN_TICKS + DATA_TICKS

# Guinea worm disease biology
EXP_SHAPE = 12       # gamma shape for pre-patent period
EXP_SCALE = 30.0     # gamma scale → mean = 360 days
INF_MEAN = 21        # worm emergence duration (days)
INF_SIGMA = 5

# Spatial coupling
GRAVITY_B = 0.3      # destination population exponent
GRAVITY_C = 2.0      # distance decay exponent

# Vital dynamics (per-1000/year — required by BirthsByCBR/MortalityByCDR)
HUMAN_CBR = 44.0
HUMAN_CDR = 13.0
DOG_CBR = 80.0
DOG_CDR = 100.0


# ============================================================================
# Parameter Space
# ============================================================================

PARAM_BOUNDS = {
    "beta_human":             (0.0005, 0.005),
    "beta_dog":               (0.005,  0.05),
    "cross_species_coupling": (0.05,   0.5),
    "seasonal_amplitude":     (0.3,    0.8),
    "containment_efficacy":   (0.5,    0.9),
}


# ============================================================================
# Section 1: Scoring Functions
# ============================================================================

def log_case_mse(simulated, observed):
    """Log-scale mean squared error for annual case counts.

    Uses log(x + 1) to handle zeros and give appropriate weight to
    small counts. A factor-of-2 error at 10 cases contributes the
    same as a factor-of-2 error at 1000 cases.

    Args:
        simulated: array of simulated annual counts (length N_DATA_YEARS)
        observed: array of observed annual counts (length N_DATA_YEARS)

    Returns:
        float: mean squared error on log(x+1) scale
    """
    log_sim = np.log(np.asarray(simulated, dtype=np.float64) + 1)
    log_obs = np.log(np.asarray(observed, dtype=np.float64) + 1)
    return float(np.mean((log_sim - log_obs) ** 2))


def ratio_penalty(human_annual, dog_annual, target_low=50, target_high=100):
    """Penalty for dog-to-human ratio outside target range.

    Carter Center data shows ~50:1 to ~100:1 dog-to-human ratio in Chad.
    No penalty when ratio is within range; quadratic log-scale penalty
    outside.

    Args:
        human_annual: array of annual human case counts
        dog_annual: array of annual dog infection counts
        target_low: lower bound of acceptable ratio (default 50)
        target_high: upper bound of acceptable ratio (default 100)

    Returns:
        float: penalty score (0.0 if ratio is in range)
    """
    total_human = max(float(np.sum(human_annual)), 1.0)
    total_dog = float(np.sum(dog_annual))
    ratio = total_dog / total_human

    if target_low <= ratio <= target_high:
        return 0.0

    log_ratio = np.log(ratio + 1)
    log_low = np.log(target_low + 1)
    log_high = np.log(target_high + 1)

    if ratio < target_low:
        return float((log_low - log_ratio) ** 2)
    return float((log_ratio - log_high) ** 2)


def trend_penalty(series):
    """Penalty for failing to show a declining trend.

    Fits a linear regression on log-counts over years. Only penalizes
    positive slopes (increasing trend). A declining series scores 0.

    Args:
        series: array of annual counts (should be declining)

    Returns:
        float: slope^2 if positive, 0.0 if declining or flat
    """
    n = len(series)
    if n < 2:
        return 0.0

    x = np.arange(n, dtype=np.float64)
    y = np.log(np.asarray(series, dtype=np.float64) + 1)
    x_mean, y_mean = x.mean(), y.mean()
    denom = np.sum((x - x_mean) ** 2)
    if denom < 1e-10:
        return 0.0
    slope = np.sum((x - x_mean) * (y - y_mean)) / denom
    return float(max(slope, 0.0) ** 2)


def spatial_concentration_score(annual_by_patch, endemic_indices,
                                threshold=0.80):
    """Penalty if endemic patches don't contain >= threshold of cases.

    Carter Center data shows cases cluster in 3-4 of 8 districts.
    Penalizes uniform spatial distribution.

    Args:
        annual_by_patch: 2D array (years x patches) of annual counts
        endemic_indices: list of patch indices expected to hold most cases
        threshold: minimum fraction in endemic patches (default 0.80)

    Returns:
        float: penalty (0.0 if spatial concentration is sufficient)
    """
    total = annual_by_patch.sum()
    if total == 0:
        return 10.0  # heavy penalty for no cases at all

    endemic_total = annual_by_patch[:, endemic_indices].sum()
    frac = endemic_total / total

    if frac >= threshold:
        return 0.0
    return float((threshold - frac) ** 2)


def combined_fitness(sim_human_annual, sim_dog_annual,
                     sim_human_by_patch, sim_dog_by_patch):
    """Multi-objective fitness combining all scoring components.

    Weights are set so each component contributes meaningfully when
    the model is moderately wrong. Lower total score is better.

    Args:
        sim_human_annual: (6,) simulated annual human cases
        sim_dog_annual: (6,) simulated annual dog infections
        sim_human_by_patch: (6, 8) human cases by year x patch
        sim_dog_by_patch: (6, 8) dog infections by year x patch

    Returns:
        total_loss: float, combined weighted loss
        scores: dict of individual component scores
    """
    obs_human = OBSERVED["human_cases"].values
    obs_dog = OBSERVED["dog_infections"].values

    scores = {
        "human_log_mse":    log_case_mse(sim_human_annual, obs_human),
        "dog_log_mse":      log_case_mse(sim_dog_annual, obs_dog),
        "ratio":            ratio_penalty(sim_human_annual, sim_dog_annual),
        "human_trend":      trend_penalty(sim_human_annual),
        "dog_trend":        trend_penalty(sim_dog_annual),
        "human_spatial":    spatial_concentration_score(
                                sim_human_by_patch, ENDEMIC_PATCHES),
        "dog_spatial":      spatial_concentration_score(
                                sim_dog_by_patch, ENDEMIC_PATCHES),
    }

    weights = {
        "human_log_mse": 2.0,   # primary target
        "dog_log_mse":   2.0,   # primary target
        "ratio":         1.0,   # species balance
        "human_trend":   1.5,   # intervention effect
        "dog_trend":     1.5,   # intervention effect
        "human_spatial": 1.0,   # geographic pattern
        "dog_spatial":   1.0,   # geographic pattern
    }

    total = sum(weights[k] * scores[k] for k in scores)
    return total, scores


# ============================================================================
# Section 2: Model Infrastructure
# ============================================================================

def build_scenario(pop_col, seed_patches, seed_counts):
    """Build GeoDataFrame scenario for one species."""
    geometry = [Point(lon, lat) for lat, lon
                in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(
        DISTRICTS.copy(), geometry=geometry, crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))
    scenario["population"] = scenario[pop_col].astype(np.uint32)
    scenario["E"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["I"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["R"] = np.zeros(len(scenario), dtype=np.uint32)

    for idx, n_inf in zip(seed_patches, seed_counts):
        scenario.loc[idx, "I"] = n_inf

    scenario["S"] = (scenario["population"] - scenario["E"]
                     - scenario["I"] - scenario["R"]).astype(np.uint32)

    assert (scenario.S + scenario.E + scenario.I + scenario.R
            == scenario.population).all(), \
        "Initial S+E+I+R must equal population in every patch"
    return scenario


def build_distance_matrix(scenario):
    """Pairwise haversine distance matrix (km)."""
    nnodes = len(scenario)
    lats, lons = np.array(scenario.lat), np.array(scenario.lon)
    d = np.zeros((nnodes, nnodes))
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                d[i, j] = haversine_distance(lats[i], lons[i],
                                             lats[j], lons[j])
    return d


def build_gravity_network(scenario, dist_matrix, gravity_k):
    """Build gravity-model migration network."""
    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(pops, dist_matrix, 1, 0, GRAVITY_B, GRAVITY_C)
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * gravity_k
    network = row_normalizer(network, 0.15)

    assert network.sum() > 0, "Network is all zeros — check gravity_k"
    assert network.sum(axis=1).max() < 0.3, \
        f"Max row sum = {network.sum(axis=1).max():.3f} — too high"
    return network


def build_age_pyramid(mean_life):
    """Simple exponential age pyramid."""
    ages = np.arange(100) if mean_life > 20 else np.arange(20)
    stable_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_dist = np.maximum(stable_dist, 1)
    return AliasedDistribution(stable_dist), stable_dist


def initialize_ages(model, pyramid_tuple):
    """Set initial agent DOBs from age pyramid."""
    count = model.people.count
    max_years = len(pyramid_tuple[1])
    ages_years = pyramid_tuple[0].sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, max_years - 1)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def extract_data_period_annual(model, burnin_ticks):
    """Extract annual new infections from the data period only.

    Returns:
        (N_DATA_YEARS, nnodes) array of annual infections per patch
    """
    nnodes = model.nodes.count
    annual = np.zeros((N_DATA_YEARS, nnodes), dtype=np.int64)
    for yr in range(N_DATA_YEARS):
        t0 = burnin_ticks + yr * 365
        t1 = burnin_ticks + (yr + 1) * 365
        annual[yr] = model.nodes.newly_infected[t0:t1, :].sum(axis=0)
    return annual


# ============================================================================
# Section 3: Pre-computed Context and Trial Runner
# ============================================================================

class CalibrationContext:
    """Pre-computes data shared across all calibration trials.

    Distance matrix, age pyramids, duration distributions, and the
    base seasonal profile are parameter-independent and only need
    to be built once.
    """

    def __init__(self):
        # Seed enough initial infections for burn-in equilibrium
        seed_patches = [0, 1, 2]

        self.human_scenario = build_scenario(
            "human_pop", seed_patches, [5, 5, 3])
        self.dog_scenario = build_scenario(
            "dog_pop", seed_patches, [30, 30, 20])

        self.dist_matrix = build_distance_matrix(self.human_scenario)
        self.human_pyramid = build_age_pyramid(55.0)
        self.dog_pyramid = build_age_pyramid(8.0)

        self.expdurdist = dists.gamma(shape=EXP_SHAPE, scale=EXP_SCALE)
        self.infdurdist = dists.normal(loc=INF_MEAN, scale=INF_SIGMA)

        # Cache base 365-day seasonal profile (normalized to mean=1)
        _, self.base_season_365 = build_chad_dry_season(365, 1)

        self.nnodes = len(DISTRICTS)


def run_trial(params, ctx, seed=42):
    """Run one calibration trial with the given parameters.

    Builds the dual-host SEIS model, runs a 4-year burn-in (no
    interventions) followed by a 6-year data period (containment
    ramps up linearly), and extracts annual case counts.

    Args:
        params: dict with keys matching PARAM_BOUNDS
        ctx: CalibrationContext with pre-computed shared data
        seed: random seed for reproducibility

    Returns:
        dict with keys:
            human_annual:   (6,) total human cases per year
            dog_annual:     (6,) total dog infections per year
            human_by_patch: (6, 8) human cases by year x patch
            dog_by_patch:   (6, 8) dog infections by year x patch
    """
    nnodes = ctx.nnodes

    # --- Unpack calibration parameters ---
    beta_human = params["beta_human"]
    beta_dog = params["beta_dog"]
    coupling = params["cross_species_coupling"]
    amplitude = params["seasonal_amplitude"]
    containment_eff = params["containment_efficacy"]

    # Asymmetric cross-species coupling (dogs dominate water contamination)
    dog_to_human = coupling
    human_to_dog = coupling * 0.33

    # --- Parameterized seasonality ---
    # Scale the Chad dry-season deviations by amplitude
    scaled_season = 1.0 + amplitude * (ctx.base_season_365 - 1.0)
    scaled_season /= scaled_season.mean()
    assert abs(scaled_season.mean() - 1.0) < 0.01
    season_tiled = np.tile(scaled_season, TOTAL_TICKS // 365 + 1)[:TOTAL_TICKS]

    # --- PropertySets ---
    human_ps = PropertySet({
        "prng_seed": seed, "nticks": TOTAL_TICKS,
        "beta": beta_human, "cbr": HUMAN_CBR,
        "gravity_k": 0.005, "gravity_b": GRAVITY_B, "gravity_c": GRAVITY_C,
        "capacity_safety_factor": 1.5,
    })
    dog_ps = PropertySet({
        "prng_seed": seed, "nticks": TOTAL_TICKS,
        "beta": beta_dog, "cbr": DOG_CBR,
        "gravity_k": 0.003, "gravity_b": GRAVITY_B, "gravity_c": GRAVITY_C,
        "capacity_safety_factor": 1.5,
    })

    # --- Birth/death rate arrays (per-1000/year) ---
    human_br = np.full((TOTAL_TICKS, nnodes), HUMAN_CBR, dtype=np.float32)
    human_dr = np.full((TOTAL_TICKS, nnodes), HUMAN_CDR, dtype=np.float32)
    dog_br = np.full((TOTAL_TICKS, nnodes), DOG_CBR, dtype=np.float32)
    dog_dr = np.full((TOTAL_TICKS, nnodes), DOG_CDR, dtype=np.float32)

    # --- Build LASER models ---
    human_model = Model(ctx.human_scenario.copy(), human_ps,
                        birthrates=human_br)
    dog_model = Model(ctx.dog_scenario.copy(), dog_ps,
                      birthrates=dog_br)

    # --- Gravity networks ---
    human_model.network = build_gravity_network(
        ctx.human_scenario, ctx.dist_matrix, 0.005)
    dog_model.network = build_gravity_network(
        ctx.dog_scenario, ctx.dist_matrix, 0.003)

    # --- Initialize agent ages ---
    initialize_ages(human_model, ctx.human_pyramid)
    initialize_ages(dog_model, ctx.dog_pyramid)

    # --- Intervention components (detection rate set dynamically) ---
    # CaseContainment.detection_rate = how many cases get detected
    # The per-case FOI reduction is fixed at 0.90
    h_containment = CaseContainment(
        human_model, detection_rate=0.0, containment_efficacy=0.90)
    d_tethering = DogTethering(
        dog_model, village_coverage=0.0, tether_efficacy=0.85)

    # --- Assemble components ---
    # Order: state propagation → SEIS recovery → interventions →
    #        vital dynamics → transmission
    human_model.components = [
        SEIR.Susceptible(human_model),
        SEIR.Exposed(human_model, ctx.expdurdist, ctx.infdurdist),
        SEIR.Infectious(human_model, ctx.infdurdist),
        SEIR.Recovered(human_model),
        SEISRecovery(human_model),
        h_containment,
        BirthsByCBR(human_model, birthrates=human_br,
                    pyramid=ctx.human_pyramid[0]),
        MortalityByCDR(human_model, mortalityrates=human_dr),
        IntervenedWaterTransmission(
            human_model, ctx.expdurdist,
            seasonality_values=season_tiled,
            containment_efficacy=0.90, filter_efficacy=0.0),
    ]

    dog_model.components = [
        SEIR.Susceptible(dog_model),
        SEIR.Exposed(dog_model, ctx.expdurdist, ctx.infdurdist),
        SEIR.Infectious(dog_model, ctx.infdurdist),
        SEIR.Recovered(dog_model),
        SEISRecovery(dog_model),
        d_tethering,
        BirthsByCBR(dog_model, birthrates=dog_br,
                    pyramid=ctx.dog_pyramid[0]),
        MortalityByCDR(dog_model, mortalityrates=dog_dr),
        IntervenedWaterTransmission(
            dog_model, ctx.expdurdist,
            seasonality_values=season_tiled,
            tether_efficacy=0.85, filter_efficacy=0.0),
    ]

    # --- Cross-species coupling ---
    cross_host = IntervenedDualHostTransmission(
        human_model, dog_model, ctx.expdurdist,
        dog_to_human=dog_to_human,
        human_to_dog=human_to_dog,
        seasonality_values=season_tiled,
        containment_efficacy=0.90, tether_efficacy=0.85,
        filter_efficacy=0.0,
    )

    # --- Importation (prevents stochastic extinction) ---
    endemic = [0, 1, 2]
    h_import = PatchImportation(
        human_model, ctx.expdurdist, endemic,
        period=180, count=1, end_tick=TOTAL_TICKS)
    d_import = PatchImportation(
        dog_model, ctx.expdurdist, endemic,
        period=60, count=2, end_tick=TOTAL_TICKS)

    # --- Run coupled simulation ---
    lcr.seed(seed)
    np.random.seed(seed)

    for tick in range(TOTAL_TICKS):
        # Time-varying containment: off during burn-in, ramps up during data
        if tick < BURNIN_TICKS:
            h_containment.detection_rate = 0.0
            d_tethering.village_coverage = 0.0
        else:
            progress = (tick - BURNIN_TICKS) / DATA_TICKS
            # Human case detection ramps from 20% to containment_efficacy
            h_containment.detection_rate = (
                0.2 + progress * (containment_eff - 0.2))
            # Dog tethering ramps from 15% to 70% of containment_efficacy
            d_tethering.village_coverage = (
                0.15 + progress * (containment_eff * 0.7 - 0.15))

        # Initialize flows: copy S[t]→S[t+1], etc.
        human_model._initialize_flows(tick)
        dog_model._initialize_flows(tick)

        # Step each model's components
        for c in human_model.components:
            c.step(tick)
        for c in dog_model.components:
            c.step(tick)

        # Cross-species transmission via shared water sources
        cross_host.step(tick)

        # Background importation
        h_import.step(tick)
        d_import.step(tick)

    # --- Extract data-period results ---
    h_by_patch = extract_data_period_annual(human_model, BURNIN_TICKS)
    d_by_patch = extract_data_period_annual(dog_model, BURNIN_TICKS)

    return {
        "human_annual": h_by_patch.sum(axis=1),
        "dog_annual": d_by_patch.sum(axis=1),
        "human_by_patch": h_by_patch,
        "dog_by_patch": d_by_patch,
    }


# ============================================================================
# Section 4: Calibration Loop
# ============================================================================

def sample_params(n, method="sobol", seed=0):
    """Generate n parameter sets from the calibration space.

    Args:
        n: number of parameter sets
        method: "sobol" (quasi-random, space-filling) or "random"
        seed: random seed

    Returns:
        list of dicts, each with keys matching PARAM_BOUNDS
    """
    names = list(PARAM_BOUNDS.keys())
    lowers = np.array([PARAM_BOUNDS[k][0] for k in names])
    uppers = np.array([PARAM_BOUNDS[k][1] for k in names])

    if method == "sobol":
        try:
            from scipy.stats.qmc import Sobol
            sampler = Sobol(d=len(names), scramble=True, seed=seed)
            m = int(np.ceil(np.log2(max(n, 2))))
            raw = sampler.random_base2(m)[:n]
        except ImportError:
            print("scipy.stats.qmc not available, falling back to random")
            rng = np.random.default_rng(seed)
            raw = rng.uniform(size=(n, len(names)))
    else:
        rng = np.random.default_rng(seed)
        raw = rng.uniform(size=(n, len(names)))

    samples = []
    for i in range(n):
        p = {}
        for j, name in enumerate(names):
            p[name] = float(lowers[j] + raw[i, j] * (uppers[j] - lowers[j]))
        samples.append(p)
    return samples


def calibrate(n_trials=64, seed=0):
    """Run calibration sweep and return results sorted by fitness.

    Uses Sobol quasi-random sampling to explore the 5D parameter
    space, evaluates each parameter set by running the full dual-host
    model, and ranks results by combined multi-objective fitness.

    Args:
        n_trials: number of parameter sets to evaluate (default 64)
        seed: random seed for reproducibility

    Returns:
        list of result dicts, sorted by loss (best first)
    """
    print("Initializing calibration context...")
    ctx = CalibrationContext()

    print(f"Generating {n_trials} parameter sets (Sobol)...")
    param_sets = sample_params(n_trials, method="sobol", seed=seed)

    results = []
    best_loss = float("inf")

    header = (f"{'Trial':>5} {'Loss':>10} {'H_MSE':>8} {'D_MSE':>8} "
              f"{'Ratio':>8} {'H_Trnd':>8} {'D_Trnd':>8} "
              f"{'H_Sptl':>8} {'D_Sptl':>8}  {'Time':>6}")
    print(f"\n{header}")
    print("-" * len(header))

    for i, params in enumerate(param_sets):
        t0 = time.time()

        try:
            outputs = run_trial(params, ctx, seed=seed + i)
            loss, scores = combined_fitness(
                outputs["human_annual"], outputs["dog_annual"],
                outputs["human_by_patch"], outputs["dog_by_patch"],
            )
            elapsed = time.time() - t0

            is_best = loss < best_loss
            if is_best:
                best_loss = loss

            result = {
                "trial": i,
                "params": params,
                "loss": loss,
                "scores": scores,
                "human_annual": outputs["human_annual"].tolist(),
                "dog_annual": outputs["dog_annual"].tolist(),
                "elapsed": elapsed,
            }
            results.append(result)

            status = " **" if is_best else ""
            print(f"{i:>5} {loss:>10.4f} "
                  f"{scores['human_log_mse']:>8.4f} "
                  f"{scores['dog_log_mse']:>8.4f} "
                  f"{scores['ratio']:>8.4f} "
                  f"{scores['human_trend']:>8.4f} "
                  f"{scores['dog_trend']:>8.4f} "
                  f"{scores['human_spatial']:>8.4f} "
                  f"{scores['dog_spatial']:>8.4f}"
                  f"  {elapsed:>5.1f}s{status}")

        except Exception as e:
            elapsed = time.time() - t0
            print(f"{i:>5} {'ERROR':>10}  {str(e)[:60]}")
            results.append({
                "trial": i, "params": params, "loss": float("inf"),
                "scores": {}, "elapsed": elapsed,
                "human_annual": [], "dog_annual": [],
            })

    # Sort by loss
    results.sort(key=lambda r: r["loss"])

    # --- Print summary ---
    print(f"\n{'=' * 70}")
    print(f"CALIBRATION COMPLETE — {n_trials} trials")
    print(f"{'=' * 70}")

    valid = [r for r in results if r["loss"] < float("inf")]
    if valid:
        best = valid[0]
        print(f"\nBest trial #{best['trial']} (loss = {best['loss']:.4f}):")
        print(f"\n  Parameters:")
        for k, v in best["params"].items():
            lo, hi = PARAM_BOUNDS[k]
            print(f"    {k:>25s} = {v:.6f}  "
                  f"[{lo:.4f}, {hi:.4f}]")

        print(f"\n  Annual comparison (simulated vs observed):")
        print(f"  {'Year':>6} {'Sim_H':>8} {'Obs_H':>8} "
              f"{'Sim_D':>8} {'Obs_D':>8}")
        for yr_idx in range(N_DATA_YEARS):
            yr = OBSERVED["year"].iloc[yr_idx]
            sh = best["human_annual"][yr_idx]
            oh = OBSERVED["human_cases"].iloc[yr_idx]
            sd = best["dog_annual"][yr_idx]
            od = OBSERVED["dog_infections"].iloc[yr_idx]
            print(f"  {yr:>6} {sh:>8} {oh:>8} {sd:>8} {od:>8}")

        h_total = sum(best["human_annual"])
        d_total = sum(best["dog_annual"])
        ratio = d_total / max(h_total, 1)
        print(f"\n  Dog:Human ratio = {ratio:.1f}:1 (target: 50-100:1)")

        # Top 5 parameter sets
        print(f"\n  Top 5 parameter sets:")
        print(f"  {'#':>3} {'Loss':>10} {'beta_h':>10} {'beta_d':>10} "
              f"{'coupling':>10} {'season':>10} {'contain':>10}")
        for r in valid[:5]:
            p = r["params"]
            print(f"  {r['trial']:>3} {r['loss']:>10.4f} "
                  f"{p['beta_human']:>10.5f} "
                  f"{p['beta_dog']:>10.5f} "
                  f"{p['cross_species_coupling']:>10.4f} "
                  f"{p['seasonal_amplitude']:>10.4f} "
                  f"{p['containment_efficacy']:>10.4f}")

    return results


# ============================================================================
# Section 5: Diagnostics and Plotting
# ============================================================================

def plot_calibration_results(results, outdir=None):
    """Generate diagnostic plots from calibration results.

    Produces a 4-panel figure:
        1. Loss convergence (trial order + cumulative best)
        2. Best-fit human cases vs observed
        3. Best-fit dog infections vs observed
        4. Parameter sensitivity (normalized params vs loss)

    Also saves top-20 results as JSON.
    """
    if outdir is None:
        outdir = Path(__file__).parent.parent / "eval" / "guinea-worm" / "outputs"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    valid = [r for r in results if r["loss"] < float("inf")]
    if not valid:
        print("No valid results to plot.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Guinea Worm Calibration — Chad 2018-2023",
                 fontsize=14, fontweight="bold")

    # 1. Loss convergence
    ax = axes[0, 0]
    trial_order = [r["trial"] for r in valid]
    losses = [r["loss"] for r in valid]
    # Re-sort by trial order for convergence plot
    order = np.argsort(trial_order)
    ordered_losses = [losses[i] for i in order]
    cummin = np.minimum.accumulate(ordered_losses)
    ax.plot(ordered_losses, ".", alpha=0.4, color="steelblue",
            label="Trial loss")
    ax.plot(cummin, "r-", linewidth=2, label="Best so far")
    ax.set_xlabel("Trial (Sobol order)")
    ax.set_ylabel("Combined Loss")
    ax.set_title("Calibration Convergence")
    ax.legend()
    ax.set_yscale("log")

    # 2. Best fit: human cases
    ax = axes[0, 1]
    best = valid[0]
    years = OBSERVED["year"].values
    ax.bar(years - 0.2, OBSERVED["human_cases"].values, 0.35,
           label="Observed", color="steelblue", alpha=0.8)
    ax.bar(years + 0.2, best["human_annual"], 0.35,
           label="Simulated", color="coral", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual Human Cases")
    ax.set_title("Human Cases: Best Fit vs Observed")
    ax.legend()

    # 3. Best fit: dog infections
    ax = axes[1, 0]
    ax.bar(years - 0.2, OBSERVED["dog_infections"].values, 0.35,
           label="Observed", color="steelblue", alpha=0.8)
    ax.bar(years + 0.2, best["dog_annual"], 0.35,
           label="Simulated", color="coral", alpha=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual Dog Infections")
    ax.set_title("Dog Infections: Best Fit vs Observed")
    ax.legend()

    # 4. Parameter sensitivity
    ax = axes[1, 1]
    param_names = list(PARAM_BOUNDS.keys())
    short_names = ["beta_h", "beta_d", "coupling", "season", "contain"]
    losses_arr = np.array([r["loss"] for r in valid])
    for j, (pname, sname) in enumerate(zip(param_names, short_names)):
        vals = np.array([r["params"][pname] for r in valid])
        lo, hi = PARAM_BOUNDS[pname]
        normalized = (vals - lo) / (hi - lo)
        ax.scatter(normalized, losses_arr, alpha=0.3, s=12, label=sname)
    ax.set_xlabel("Normalized Parameter Value (0=lower, 1=upper)")
    ax.set_ylabel("Loss")
    ax.set_title("Parameter Sensitivity")
    ax.set_yscale("log")
    ax.legend(fontsize=8, ncol=2)

    plt.tight_layout()
    fig_path = outdir / "calibration_results.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to {fig_path}")
    plt.close(fig)

    # Save top results as JSON
    serializable = []
    for r in valid[:20]:
        serializable.append({
            "trial": r["trial"],
            "loss": r["loss"],
            "scores": r["scores"],
            "params": r["params"],
            "human_annual": r["human_annual"],
            "dog_annual": r["dog_annual"],
        })
    results_path = outdir / "calibration_results.json"
    with open(results_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"Results saved to {results_path}")


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Guinea Worm Calibration — Chad Dual-Host SEIS Model")
    parser.add_argument(
        "--n-trials", type=int, default=64,
        help="Number of calibration trials (default 64, use power of 2 "
             "for Sobol)")
    parser.add_argument(
        "--seed", type=int, default=0,
        help="Random seed (default 0)")
    args = parser.parse_args()

    results = calibrate(n_trials=args.n_trials, seed=args.seed)
    plot_calibration_results(results)
