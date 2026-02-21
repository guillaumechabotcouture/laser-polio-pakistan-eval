#!/usr/bin/env python3
"""
Calibration framework for Pakistan polio SEIRV model.

Samples a 3D parameter space (beta, gravity_k, seasonal_amplitude),
runs LASER SEIRV simulations, and scores against AFP surveillance data
using a CCS-based fitness metric (proportion of zero-incidence weeks
per district vs. population size).

The SEIRV model tracks vaccine-derived immunity (V) separately from
natural infection immunity (R). OPV mucosal immunity wanes (~3 years),
while natural immunity is permanent.

Fitness function:
    1. Derives expected weekly infection rates per district from MMWR AFP
       surveillance data using the 1:200 paralysis-to-infection ratio.
    2. Computes the proportion of zero-infection weeks for each district
       (model) and the Poisson-expected proportion (observed).
    3. Scores as SSE between model and expected zero-week proportions,
       plus a log-RMSE AFP match term for province-level annual totals.
    4. Ranks all simulations by combined fitness (lower = better).

Usage:
    /opt/anaconda3/bin/python3 scripts/polio_calibration.py
    /opt/anaconda3/bin/python3 scripts/polio_calibration.py --n-samples 64 --seed 7
"""

import sys
import argparse
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer
from laser.core.migration import distance as haversine_distance

# Add scripts dir to path for custom_components and main model
sys.path.insert(0, str(Path(__file__).parent))
from custom_components import (VaccinatedCompartment, PerPatchVaccinationSEIRV,
                               PatchImportation, VACCINATED)
from polio_seir_10patch import (
    DISTRICTS, BURNIN_YEARS, PARALYSIS_RATIO,
    OBSERVED_WPV1, DISTRICT_TO_PROVINCE,
    build_scenario, build_pakistan_age_pyramid, initialize_ages,
)

# ============================================================================
# Parameter Space
# ============================================================================

PARAM_RANGES = {
    "beta":                {"low": 0.15, "high": 0.30, "scale": "linear"},
    "gravity_k":           {"low": -3.0, "high": -1.0, "scale": "log10"},
    "seasonal_amplitude":  {"low": 1.0,  "high": 1.5,  "scale": "linear"},
}

# Fixed parameters (not calibrated)
FIXED = {
    "exp_shape": 3,
    "exp_scale": 1.0,
    "inf_mean": 28,
    "inf_sigma": 3,
    "cbr": 29.0,
    "cdr": 7.0,
    "gravity_b": 0.5,
    "gravity_c": 1.5,
    "capacity_safety_factor": 2.0,
    "nticks": 20 * 365,
    "sia_period": 180,
    "import_period": 30,
    "import_count": 3,
}


# ============================================================================
# Observed Data / Calibration Targets
# ============================================================================

def derive_observed_targets(scenario):
    """Derive per-district calibration targets from MMWR AFP surveillance.

    Uses the 1:200 paralysis-to-infection ratio to convert observed AFP
    counts into expected infection rates, then computes expected proportion
    of zero-infection weeks per district (Poisson assumption).

    Returns:
        dict with district_afp_yr, expected_prop_zero, prov_afp_yr
    """
    names = list(scenario.name)
    nnodes = len(scenario)

    # Province-level average AFP/year (2020-2025)
    prov_afp_yr = {}
    for prov in ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]:
        prov_afp_yr[prov] = float(OBSERVED_WPV1[prov].mean())

    # Model population per province (sum of modeled districts)
    prov_model_pop = {}
    for i in range(nnodes):
        prov = DISTRICT_TO_PROVINCE[names[i]]
        prov_model_pop[prov] = prov_model_pop.get(prov, 0) + int(scenario.population.iloc[i])

    # Distribute AFP to districts proportional to population within province
    district_afp_yr = np.zeros(nnodes)
    for i in range(nnodes):
        prov = DISTRICT_TO_PROVINCE[names[i]]
        pop_frac = scenario.population.iloc[i] / prov_model_pop[prov]
        district_afp_yr[i] = prov_afp_yr[prov] * pop_frac

    # Convert AFP → infections → weekly rate → expected P(zero) via Poisson
    district_infections_wk = district_afp_yr / PARALYSIS_RATIO / 52.0
    expected_prop_zero = np.exp(-district_infections_wk)

    return {
        "district_afp_yr": district_afp_yr,
        "expected_prop_zero": expected_prop_zero,
        "district_infections_wk": district_infections_wk,
        "prov_afp_yr": prov_afp_yr,
    }


# ============================================================================
# Scoring Functions
# ============================================================================

def logistic(x, x0, k):
    """Logistic S-curve for CCS fit: transitions from 1→0 at x0."""
    return 1.0 / (1.0 + np.exp(k * (x - x0)))


def compute_zero_week_proportions(weekly_incidence):
    """Proportion of weeks with zero infections per district.

    Args:
        weekly_incidence: 2D array (weeks, districts)
    Returns:
        1D array of shape (districts,)
    """
    return np.mean(weekly_incidence == 0, axis=0)


def score_ccs_pattern(model_prop_zero, target_prop_zero, log_populations):
    """Score how well the model's zero-week profile matches observed targets.

    Two components:
        1. Direct SSE between model and target prop_zero per district.
        2. Quality of logistic CCS fit (does the model produce a clean
           population-size-dependent extinction gradient?).

    Returns:
        float (lower is better)
    """
    # Component 1: Direct district-level comparison
    direct_sse = float(np.sum((model_prop_zero - target_prop_zero) ** 2))

    # Component 2: Logistic CCS fit quality
    try:
        p0 = [np.median(log_populations), 3.0]
        bounds = ([2.0, 0.01], [8.0, 20.0])
        popt, _ = curve_fit(
            logistic, log_populations, model_prop_zero,
            p0=p0, bounds=bounds, maxfev=5000,
        )
        fitted = logistic(log_populations, *popt)
        fit_residual = float(np.sum((model_prop_zero - fitted) ** 2))

        # Penalise if CCS midpoint is outside plausible range for polio
        # (log10 pop 4.5–6.5, i.e. ~30K–3M given vaccination)
        ccs_x0 = popt[0]
        ccs_penalty = 0.0
        if ccs_x0 < 4.5 or ccs_x0 > 6.5:
            ccs_penalty = (max(4.5 - ccs_x0, ccs_x0 - 6.5)) ** 2

        logistic_score = fit_residual + 0.1 * ccs_penalty
    except (RuntimeError, ValueError):
        logistic_score = 10.0

    return direct_sse + 0.3 * logistic_score


def score_afp_match(model_annual_infections, scenario, prov_afp_targets):
    """Province-level annual AFP comparison via log-RMSE.

    Converts model infections to expected AFP (÷200) and compares to
    MMWR observed annual averages on a log scale.

    Returns:
        float RMSLE (lower is better)
    """
    names = list(scenario.name)
    nnodes = len(scenario)

    model_prov_infections = {}
    for i in range(nnodes):
        prov = DISTRICT_TO_PROVINCE[names[i]]
        yearly_mean = float(model_annual_infections[:, i].mean())
        model_prov_infections[prov] = model_prov_infections.get(prov, 0.0) + yearly_mean

    log_diffs_sq = []
    for prov in ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]:
        model_afp = model_prov_infections.get(prov, 0.0) * PARALYSIS_RATIO
        obs_afp = prov_afp_targets[prov]
        diff = np.log(model_afp + 1) - np.log(obs_afp + 1)
        log_diffs_sq.append(diff ** 2)

    return float(np.sqrt(np.mean(log_diffs_sq)))


def combined_fitness(weekly_incidence, annual_infections, scenario, targets,
                     ccs_weight=0.6, afp_weight=0.4):
    """Compute combined fitness: CCS zero-week pattern + AFP match.

    Args:
        weekly_incidence: (weeks, districts) infection counts
        annual_infections: (years, districts) infection counts
        scenario: GeoDataFrame
        targets: from derive_observed_targets()
        ccs_weight / afp_weight: relative weighting

    Returns:
        dict with total_score, ccs_score, afp_score, prop_zero
    """
    prop_zero = compute_zero_week_proportions(weekly_incidence)
    log_pops = np.log10(np.array(scenario.population, dtype=np.float64))

    ccs = score_ccs_pattern(prop_zero, targets["expected_prop_zero"], log_pops)
    afp = score_afp_match(annual_infections, scenario, targets["prov_afp_yr"])

    total = ccs_weight * ccs + afp_weight * afp
    return {
        "total_score": total,
        "ccs_score": ccs,
        "afp_score": afp,
        "prop_zero": prop_zero,
    }


# ============================================================================
# Model Runner
# ============================================================================

def _precompute_shared(scenario):
    """Precompute expensive data that is constant across all parameter sets."""
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

    pyramid, _ = build_pakistan_age_pyramid()

    name_to_idx = {n: i for i, n in enumerate(scenario.name)}
    endemic_names = [
        "Quetta", "Peshawar", "Karachi", "Hyderabad",
        "Bannu", "N_Waziristan", "S_Waziristan", "DI_Khan",
    ]
    endemic_patches = [name_to_idx[n] for n in endemic_names if n in name_to_idx]

    return {
        "dist_matrix": dist_matrix,
        "pyramid": pyramid,
        "endemic_patches": endemic_patches,
    }


def run_single_simulation(params, scenario, shared, seed=42):
    """Build and run one LASER simulation with given calibration parameters.

    Args:
        params: dict with beta, gravity_k, seasonal_amplitude
        scenario: GeoDataFrame from build_scenario()
        shared: from _precompute_shared()
        seed: random seed

    Returns:
        dict with weekly_incidence (weeks, districts),
                  annual_infections (years, districts)
    """
    np.random.seed(seed)

    nnodes = len(scenario)
    nticks = FIXED["nticks"]
    burnin_ticks = BURNIN_YEARS * 365
    nyears = nticks // 365
    analysis_years = nyears - BURNIN_YEARS

    # --- Seasonal forcing ---
    # seasonal_amplitude in [1.0, 1.5]: cosine amplitude = seasonal_amplitude - 1
    # At 1.0 → flat; at 1.5 → profile ranges [0.5, 1.5] with mean 1.0
    amplitude = params["seasonal_amplitude"] - 1.0
    days = np.arange(365)
    peak_day = 245  # Early September (monsoon)
    season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
    season_tiled = np.tile(season_365, nyears + 1)[:nticks]
    seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

    # --- Duration distributions ---
    expdurdist = dists.gamma(shape=FIXED["exp_shape"], scale=FIXED["exp_scale"])
    infdurdist = dists.normal(loc=FIXED["inf_mean"], scale=FIXED["inf_sigma"])

    # --- Vital dynamics arrays ---
    birthrate_array = np.full((nticks, nnodes), FIXED["cbr"], dtype=np.float32)
    mortality_array = np.full((nticks, nnodes), FIXED["cdr"], dtype=np.float32)

    # --- LASER PropertySet ---
    laser_params = PropertySet({
        "prng_seed": seed,
        "nticks": nticks,
        "beta": params["beta"],
        "exp_shape": FIXED["exp_shape"],
        "exp_scale": FIXED["exp_scale"],
        "inf_mean": FIXED["inf_mean"],
        "inf_sigma": FIXED["inf_sigma"],
        "cbr": FIXED["cbr"],
        "gravity_k": params["gravity_k"],
        "gravity_b": FIXED["gravity_b"],
        "gravity_c": FIXED["gravity_c"],
        "capacity_safety_factor": FIXED["capacity_safety_factor"],
    })

    # --- Build model with V compartment ---
    model = Model(scenario, laser_params, birthrates=birthrate_array,
                  additional_states=["V"])

    # --- Gravity migration network ---
    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(
        pops, shared["dist_matrix"],
        1, 0, FIXED["gravity_b"], FIXED["gravity_c"],
    )
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * params["gravity_k"]
    network = row_normalizer(network, 0.2)
    model.network = network

    # --- OPV waning distribution (mean ~3 years) ---
    vax_waning_dist = dists.gamma(shape=3, scale=365.0)

    # --- Vaccination (correlated missedness) ---
    unreachable_frac = np.array(scenario.unreachable_frac, dtype=np.float32)

    # --- Assemble components (order matters) ---
    susceptible = SEIR.Susceptible(model)
    exposed = SEIR.Exposed(model, expdurdist, infdurdist)
    infectious = SEIR.Infectious(model, infdurdist)
    recovered = SEIR.Recovered(model)

    vax_compartment = VaccinatedCompartment(model, vax_waning_dist, wandurmin=365)
    vaccination = PerPatchVaccinationSEIRV(
        model, vax_compartment, unreachable_frac,
        ri_coverage=0.80, sia_coverage=0.90,
        ri_age=42, sia_period=FIXED["sia_period"], sia_max_age=1825,
    )

    model.components = [
        susceptible,
        exposed,
        infectious,
        recovered,
        vax_compartment,
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=shared["pyramid"]),
        MortalityByCDR(model, mortalityrates=mortality_array, mappings=[
            (SEIR.State.SUSCEPTIBLE.value, "S"),
            (SEIR.State.EXPOSED.value, "E"),
            (SEIR.State.INFECTIOUS.value, "I"),
            (SEIR.State.RECOVERED.value, "R"),
            (VACCINATED, "V"),
        ]),
        vaccination,
        PatchImportation(
            model, infdurdist, shared["endemic_patches"],
            period=FIXED["import_period"], count=FIXED["import_count"],
        ),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    ]

    initialize_ages(model, shared["pyramid"])
    model.run("Calibration")

    # --- Extract post-burn-in results ---
    incidence = model.nodes.newly_infected

    # Weekly incidence
    post_inc = incidence[burnin_ticks:]
    nweeks = post_inc.shape[0] // 7
    weekly = post_inc[: nweeks * 7].reshape(nweeks, 7, nnodes).sum(axis=1)

    # Annual incidence
    annual = np.zeros((analysis_years, nnodes))
    for y in range(analysis_years):
        yr = BURNIN_YEARS + y
        annual[y] = incidence[yr * 365 : (yr + 1) * 365].sum(axis=0)

    return {"weekly_incidence": weekly, "annual_infections": annual}


# ============================================================================
# Parameter Sampling
# ============================================================================

def sample_parameters(n_samples, seed=42):
    """Latin Hypercube Sampling over the 3D parameter space.

    Handles log-uniform distribution for gravity_k (10^U(-3, -1)).

    Returns:
        list of dicts with beta, gravity_k, seasonal_amplitude
    """
    try:
        from scipy.stats.qmc import LatinHypercube
        sampler = LatinHypercube(d=3, seed=seed)
        unit = sampler.random(n_samples)
    except ImportError:
        # Stratified-random fallback
        rng = np.random.default_rng(seed)
        unit = np.zeros((n_samples, 3))
        for j in range(3):
            perm = rng.permutation(n_samples)
            for i in range(n_samples):
                unit[i, j] = (perm[i] + rng.random()) / n_samples

    param_sets = []
    for i in range(n_samples):
        u = unit[i]
        beta = 0.15 + u[0] * 0.15                       # U(0.15, 0.30)
        gravity_k = 10 ** (-3.0 + u[1] * 2.0)           # 10^U(-3, -1)
        seasonal_amplitude = 1.0 + u[2] * 0.5            # U(1.0, 1.5)
        param_sets.append({
            "beta": beta,
            "gravity_k": gravity_k,
            "seasonal_amplitude": seasonal_amplitude,
        })
    return param_sets


# ============================================================================
# Calibration Loop
# ============================================================================

def run_calibration(n_samples=32, seed=42):
    """Run full calibration sweep.

    Returns:
        results_df: DataFrame sorted by total_score (best first)
        scenario: GeoDataFrame
        targets: observed calibration targets
    """
    print("=" * 70)
    print("POLIO SEIRV CALIBRATION FRAMEWORK")
    print("=" * 70)

    scenario = build_scenario()
    shared = _precompute_shared(scenario)
    targets = derive_observed_targets(scenario)
    names = list(scenario.name)

    print(f"\nCalibration targets (MMWR AFP, 1:{int(1/PARALYSIS_RATIO)} ratio):")
    for i, name in enumerate(names):
        print(
            f"  {name:<14} pop={scenario.population.iloc[i]:>9,}  "
            f"AFP/yr={targets['district_afp_yr'][i]:5.2f}  "
            f"inf/wk={targets['district_infections_wk'][i]:6.1f}  "
            f"E[prop_zero]={targets['expected_prop_zero'][i]:.4f}"
        )

    param_sets = sample_parameters(n_samples, seed=seed)

    print(f"\n{n_samples} parameter sets (Latin Hypercube Sampling):")
    print(f"  beta               : U(0.15, 0.30)")
    print(f"  gravity_k          : 10^U(-3, -1)  [{1e-3:.4f}, {1e-1:.1f}]")
    print(f"  seasonal_amplitude : U(1.0, 1.5)")
    print()

    results = []
    for i, params in enumerate(param_sets):
        t0 = time.time()

        try:
            out = run_single_simulation(params, scenario, shared, seed=seed + i)

            fit = combined_fitness(
                out["weekly_incidence"], out["annual_infections"],
                scenario, targets,
            )

            # Province-level AFP summary for diagnostics
            prov_inf = {}
            for j in range(len(names)):
                prov = DISTRICT_TO_PROVINCE[names[j]]
                prov_inf[prov] = prov_inf.get(prov, 0.0) + float(
                    out["annual_infections"][:, j].mean()
                )
            total_afp = sum(v * PARALYSIS_RATIO for v in prov_inf.values())

            elapsed = time.time() - t0

            row = {
                "sample_id": i,
                "beta": params["beta"],
                "gravity_k": params["gravity_k"],
                "seasonal_amplitude": params["seasonal_amplitude"],
                "total_score": fit["total_score"],
                "ccs_score": fit["ccs_score"],
                "afp_score": fit["afp_score"],
                "total_afp_yr": total_afp,
                "elapsed_sec": elapsed,
            }
            for j, name in enumerate(names):
                row[f"pz_{name}"] = fit["prop_zero"][j]
            results.append(row)

            print(
                f"  [{i+1:>3}/{n_samples}] "
                f"beta={params['beta']:.4f}  "
                f"gk={params['gravity_k']:.5f}  "
                f"amp={params['seasonal_amplitude']:.2f}  "
                f"=> score={fit['total_score']:.4f} "
                f"(CCS={fit['ccs_score']:.3f} AFP={fit['afp_score']:.3f}) "
                f"AFP/yr={total_afp:5.1f}  "
                f"[{elapsed:.0f}s]"
            )

        except Exception as e:
            elapsed = time.time() - t0
            print(
                f"  [{i+1:>3}/{n_samples}] "
                f"beta={params['beta']:.4f}  "
                f"gk={params['gravity_k']:.5f}  "
                f"amp={params['seasonal_amplitude']:.2f}  "
                f"=> FAILED: {e}  [{elapsed:.0f}s]"
            )
            results.append({
                "sample_id": i,
                "beta": params["beta"],
                "gravity_k": params["gravity_k"],
                "seasonal_amplitude": params["seasonal_amplitude"],
                "total_score": float("inf"),
                "ccs_score": float("inf"),
                "afp_score": float("inf"),
                "total_afp_yr": 0.0,
                "elapsed_sec": elapsed,
            })

    df = pd.DataFrame(results)
    df = df.sort_values("total_score").reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df, scenario, targets


# ============================================================================
# Results Reporting
# ============================================================================

def print_results(results_df, scenario, targets):
    """Print ranked calibration results and diagnostics."""
    names = list(scenario.name)

    print("\n" + "=" * 70)
    print("CALIBRATION RESULTS  (ranked by combined fitness)")
    print("=" * 70)

    top_n = min(10, len(results_df))
    print(f"\nTop {top_n} parameter sets:")
    print(
        f"{'Rank':<5} {'beta':<8} {'gravity_k':<11} {'ampl':<7} "
        f"{'Score':<9} {'CCS':<8} {'AFP':<8} {'AFP/yr':<8}"
    )
    print("-" * 64)
    for _, row in results_df.head(top_n).iterrows():
        print(
            f"{int(row['rank']):<5} "
            f"{row['beta']:<8.4f} "
            f"{row['gravity_k']:<11.6f} "
            f"{row['seasonal_amplitude']:<7.3f} "
            f"{row['total_score']:<9.4f} "
            f"{row['ccs_score']:<8.4f} "
            f"{row['afp_score']:<8.4f} "
            f"{row['total_afp_yr']:<8.1f}"
        )

    best = results_df.iloc[0]
    obs_total = float(
        OBSERVED_WPV1[["Balochistan", "KP", "Sindh", "Punjab", "ICT"]]
        .sum(axis=1).mean()
    )

    print(f"\nBest parameters:")
    print(f"  beta               = {best['beta']:.4f}")
    print(f"  gravity_k          = {best['gravity_k']:.6f}")
    print(f"  seasonal_amplitude = {best['seasonal_amplitude']:.3f}")
    print(f"  Total fitness      = {best['total_score']:.4f}")
    print(f"  Expected AFP/yr    = {best['total_afp_yr']:.1f}  "
          f"(observed: {obs_total:.1f})")

    # Zero-week comparison for best run
    pz_cols = [c for c in results_df.columns if c.startswith("pz_")]
    if pz_cols:
        print(f"\nProportion of zero-infection weeks (best run vs target):")
        print(f"  {'District':<14} {'Pop':>9} {'Model':>7} {'Target':>7} {'Diff':>8}")
        print(f"  {'-' * 45}")
        for j, name in enumerate(names):
            col = f"pz_{name}"
            if col in results_df.columns:
                m = best[col]
                t = targets["expected_prop_zero"][j]
                pop = scenario.population.iloc[j]
                print(
                    f"  {name:<14} {pop:>9,} {m:>7.4f} {t:>7.4f} {m - t:>+8.4f}"
                )

    # Save
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "calibration_results.csv"
    results_df.to_csv(outpath, index=False, float_format="%.6f")
    print(f"\nResults saved to {outpath}")


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate polio SEIRV model against AFP surveillance data"
    )
    parser.add_argument(
        "--n-samples", type=int, default=32,
        help="Number of LHS parameter sets (default: 32)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for sampling and simulations (default: 42)",
    )
    args = parser.parse_args()

    results_df, scenario, targets = run_calibration(
        n_samples=args.n_samples, seed=args.seed,
    )
    print_results(results_df, scenario, targets)
