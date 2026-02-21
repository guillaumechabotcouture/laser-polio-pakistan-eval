#!/usr/bin/env python3
"""
Calibration framework for Pakistan Polio SEIRV model.

Samples parameter space using Latin Hypercube Sampling, runs the LASER model
for each parameter set, and scores against AFP surveillance data using two
fitness metrics:

1. AFP Provincial Comparison (AFP loss):
   Model infections × 1/200 vs MMWR observed paralytic cases.
   Log-scale MSE across provinces.

2. Zero-AFP-Week Proportion (CCS loss):
   For each district, computes the expected proportion of weeks with zero
   AFP detections given simulated infection counts and the 1:200 paralysis-
   to-infection ratio. Compared to observed targets derived from MMWR annual
   data via Poisson approximation. This metric captures the relationship
   between population size and disease persistence — the hallmark of
   critical community size (CCS) analysis.

Parameter ranges:
    beta:                U(0.15, 0.30)       daily transmission rate
    gravity_k:           10^U(-3, -1)        spatial coupling strength
    seasonal_amplitude:  U(1.0, 1.5)         monsoon forcing half-amplitude

Usage:
    /opt/anaconda3/bin/python3 scripts/calibrate_polio.py --n-samples 32 --seed 42
"""

import sys
import argparse
import time
import gc
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add scripts dir for imports
sys.path.insert(0, str(Path(__file__).parent))

from polio_seir_10patch import (
    build_scenario,
    build_pakistan_age_pyramid,
    initialize_ages,
    OBSERVED_WPV1,
    PARALYSIS_RATIO,
    DISTRICT_TO_PROVINCE,
    BURNIN_YEARS,
)
from custom_components import (VaccinatedCompartment, PerPatchVaccinationSEIRV,
                               PatchImportation)

from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer
from laser.core.migration import distance as haversine_distance


# ============================================================================
# Parameter Space
# ============================================================================

PARAM_RANGES = {
    "beta":                {"low": 0.15, "high": 0.30},
    "gravity_k_log10":     {"low": -3.0, "high": -1.0},   # 10^U(-3,-1)
    "seasonal_amplitude":  {"low": 1.0,  "high": 1.5},
}

# Fixed parameters (not calibrated)
FIXED = {
    "nticks": 20 * 365,
    "exp_shape": 3,
    "exp_scale": 1.0,
    "inf_mean": 28,
    "inf_sigma": 3,
    "cbr": 29.0,
    "cdr": 7.0,
    "gravity_b": 0.5,
    "gravity_c": 1.5,
    "capacity_safety_factor": 2.0,
    "ri_coverage": 0.80,       # OPV RI at 6 weeks
    "sia_coverage": 0.90,      # SIA biannual 0-5yr
    "ri_age": 42,              # 6 weeks = 42 days
    "sia_period": 180,
    "sia_max_age": 1825,       # 5 years
    "vax_waning_shape": 3,     # Gamma shape for OPV waning
    "vax_waning_scale": 365.0, # Gamma scale → mean ~3yr (1095 days)
    "vax_waning_min": 365,     # Minimum waning duration (days)
    "import_period": 30,
    "import_count": 3,
}


# ============================================================================
# Observed Targets (from MMWR AFP Surveillance)
# ============================================================================

def compute_observed_targets(scenario):
    """Estimate expected zero-AFP-week proportions from MMWR surveillance data.

    For each model district:
    1. Distributes provincial AFP cases proportionally by population
    2. Computes expected weekly AFP rate
    3. Estimates P(zero AFP in a week) via Poisson: exp(-weekly_rate)

    Returns:
        obs_prop_zero_afp: (nnodes,) expected proportion of zero-AFP weeks
        obs_expected_para: (nnodes,) expected annual paralytic cases
    """
    nnodes = len(scenario)
    names = list(scenario.name)
    pops = np.array(scenario.population, dtype=np.float64)

    obs_avg = OBSERVED_WPV1[["Balochistan", "KP", "Sindh", "Punjab", "ICT"]].mean()

    # Province-level model population totals
    prov_model_pop = {}
    for i in range(nnodes):
        prov = DISTRICT_TO_PROVINCE[names[i]]
        prov_model_pop.setdefault(prov, 0.0)
        prov_model_pop[prov] += pops[i]

    # Distribute provincial AFP cases to districts by population share
    obs_expected_para = np.zeros(nnodes)
    for i in range(nnodes):
        prov = DISTRICT_TO_PROVINCE[names[i]]
        obs_expected_para[i] = obs_avg[prov] * pops[i] / prov_model_pop[prov]

    # Poisson approximation for zero-AFP weeks
    weekly_afp_rate = obs_expected_para / 52.0
    obs_prop_zero_afp = np.exp(-weekly_afp_rate)

    return obs_prop_zero_afp, obs_expected_para


# ============================================================================
# CCS Logistic Fit (for diagnostics)
# ============================================================================

def logistic(x, x0, k):
    """Logistic function bounded [0,1], transitioning from 1 to 0."""
    return 1.0 / (1.0 + np.exp(k * (x - x0)))


def fit_ccs_logistic(log_pop, prop_zero):
    """Fit logistic to proportion-zero-AFP vs log10(population).

    Returns (x0, k) or None if fit fails.
    """
    try:
        p0 = [np.median(log_pop), 5.0]
        bounds = ([-np.inf, 0.01], [np.inf, 50.0])
        popt, _ = curve_fit(logistic, log_pop, prop_zero, p0=p0,
                            bounds=bounds, maxfev=10000)
        return popt
    except (RuntimeError, ValueError):
        return None


# ============================================================================
# Scoring Functions
# ============================================================================

def score_afp_provincial(model, scenario):
    """AFP comparison: model expected paralytic vs MMWR observed, by province.

    Computes mean annual model infections per province (post-burn-in),
    converts to expected paralytic cases (×1/200), then calculates
    log-scale MSE against MMWR observed averages.

    Returns:
        float: Mean squared log-error across 5 provinces (lower = better)
    """
    nticks = model.params.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    nyears = nticks // 365
    incidence = model.nodes.newly_infected

    # Aggregate post-burn-in infections by province
    model_by_province = {}
    for i in range(nnodes):
        annual = np.array([
            incidence[y * 365:(y + 1) * 365, i].sum()
            for y in range(BURNIN_YEARS, nyears)
        ])
        prov = DISTRICT_TO_PROVINCE[names[i]]
        model_by_province.setdefault(prov, np.zeros(len(annual)))
        model_by_province[prov] += annual

    obs_avg = OBSERVED_WPV1[["Balochistan", "KP", "Sindh", "Punjab", "ICT"]].mean()

    loss = 0.0
    n_provinces = 0
    for prov in ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]:
        model_para = model_by_province.get(prov, np.zeros(1)).mean() * PARALYSIS_RATIO
        obs_para = obs_avg[prov]
        loss += (np.log(model_para + 1) - np.log(obs_para + 1)) ** 2
        n_provinces += 1

    return loss / n_provinces


def score_ccs_zero_weeks(model, scenario, obs_prop_zero_afp):
    """CCS fitness: proportion of zero-AFP-detection weeks per district.

    For each week and district, computes:
        P(0 AFP cases | W infections) = (199/200)^W

    This accounts for the 1:200 paralysis-to-infection ratio — a week
    with many infections can still produce zero detectable AFP cases.

    Compares model zero-AFP-week proportions to observed targets derived
    from MMWR annual case counts.

    Returns:
        float: Mean squared difference across districts (lower = better)
    """
    nticks = model.params.nticks
    nnodes = len(scenario)
    burnin_ticks = BURNIN_YEARS * 365
    incidence = model.nodes.newly_infected[burnin_ticks:]

    nweeks = incidence.shape[0] // 7
    weekly = incidence[:nweeks * 7].reshape(nweeks, 7, nnodes).sum(axis=1)

    # P(0 AFP in a week) = (1 - 1/200)^W for W weekly infections
    p_zero_afp = (1.0 - PARALYSIS_RATIO) ** weekly
    model_prop_zero_afp = p_zero_afp.mean(axis=0)

    return np.sum((model_prop_zero_afp - obs_prop_zero_afp) ** 2) / nnodes


def score_combined(model, scenario, obs_prop_zero_afp,
                   afp_weight=0.5, ccs_weight=0.5):
    """Combined fitness score (lower = better).

    Returns dict with component and weighted total scores.
    """
    afp = score_afp_provincial(model, scenario)
    ccs = score_ccs_zero_weeks(model, scenario, obs_prop_zero_afp)

    return {
        "afp_loss": afp,
        "ccs_loss": ccs,
        "total_loss": afp_weight * afp + ccs_weight * ccs,
    }


def extract_district_diagnostics(model, scenario, obs_prop_zero_afp):
    """Extract per-district diagnostics for a single run."""
    nticks = model.params.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    nyears = nticks // 365
    burnin_ticks = BURNIN_YEARS * 365
    incidence = model.nodes.newly_infected

    # Weekly incidence (post-burn-in)
    post_inc = incidence[burnin_ticks:]
    nweeks = post_inc.shape[0] // 7
    weekly = post_inc[:nweeks * 7].reshape(nweeks, 7, nnodes).sum(axis=1)

    # Zero-AFP-week proportions
    p_zero_afp = (1.0 - PARALYSIS_RATIO) ** weekly
    model_prop_zero_afp = p_zero_afp.mean(axis=0)

    rows = []
    for i in range(nnodes):
        annual = np.array([
            incidence[y * 365:(y + 1) * 365, i].sum()
            for y in range(BURNIN_YEARS, nyears)
        ])
        rows.append({
            "district": names[i],
            "population": int(scenario.population.iloc[i]),
            "log10_pop": np.log10(scenario.population.iloc[i]),
            "mean_annual_infections": annual.mean(),
            "mean_annual_paralytic": annual.mean() * PARALYSIS_RATIO,
            "model_prop_zero_afp": model_prop_zero_afp[i],
            "obs_prop_zero_afp": obs_prop_zero_afp[i],
        })

    return pd.DataFrame(rows)


# ============================================================================
# Model Runner
# ============================================================================

def precompute_distance_matrix(scenario):
    """Compute pairwise haversine distance matrix (km) between districts."""
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


def build_and_run(beta, gravity_k, seasonal_amplitude, dist_matrix, seed=42):
    """Build and run the Pakistan polio SEIRV model with given parameters.

    Parameters:
        beta: Daily transmission rate (R0 ≈ beta × infectious_period)
        gravity_k: Spatial coupling strength (higher = more inter-district mixing)
        seasonal_amplitude: Monsoon forcing half-amplitude.
            The seasonal multiplier = 1 + amplitude × cos(2π(day-245)/365),
            floored at 0. At amplitude=1.0, transmission ranges from 0 to 2×
            baseline; at 1.5, stronger monsoon peak with near-zero off-season.
        dist_matrix: Pre-computed pairwise distance matrix (km)
        seed: Random seed for reproducibility

    Returns:
        (model, scenario) tuple
    """
    np.random.seed(seed)

    scenario = build_scenario()
    nnodes = len(scenario)

    params = PropertySet({
        "prng_seed": seed,
        "nticks": FIXED["nticks"],
        "beta": beta,
        "exp_shape": FIXED["exp_shape"],
        "exp_scale": FIXED["exp_scale"],
        "inf_mean": FIXED["inf_mean"],
        "inf_sigma": FIXED["inf_sigma"],
        "cbr": FIXED["cbr"],
        "gravity_k": gravity_k,
        "gravity_b": FIXED["gravity_b"],
        "gravity_c": FIXED["gravity_c"],
        "capacity_safety_factor": FIXED["capacity_safety_factor"],
    })

    expdurdist = dists.gamma(shape=params.exp_shape, scale=params.exp_scale)
    infdurdist = dists.normal(loc=params.inf_mean, scale=params.inf_sigma)

    # Monsoon seasonal forcing (peak early September, day 245)
    days = np.arange(365)
    peak_day = 245
    season_365 = 1.0 + seasonal_amplitude * np.cos(
        2 * np.pi * (days - peak_day) / 365
    )
    season_365 = np.maximum(season_365, 0.0)  # Floor at zero
    season_tiled = np.tile(season_365, params.nticks // 365 + 1)[:params.nticks]
    seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

    # Vital dynamics
    birthrate_array = np.full(
        (params.nticks, nnodes), params.cbr, dtype=np.float32
    )
    pyramid, _ = build_pakistan_age_pyramid()
    mortality_array = np.full(
        (params.nticks, nnodes), FIXED["cdr"], dtype=np.float32
    )

    # Build LASER model with V compartment for vaccine tracking
    model = Model(scenario, params, birthrates=birthrate_array,
                  additional_states=["V"])

    # Gravity migration network
    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(
        pops, dist_matrix, 1, 0, FIXED["gravity_b"], FIXED["gravity_c"]
    )
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * gravity_k
    network = row_normalizer(network, 0.2)
    model.network = network

    # Per-patch unreachable fractions (correlated missedness)
    unreach = np.array(scenario.unreachable_frac, dtype=np.float32)

    # OPV vaccine waning distribution (~3 year mean)
    vax_waning_dist = dists.gamma(
        shape=FIXED["vax_waning_shape"],
        scale=FIXED["vax_waning_scale"],
    )

    # Endemic corridor patches
    name_to_idx = {n: i for i, n in enumerate(scenario.name)}
    endemic_names = [
        "Quetta", "Peshawar", "Karachi", "Hyderabad",
        "Bannu", "N_Waziristan", "S_Waziristan", "DI_Khan",
    ]
    endemic_patches = [name_to_idx[n] for n in endemic_names]

    # Assemble SEIRV components (same structure as polio_seir_10patch.py)
    susceptible = SEIR.Susceptible(model)
    exposed = SEIR.Exposed(model, expdurdist, infdurdist)
    infectious = SEIR.Infectious(model, infdurdist)
    recovered = SEIR.Recovered(model)

    vax_compartment = VaccinatedCompartment(
        model, vax_waning_dist, wandurmin=FIXED["vax_waning_min"]
    )
    vaccination = PerPatchVaccinationSEIRV(
        model, vax_compartment, unreach,
        ri_coverage=FIXED["ri_coverage"],
        sia_coverage=FIXED["sia_coverage"],
        ri_age=FIXED["ri_age"],
        sia_period=FIXED["sia_period"],
        sia_max_age=FIXED["sia_max_age"],
    )

    model.components = [
        susceptible,
        exposed,
        infectious,
        recovered,
        vax_compartment,
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
        MortalityByCDR(model, mortalityrates=mortality_array),
        vaccination,
        PatchImportation(model, infdurdist, endemic_patches,
                         period=FIXED["import_period"],
                         count=FIXED["import_count"]),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    ]

    initialize_ages(model, pyramid)
    model.run("Calibration")

    return model, scenario


# ============================================================================
# Latin Hypercube Sampling
# ============================================================================

def generate_samples(n_samples, seed=42):
    """Generate space-filling parameter samples using Latin Hypercube Sampling.

    gravity_k is sampled log-uniformly: 10^U(-3, -1).

    Returns list of dicts with keys: beta, gravity_k, seasonal_amplitude.
    """
    try:
        from scipy.stats.qmc import LatinHypercube
        sampler = LatinHypercube(d=3, seed=seed)
        unit = sampler.random(n=n_samples)
    except ImportError:
        # Fallback to stratified random if scipy.stats.qmc unavailable
        rng = np.random.RandomState(seed)
        unit = rng.random((n_samples, 3))

    samples = []
    for i in range(n_samples):
        beta = (PARAM_RANGES["beta"]["low"]
                + unit[i, 0] * (PARAM_RANGES["beta"]["high"]
                                - PARAM_RANGES["beta"]["low"]))
        gk_log = (PARAM_RANGES["gravity_k_log10"]["low"]
                  + unit[i, 1] * (PARAM_RANGES["gravity_k_log10"]["high"]
                                  - PARAM_RANGES["gravity_k_log10"]["low"]))
        seasonal_amp = (PARAM_RANGES["seasonal_amplitude"]["low"]
                        + unit[i, 2] * (PARAM_RANGES["seasonal_amplitude"]["high"]
                                        - PARAM_RANGES["seasonal_amplitude"]["low"]))

        samples.append({
            "beta": round(beta, 5),
            "gravity_k": round(10 ** gk_log, 6),
            "seasonal_amplitude": round(seasonal_amp, 4),
        })

    return samples


# ============================================================================
# Calibration Loop
# ============================================================================

def run_calibration(n_samples=32, seed=42):
    """Run calibration: sample parameter space, simulate, score, rank.

    Returns:
        results_df: DataFrame of all samples ranked by total_loss
    """
    print("=" * 70)
    print("PAKISTAN POLIO MODEL CALIBRATION")
    print("=" * 70)
    print(f"\nParameter ranges:")
    for name, bounds in PARAM_RANGES.items():
        if name == "gravity_k_log10":
            print(f"  gravity_k:           10^U({bounds['low']}, {bounds['high']})"
                  f"  = [{10**bounds['low']:.4f}, {10**bounds['high']:.2f}]")
        else:
            print(f"  {name + ':':<23} U({bounds['low']}, {bounds['high']})")
    print(f"\nFixed parameters:")
    print(f"  nticks={FIXED['nticks']}  cbr={FIXED['cbr']}  cdr={FIXED['cdr']}")
    print(f"  inf_mean={FIXED['inf_mean']}d  gravity_b={FIXED['gravity_b']}"
          f"  gravity_c={FIXED['gravity_c']}")
    print(f"\nSamples: {n_samples}  |  Seed: {seed}")
    print(f"Burn-in: {BURNIN_YEARS} years  |  Analysis: years {BURNIN_YEARS}"
          f"–{FIXED['nticks'] // 365}")

    # Pre-compute shared data
    scenario_template = build_scenario()
    dist_matrix = precompute_distance_matrix(scenario_template)
    obs_prop_zero_afp, obs_expected_para = compute_observed_targets(
        scenario_template
    )
    obs_avg = OBSERVED_WPV1[["Balochistan", "KP", "Sindh", "Punjab", "ICT"]].mean()
    obs_total_para = obs_avg.sum()

    print(f"\nObserved targets (MMWR 2020-2025 average):")
    print(f"  National: {obs_total_para:.1f} paralytic cases/year")
    names = list(scenario_template.name)
    print(f"\n  {'District':<14} {'Exp AFP/yr':<12} {'Exp prop_zero':<14}")
    print(f"  {'-'*40}")
    for i, name in enumerate(names):
        print(f"  {name:<14} {obs_expected_para[i]:<12.2f}"
              f" {obs_prop_zero_afp[i]:<14.4f}")

    # Generate parameter samples
    samples = generate_samples(n_samples, seed)

    # Main calibration loop
    results = []
    best_diagnostics = None
    best_loss = float("inf")
    t_total = time.time()

    for idx, params in enumerate(samples):
        print(f"\n{'─' * 60}")
        print(f"Sample {idx + 1}/{n_samples}  "
              f"beta={params['beta']:.4f}  "
              f"gravity_k={params['gravity_k']:.5f}  "
              f"seasonal_amp={params['seasonal_amplitude']:.3f}")

        t0 = time.time()
        try:
            model, scenario = build_and_run(
                params["beta"],
                params["gravity_k"],
                params["seasonal_amplitude"],
                dist_matrix,
                seed=seed + idx,
            )

            scores = score_combined(model, scenario, obs_prop_zero_afp)

            # Summary diagnostics
            nticks = model.params.nticks
            nyears = nticks // 365
            analysis_years = nyears - BURNIN_YEARS
            incidence = model.nodes.newly_infected
            total_inf = incidence[BURNIN_YEARS * 365:].sum()
            avg_annual_inf = total_inf / analysis_years
            avg_annual_para = avg_annual_inf * PARALYSIS_RATIO

            elapsed = time.time() - t0

            result = {
                "sample_id": idx,
                **params,
                "afp_loss": scores["afp_loss"],
                "ccs_loss": scores["ccs_loss"],
                "total_loss": scores["total_loss"],
                "avg_annual_infections": avg_annual_inf,
                "avg_annual_paralytic": avg_annual_para,
                "elapsed_sec": round(elapsed, 1),
            }
            results.append(result)

            # Track best for detailed diagnostics
            if scores["total_loss"] < best_loss:
                best_loss = scores["total_loss"]
                best_diagnostics = extract_district_diagnostics(
                    model, scenario, obs_prop_zero_afp
                )
                best_params = params.copy()

            print(f"  AFP loss: {scores['afp_loss']:.4f}  "
                  f"CCS loss: {scores['ccs_loss']:.4f}  "
                  f"Total: {scores['total_loss']:.4f}")
            print(f"  Avg paralytic/yr: {avg_annual_para:.1f}  "
                  f"(target: ~{obs_total_para:.0f}/yr)")
            print(f"  Elapsed: {elapsed:.0f}s")

        except Exception as e:
            elapsed = time.time() - t0
            print(f"  FAILED after {elapsed:.0f}s: {e}")
            results.append({
                "sample_id": idx,
                **params,
                "afp_loss": float("inf"),
                "ccs_loss": float("inf"),
                "total_loss": float("inf"),
                "avg_annual_infections": 0,
                "avg_annual_paralytic": 0,
                "elapsed_sec": round(elapsed, 1),
            })

        # Free memory
        del model, scenario
        gc.collect()

    total_elapsed = time.time() - t_total

    # ================================================================
    # Rank Results
    # ================================================================

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("total_loss").reset_index(drop=True)
    results_df["rank"] = range(1, len(results_df) + 1)

    print("\n" + "=" * 70)
    print("CALIBRATION RESULTS — Ranked by Total Loss")
    print("=" * 70)

    display_cols = [
        "rank", "beta", "gravity_k", "seasonal_amplitude",
        "afp_loss", "ccs_loss", "total_loss", "avg_annual_paralytic",
    ]
    print(results_df[display_cols].to_string(
        index=False,
        formatters={
            "rank": "{:>4d}".format,
            "beta": "{:.4f}".format,
            "gravity_k": "{:.5f}".format,
            "seasonal_amplitude": "{:.3f}".format,
            "afp_loss": "{:.4f}".format,
            "ccs_loss": "{:.4f}".format,
            "total_loss": "{:.4f}".format,
            "avg_annual_paralytic": "{:.1f}".format,
        },
    ))

    # Best parameters
    best = results_df.iloc[0]
    print(f"\n{'─' * 50}")
    print(f"BEST PARAMETERS (rank 1):")
    print(f"  beta:                {best['beta']:.4f}")
    print(f"  gravity_k:           {best['gravity_k']:.5f}")
    print(f"  seasonal_amplitude:  {best['seasonal_amplitude']:.3f}")
    print(f"  Total loss:          {best['total_loss']:.4f}")
    print(f"    AFP component:     {best['afp_loss']:.4f}")
    print(f"    CCS component:     {best['ccs_loss']:.4f}")
    print(f"  Avg paralytic/yr:    {best['avg_annual_paralytic']:.1f}"
          f"  (target: ~{obs_total_para:.0f})")

    # Best-run district diagnostics
    if best_diagnostics is not None:
        print(f"\nBest run — per-district diagnostics:")
        print(f"  {'District':<14} {'Pop':<10} {'Inf/yr':<10} {'Para/yr':<10}"
              f" {'Model p0':<10} {'Obs p0':<10}")
        print(f"  {'-'*64}")
        for _, row in best_diagnostics.iterrows():
            print(f"  {row['district']:<14} {row['population']:<10}"
                  f" {row['mean_annual_infections']:<10.0f}"
                  f" {row['mean_annual_paralytic']:<10.1f}"
                  f" {row['model_prop_zero_afp']:<10.4f}"
                  f" {row['obs_prop_zero_afp']:<10.4f}")

    print(f"\nTotal calibration time: {total_elapsed / 60:.1f} minutes"
          f" ({total_elapsed / n_samples:.0f}s per sample)")

    # Save results
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "calibration_results.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")

    if best_diagnostics is not None:
        diag_path = outdir / "calibration_best_diagnostics.csv"
        best_diagnostics.to_csv(diag_path, index=False)
        print(f"Best-run diagnostics saved to {diag_path}")

    # Generate summary plot
    plot_calibration_summary(
        results_df, best_diagnostics, obs_prop_zero_afp,
        scenario_template, outdir
    )

    return results_df


# ============================================================================
# Diagnostic Plots
# ============================================================================

def plot_calibration_summary(results_df, best_diagnostics, obs_prop_zero_afp,
                             scenario, outdir):
    """Generate 4-panel calibration summary figure.

    (A) Parameter samples colored by total loss
    (B) CCS curve: prop_zero_afp vs log10(pop) — best run vs observed
    (C) Provincial AFP comparison for best run
    (D) Loss distribution
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    valid = results_df[results_df["total_loss"] < float("inf")]
    if len(valid) == 0:
        plt.close(fig)
        return

    # --- (A) Parameter space: beta vs gravity_k, colored by loss ---
    ax = axes[0, 0]
    sc = ax.scatter(
        valid["beta"],
        valid["gravity_k"],
        c=valid["total_loss"],
        cmap="viridis_r",
        s=60,
        edgecolor="black",
        linewidth=0.5,
    )
    # Mark best
    best = valid.iloc[0]
    ax.scatter(best["beta"], best["gravity_k"],
               c="red", s=150, marker="*", zorder=5, label="Best")
    ax.set_xlabel("beta (transmission rate)")
    ax.set_ylabel("gravity_k (spatial coupling)")
    ax.set_yscale("log")
    ax.set_title("(A) Parameter Space (color = total loss)")
    ax.legend(loc="upper right")
    fig.colorbar(sc, ax=ax, shrink=0.8, label="Total loss")

    # --- (B) CCS curve: zero-AFP-week proportion vs population ---
    ax = axes[0, 1]
    if best_diagnostics is not None:
        log_pop = best_diagnostics["log10_pop"].values
        model_pz = best_diagnostics["model_prop_zero_afp"].values
        obs_pz = best_diagnostics["obs_prop_zero_afp"].values
        district_names = best_diagnostics["district"].values

        ax.scatter(log_pop, obs_pz, c="black", s=80, marker="s",
                   zorder=5, label="Observed (MMWR-derived)")
        ax.scatter(log_pop, model_pz, c="steelblue", s=60,
                   zorder=4, label="Model (best run)")

        # Logistic fits
        obs_fit = fit_ccs_logistic(log_pop, obs_pz)
        model_fit = fit_ccs_logistic(log_pop, model_pz)
        x_grid = np.linspace(log_pop.min() - 0.2, log_pop.max() + 0.2, 100)
        if obs_fit is not None:
            ax.plot(x_grid, logistic(x_grid, *obs_fit), "k--", alpha=0.5,
                    label=f"Obs fit (x0={obs_fit[0]:.2f})")
        if model_fit is not None:
            ax.plot(x_grid, logistic(x_grid, *model_fit), "b-", alpha=0.5,
                    label=f"Model fit (x0={model_fit[0]:.2f})")

        # Label districts
        for i, name in enumerate(district_names):
            ax.annotate(name, (log_pop[i], model_pz[i]),
                        fontsize=6, alpha=0.7,
                        xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("log10(population)")
    ax.set_ylabel("Proportion of zero-AFP weeks")
    ax.set_title("(B) CCS Analysis — Best Run vs Observed")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(-0.05, 1.05)

    # --- (C) Provincial AFP comparison (best run) ---
    ax = axes[1, 0]
    if best_diagnostics is not None:
        prov_order = ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]
        obs_avg = OBSERVED_WPV1[prov_order].mean()
        names_list = list(scenario.name)
        nnodes = len(scenario)

        model_prov_para = {}
        for _, row in best_diagnostics.iterrows():
            prov = DISTRICT_TO_PROVINCE[row["district"]]
            model_prov_para.setdefault(prov, 0.0)
            model_prov_para[prov] += row["mean_annual_paralytic"]

        x = np.arange(len(prov_order))
        width = 0.35
        obs_vals = [obs_avg[p] for p in prov_order]
        model_vals = [model_prov_para.get(p, 0) for p in prov_order]

        prov_colors = {
            "Balochistan": "#e74c3c", "KP": "#e67e22",
            "Sindh": "#3498db", "Punjab": "#2ecc71", "ICT": "#9b59b6",
        }
        bar_colors = [prov_colors[p] for p in prov_order]

        ax.bar(x - width / 2, obs_vals, width, color="black", alpha=0.7,
               label="Observed (MMWR avg)")
        ax.bar(x + width / 2, model_vals, width, color=bar_colors, alpha=0.7,
               label="Model (best, ×1/200)")
        ax.set_xticks(x)
        ax.set_xticklabels(prov_order, rotation=30, ha="right")
        ax.set_ylabel("Paralytic Cases / Year")
        ax.set_title("(C) Provincial AFP: Best Run vs Observed")
        ax.legend(fontsize=8)

    # --- (D) Loss distribution ---
    ax = axes[1, 1]
    ax.hist(valid["total_loss"], bins=min(20, len(valid)), color="steelblue",
            edgecolor="black", alpha=0.7)
    ax.axvline(valid["total_loss"].iloc[0], color="red", ls="--", lw=2,
               label=f"Best: {valid['total_loss'].iloc[0]:.4f}")
    ax.set_xlabel("Total Loss")
    ax.set_ylabel("Count")
    ax.set_title("(D) Loss Distribution Across Samples")
    ax.legend()

    plt.tight_layout()
    outpath = outdir / "calibration_summary.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Summary plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# CLI Entry Point
# ============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate Pakistan polio SEIR model against AFP surveillance"
    )
    parser.add_argument(
        "--n-samples", type=int, default=32,
        help="Number of Latin Hypercube parameter samples (default: 32)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for sampling and simulation (default: 42)"
    )
    parser.add_argument(
        "--afp-weight", type=float, default=0.5,
        help="Weight for AFP provincial loss component (default: 0.5)"
    )
    parser.add_argument(
        "--ccs-weight", type=float, default=0.5,
        help="Weight for CCS zero-week loss component (default: 0.5)"
    )
    args = parser.parse_args()

    results = run_calibration(n_samples=args.n_samples, seed=args.seed)
