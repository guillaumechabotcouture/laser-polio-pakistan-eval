#!/usr/bin/env python3
"""
Pakistan Polio SEIR Model — 20 Districts with Heterogeneous Populations

Specifications:
    - 20 districts with populations ranging 50k–2M
    - SEIR dynamics: R0≈6, 3-day latent period, 28-day infectious period
    - Gravity-model spatial coupling
    - Monsoon-season transmission forcing (peak Jul–Oct)
    - OPV routine immunization at 80% coverage
    - SIA campaigns every 6 months at 80% coverage
    - Importation: 3 infections every 60 days for the first 5 years
    - 20-year simulation with 10-year burn-in
    - Output: weekly incidence per district, proportion of zero-incidence weeks

Usage:
    /opt/anaconda3/bin/python3 scripts/polio_seir_20district.py
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
from laser.core.migration import gravity, row_normalizer

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add scripts dir to path for custom_components
sys.path.insert(0, str(Path(__file__).parent))
from custom_components import PerPatchVaccination, PatchImportation


# ============================================================================
# 20-District Configuration
#
# Heterogeneous populations from 50k (Chaman) to 2M (Karachi).
# Geographic coordinates approximate real district centroids.
#
# unreachable_frac models the zero-dose population (children in hard-to-reach,
# refusal, or inaccessible communities) who are never vaccinated by RI or SIA.
# This is the key structural parameter enabling sustained endemic transmission
# even with high nominal coverage. Districts where unreachable > 1/R0 (≈0.167)
# can sustain local transmission chains; others see only sporadic importations.
#
# Unreachable fractions reflect Pakistan's actual coverage gaps:
#   - Punjab/ICT: 3-5% zero-dose (well-served, high coverage)
#   - Sindh urban: 12-18% (urban slums, migrant communities)
#   - KP settled: 15-20% (moderate access challenges)
#   - KP tribal/border: 40-55% (security, refusal, inaccessibility)
#   - Balochistan border: 40-55% (geographic isolation, refusal)
#
# ============================================================================

DISTRICTS = pd.DataFrame([
    # name,            province,       pop,     lat,      lon,  unreach
    # Punjab/ICT — high coverage, low unreachable
    ("Lahore",         "Punjab",      1900000, 31.5204, 74.3587, 0.03),
    ("Faisalabad",     "Punjab",      1800000, 31.4504, 73.1350, 0.03),
    ("Rawalpindi",     "Punjab",      1400000, 33.5651, 73.0169, 0.03),
    ("Multan",         "Punjab",      1500000, 30.1575, 71.5249, 0.04),
    ("Bahawalpur",     "Punjab",       900000, 29.3544, 71.6911, 0.05),
    ("Islamabad",      "ICT",          450000, 33.6844, 73.0479, 0.02),
    # Sindh — moderate unreachable
    ("Karachi",        "Sindh",       2000000, 24.8607, 67.0011, 0.15),
    ("Hyderabad",      "Sindh",       1200000, 25.3960, 68.3578, 0.12),
    ("Sukkur",         "Sindh",        500000, 27.7244, 68.8228, 0.18),
    ("Jacobabad",      "Sindh",        100000, 28.2769, 68.4514, 0.25),
    # KP — settled districts
    ("Peshawar",       "KP",          1700000, 34.0150, 71.5249, 0.18),
    ("Mardan",         "KP",           700000, 34.1988, 72.0404, 0.15),
    ("Swat",           "KP",           400000, 35.2220, 72.3528, 0.20),
    # KP — border/tribal (high unreachable, endemic corridor)
    ("DI_Khan",        "KP",          1000000, 31.8320, 70.9015, 0.40),
    ("Bannu",          "KP",           600000, 32.9860, 70.6027, 0.40),
    ("N_Waziristan",   "KP",           350000, 32.3045, 69.8597, 0.55),
    ("S_Waziristan",   "KP",           300000, 32.0711, 69.5310, 0.55),
    # Balochistan — high unreachable, border districts
    ("Quetta",         "Balochistan",  800000, 30.1798, 66.9750, 0.40),
    ("Zhob",           "Balochistan",  200000, 31.3413, 69.4494, 0.45),
    ("Chaman",         "Balochistan",   50000, 30.9210, 66.4536, 0.55),
], columns=["name", "province", "population", "lat", "lon",
           "unreachable_frac"])


# ============================================================================
# Model Parameters
# ============================================================================

PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": 20 * 365,         # 20-year simulation (first 10 = burn-in)
    "beta": 6.0 / 28.0,         # R0 = 6, D_inf = 28d → beta = R0/D ≈ 0.2143
    "exp_shape": 3,              # Exposed: gamma(shape=3, scale=1) → mean = 3 days
    "exp_scale": 1.0,
    "inf_mean": 28,              # Infectious: normal(28, 3) → mean = 28 days
    "inf_sigma": 3,
    "cbr": 29.0,                 # Crude birth rate per 1000/year (Pakistan avg)
    "gravity_k": 0.005,          # Gravity coupling strength
    "gravity_b": 0.5,            # Destination population exponent
    "gravity_c": 1.5,            # Distance decay exponent
    "capacity_safety_factor": 2.5,
})

BURNIN_YEARS = 10
CDR = 7.0  # Crude death rate per 1000/year


# ============================================================================
# Monsoon Seasonal Forcing (Jul–Oct peak)
# ============================================================================

def build_monsoon_seasonality(nticks, nnodes):
    """Build 365-day seasonal profile with monsoon peak (Jul–Oct).

    Cosine-based forcing centered on day 245 (early September)
    with ±30% amplitude creating ~1.6x peak-to-trough ratio.
    """
    days = np.arange(365)
    peak_day = 245   # ~Sep 2, center of Jul–Oct monsoon window
    amplitude = 0.30
    beta_season = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
    season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
    return ValuesMap.from_timeseries(season_tiled, nnodes), beta_season


# ============================================================================
# Age Pyramid (young population, ~65-year life expectancy)
# ============================================================================

def build_age_pyramid():
    """Approximate Pakistan age distribution — heavily weighted toward youth."""
    mean_life = 65.0
    ages = np.arange(100)
    stable_age_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_age_dist = np.maximum(stable_age_dist, 1)
    return AliasedDistribution(stable_age_dist), stable_age_dist


# ============================================================================
# Build and Run Model
# ============================================================================

def build_scenario():
    """Build GeoDataFrame scenario from district configuration."""
    geometry = [Point(lon, lat) for lat, lon in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(DISTRICTS, geometry=geometry, crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))

    # Initial compartments: start near endemic equilibrium (S/N ≈ 1/R0 ≈ 17%)
    # with active infections so the burn-in starts from a steady state rather
    # than an explosive epidemic that depletes all susceptibles.
    scenario["I"] = np.maximum(
        np.round(0.001 * scenario["population"]).astype(np.uint32), 1)
    scenario["E"] = np.round(0.0003 * scenario["population"]).astype(np.uint32)
    scenario["R"] = np.round(0.83 * scenario["population"]).astype(np.uint32)
    scenario["S"] = (scenario["population"] - scenario["E"]
                     - scenario["I"] - scenario["R"]).astype(np.uint32)
    return scenario


def initialize_ages(model, pyramid):
    """Set initial agent DOBs from age pyramid to avoid all-newborn artifacts."""
    count = model.people.count
    ages_years = pyramid.sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, 99)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def run_model():
    """Build, configure, and run the 20-district Pakistan polio SEIR model."""
    scenario = build_scenario()
    nnodes = len(scenario)

    # Duration distributions
    expdurdist = dists.gamma(shape=PARAMS.exp_shape, scale=PARAMS.exp_scale)
    infdurdist = dists.normal(loc=PARAMS.inf_mean, scale=PARAMS.inf_sigma)

    # Seasonal forcing
    seasonality, beta_season_365 = build_monsoon_seasonality(PARAMS.nticks, nnodes)

    # Birth rates: per-1000/year (BirthsByCBR divides by 1000 internally)
    birthrate_array = np.full((PARAMS.nticks, nnodes), PARAMS.cbr, dtype=np.float32)
    assert np.all(birthrate_array >= 1) and np.all(birthrate_array <= 60), \
        f"Birthrates must be per-1000/year, got {birthrate_array.min():.4f}–{birthrate_array.max():.4f}"

    # Mortality: CDR per-1000/year
    mortality_array = np.full((PARAMS.nticks, nnodes), CDR, dtype=np.float32)

    # Age pyramid
    pyramid, _ = build_age_pyramid()

    # Build model
    model = Model(scenario, PARAMS, birthrates=birthrate_array)

    # --- Gravity migration network ---
    from laser.core.migration import distance as haversine_distance
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

    # --- Vaccination: 80% RI + 80% SIA every 6 months ---
    # Coverage applies AMONG REACHABLE agents only; unreachable are never vaccinated
    ri_cov = np.full(nnodes, 0.80, dtype=np.float32)
    sia_cov = np.full(nnodes, 0.80, dtype=np.float32)
    unreachable = np.array(scenario.unreachable_frac, dtype=np.float32)

    # --- Importation: 3 infections every 60 days for first 5 years ---
    # Seed in 3 endemic corridor districts (1 infection per patch per event = 3 total)
    # Runs for first 5 years only to establish endemic transmission chains
    name_to_idx = {n: i for i, n in enumerate(scenario.name)}
    import_patches = [
        name_to_idx["Quetta"],
        name_to_idx["Peshawar"],
        name_to_idx["N_Waziristan"],
    ]

    # --- Assemble components (order matters) ---
    model.components = [
        SEIR.Susceptible(model),
        SEIR.Exposed(model, expdurdist, infdurdist),
        SEIR.Infectious(model, infdurdist),
        SEIR.Recovered(model),
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
        MortalityByCDR(model, mortalityrates=mortality_array),
        PerPatchVaccination(model, ri_cov, sia_cov, unreachable,
                            ri_period=7, ri_age=42, sia_period=180),
        PatchImportation(model, infdurdist, import_patches,
                         period=60, count=1, end_tick=5 * 365),  # 3 inf/60d, first 5 years
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    ]

    # Initialize ages from pyramid (prevents all-newborn artifact)
    initialize_ages(model, pyramid)

    total_pop = scenario.population.sum()
    print(f"Running 20-year Pakistan polio SEIR simulation ({BURNIN_YEARS}yr burn-in)...")
    print(f"  Districts: {nnodes}")
    print(f"  Total population: {total_pop:,}")
    print(f"  R0 = {PARAMS.beta * PARAMS.inf_mean:.1f}, "
          f"latent = {PARAMS.exp_shape * PARAMS.exp_scale:.0f}d, "
          f"infectious = {PARAMS.inf_mean}d")
    avg_unreach = np.average(unreachable, weights=pops)
    print(f"  Vaccination: RI 80%, SIA 80% every 180d "
          f"(pop-weighted avg unreachable: {avg_unreach:.1%})")
    print(f"  Importation: 3 infections / 60 days in Quetta, Peshawar, N_Waziristan "
          f"(first 5 years only)")
    print(f"  Ticks: {PARAMS.nticks}")
    model.run("Pakistan Polio 20-District SEIR")

    return model, scenario, beta_season_365


# ============================================================================
# Output: Weekly Incidence per District
# ============================================================================

def compute_weekly_incidence(model, scenario):
    """Compute weekly incidence per district for the post-burn-in period.

    Returns a DataFrame with columns: year, week, district, incidence
    """
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    burnin_ticks = BURNIN_YEARS * 365

    incidence = model.nodes.newly_infected  # shape: (nticks, nnodes)
    post_inc = incidence[burnin_ticks:]

    # Reshape into complete weeks
    nweeks = post_inc.shape[0] // 7
    weekly = post_inc[:nweeks * 7].reshape(nweeks, 7, nnodes).sum(axis=1)

    # Build DataFrame
    rows = []
    for w in range(nweeks):
        year = BURNIN_YEARS + w * 7 / 365.0
        for d in range(nnodes):
            rows.append({
                "week": w + 1,
                "year": round(year, 3),
                "district": names[d],
                "incidence": int(weekly[w, d]),
            })
    return pd.DataFrame(rows)


# ============================================================================
# Output: Zero-Incidence Week Analysis
# ============================================================================

def analyze_zero_incidence(weekly_df, scenario):
    """Compute and print proportion of zero-incidence weeks per district."""
    names = list(scenario.name)
    pops = list(scenario.population)
    nweeks = weekly_df["week"].nunique()

    print("\n" + "=" * 75)
    print(f"ZERO-INCIDENCE WEEK ANALYSIS (post burn-in, {nweeks} weeks)")
    print("=" * 75)
    print(f"\n{'District':<16} {'Pop':>10} {'Total Inc':>10} "
          f"{'Weeks>0':>8} {'Weeks=0':>8} {'% Zero':>8}")
    print("-" * 70)

    total_zero_weeks = 0
    total_district_weeks = 0

    for i, name in enumerate(names):
        dist_data = weekly_df[weekly_df["district"] == name]
        total_inc = dist_data["incidence"].sum()
        weeks_nonzero = (dist_data["incidence"] > 0).sum()
        weeks_zero = nweeks - weeks_nonzero
        pct_zero = 100.0 * weeks_zero / nweeks

        total_zero_weeks += weeks_zero
        total_district_weeks += nweeks

        print(f"{name:<16} {pops[i]:>10,} {total_inc:>10,} "
              f"{weeks_nonzero:>8} {weeks_zero:>8} {pct_zero:>7.1f}%")

    overall_pct = 100.0 * total_zero_weeks / total_district_weeks
    print("-" * 70)
    print(f"{'OVERALL':<16} {'':>10} {'':>10} "
          f"{'':>8} {total_zero_weeks:>8} {overall_pct:>7.1f}%")

    return nweeks


# ============================================================================
# Summary Statistics
# ============================================================================

def print_summary(model, scenario):
    """Print key summary statistics for verification."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    nyears = nticks // 365
    analysis_years = nyears - BURNIN_YEARS

    incidence = model.nodes.newly_infected
    R0 = PARAMS.beta * PARAMS.inf_mean

    print("\n" + "=" * 75)
    print(f"SIMULATION SUMMARY (post burn-in: years {BURNIN_YEARS}–{nyears})")
    print("=" * 75)

    print(f"\n{'District':<16} {'Pop':>10} {'AnnInf(med)':>12} "
          f"{'AnnInf(avg)':>12} {'S/N(final)':>10}")
    print("-" * 70)

    total_infections = 0
    for i in range(nnodes):
        annual = np.array([
            incidence[y * 365:(y + 1) * 365, i].sum()
            for y in range(BURNIN_YEARS, nyears)
        ])
        median_ann = np.median(annual)
        mean_ann = annual.mean()
        total_post = annual.sum()
        total_infections += total_post

        # Final susceptible fraction
        ft = nticks - 1
        S_f = float(model.nodes.S[ft, i])
        N_f = float(model.nodes.S[ft, i] + model.nodes.E[ft, i]
                     + model.nodes.I[ft, i] + model.nodes.R[ft, i])
        sn = S_f / max(N_f, 1)

        print(f"{names[i]:<16} {scenario.population.iloc[i]:>10,} "
              f"{median_ann:>12,.0f} {mean_ann:>12,.0f} {sn:>10.3f}")

    total_pop = scenario.population.sum()
    print(f"\n  Total population: {total_pop:,}")
    print(f"  Total infections (years {BURNIN_YEARS}–{nyears}): {total_infections:,.0f}")
    print(f"  Avg annual infections: {total_infections / analysis_years:,.0f}")
    print(f"  Annual rate per million: "
          f"{total_infections / analysis_years / total_pop * 1e6:.0f}")
    print(f"  R0 = {R0:.1f}")

    # Compartment integrity check
    print("\nCompartment integrity check:")
    for tick in [0, nticks // 2, nticks - 1]:
        S = model.nodes.S[tick].sum()
        E = model.nodes.E[tick].sum()
        I = model.nodes.I[tick].sum()
        R = model.nodes.R[tick].sum()
        print(f"  Tick {tick:>5}: S={S:>10,}  E={E:>6,}  "
              f"I={I:>6,}  R={R:>10,}  N={S + E + I + R:>10,}")

    # Population trajectory
    print(f"\nPopulation trajectory:")
    for y in range(0, nyears + 1, 5):
        t = min(y * 365, nticks - 1)
        pop_t = int(model.nodes.S[t].sum() + model.nodes.E[t].sum()
                     + model.nodes.I[t].sum() + model.nodes.R[t].sum())
        print(f"  Year {y:>2}: {pop_t:>12,}")
    print(f"  Agent capacity: {model.people.capacity:,}, "
          f"Active: {model.people.count:,}")


# ============================================================================
# Diagnostic Plots
# ============================================================================

def plot_diagnostics(model, scenario, beta_season_365, outdir):
    """Generate 6-panel diagnostic figure."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    provinces = list(scenario.province)
    R0 = PARAMS.beta * PARAMS.inf_mean
    burnin_ticks = BURNIN_YEARS * 365
    nyears = nticks // 365

    prov_colors = {
        "Balochistan": "#e74c3c", "KP": "#e67e22",
        "Sindh": "#3498db", "Punjab": "#2ecc71", "ICT": "#9b59b6",
    }
    district_colors = [prov_colors.get(p, "gray") for p in provinces]

    incidence = model.nodes.newly_infected
    pops = np.array(scenario.population, dtype=np.float64)
    annual_inc = np.zeros((nyears, nnodes))
    for y in range(nyears):
        annual_inc[y] = incidence[y * 365:(y + 1) * 365].sum(axis=0)

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))

    # --- (A) Monsoon seasonal forcing ---
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, beta_season_365, "b-", lw=2)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
    ax.axvspan(182, 304, alpha=0.12, color="blue", label="Monsoon (Jul–Oct)")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Seasonal Multiplier")
    ax.set_title("(A) Monsoon Seasonal Forcing")
    ax.legend(fontsize=8)

    # --- (B) Annual incidence rate per million ---
    ax = axes[0, 1]
    years_x = np.arange(nyears) + 0.5
    for i in range(nnodes):
        rate = annual_inc[:, i] / pops[i] * 1e6
        if rate.max() > 1:
            ax.plot(years_x, rate, color=district_colors[i], alpha=0.7,
                    lw=1, label=names[i])
    ax.axvline(BURNIN_YEARS, color="black", ls=":", alpha=0.4, label="Burn-in")
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual Incidence per Million")
    ax.set_title("(B) Incidence Rate by District")
    ax.legend(fontsize=5, ncol=3, loc="upper right")

    # --- (C) R_eff over time ---
    ax = axes[0, 2]
    sample30 = np.arange(0, nticks, 30)
    time30 = sample30 / 365.0
    for i in range(nnodes):
        S = model.nodes.S[sample30, i].astype(np.float64)
        E = model.nodes.E[sample30, i].astype(np.float64)
        I = model.nodes.I[sample30, i].astype(np.float64)
        R = model.nodes.R[sample30, i].astype(np.float64)
        N = np.maximum(S + E + I + R, 1.0)
        ax.plot(time30, R0 * S / N, color=district_colors[i], alpha=0.6,
                lw=0.8, label=names[i])
    ax.axhline(1.0, color="black", ls="-", lw=2, alpha=0.4, label="R_eff=1")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Year")
    ax.set_ylabel("R_eff = R0 × S/N")
    ax.set_title("(C) Effective Reproduction Number")
    ax.set_ylim(0, 4)
    ax.legend(fontsize=4, ncol=3, loc="upper right")

    # --- (D) Heatmap: annual incidence per million (log) ---
    ax = axes[1, 0]
    log_rate = np.log10(annual_inc / pops[None, :] * 1e6 + 1)
    im = ax.imshow(log_rate.T, aspect="auto", cmap="YlOrRd", origin="lower",
                   extent=[0, nyears, -0.5, nnodes - 0.5])
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Year")
    ax.axvline(BURNIN_YEARS, color="white", ls="--", lw=1.5)
    ax.set_title("(D) Annual Incidence per Million (log10)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="log10(per million + 1)")

    # --- (E) Post-burn-in weekly incidence ---
    ax = axes[1, 1]
    post_inc = incidence[burnin_ticks:]
    post_weeks = post_inc.shape[0] // 7
    post_weekly = post_inc[:post_weeks * 7].reshape(post_weeks, 7, nnodes).sum(axis=1)
    weeks_x = np.arange(post_weeks) / 52.0 + BURNIN_YEARS
    for i in range(nnodes):
        if post_weekly[:, i].max() > 0:
            ax.plot(weeks_x, post_weekly[:, i], color=district_colors[i],
                    alpha=0.7, lw=0.6, label=names[i])
    for yr in range(BURNIN_YEARS, nyears):
        ax.axvspan(yr + 182 / 365, yr + 304 / 365, alpha=0.04, color="blue")
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly Infections")
    ax.set_title("(E) Post-Burn-in Weekly Incidence (blue=monsoon)")
    ax.legend(fontsize=4, ncol=3, loc="upper right")

    # --- (F) Zero-incidence proportion by district ---
    ax = axes[1, 2]
    zero_pcts = []
    for i in range(nnodes):
        weekly_i = post_weekly[:, i]
        pct = 100.0 * (weekly_i == 0).sum() / post_weeks
        zero_pcts.append(pct)
    bars = ax.barh(range(nnodes), zero_pcts,
                   color=[district_colors[i] for i in range(nnodes)])
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("% Weeks with Zero Incidence")
    ax.set_title("(F) Zero-Incidence Weeks (post burn-in)")
    ax.set_xlim(0, 105)
    for i, v in enumerate(zero_pcts):
        ax.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=6)

    plt.tight_layout()
    outpath = outdir / "polio_20district_diagnostics.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nDiagnostic plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    # Run simulation
    model, scenario, beta_season_365 = run_model()

    # Print summary
    print_summary(model, scenario)

    # Compute weekly incidence
    weekly_df = compute_weekly_incidence(model, scenario)

    # Save weekly incidence CSV
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "weekly_incidence_20district.csv"
    weekly_df.to_csv(csv_path, index=False)
    print(f"\nWeekly incidence saved to {csv_path}")
    print(f"  Shape: {weekly_df.shape} ({weekly_df['week'].nunique()} weeks × "
          f"{weekly_df['district'].nunique()} districts)")

    # Zero-incidence analysis
    analyze_zero_incidence(weekly_df, scenario)

    # Diagnostic plots
    plot_diagnostics(model, scenario, beta_season_365, outdir)

    print("\nDone.")
