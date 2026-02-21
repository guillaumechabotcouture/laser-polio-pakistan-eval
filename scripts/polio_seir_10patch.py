#!/usr/bin/env python3
"""
Pakistan Polio SEIRV Model — 13 Districts with Separate Vaccine Tracking

Calibrated to MMWR/WER surveillance data (Vol. 72 No. 33, Vol. 73 No. 36).
Uses LASER framework for agent-based spatial SEIRV simulation.

Key features:
    - SEIRV compartments: Vaccinated (V) tracked separately from naturally
      Recovered (R). OPV provides mucosal immunity that wanes (~3 years),
      while natural infection provides stronger, longer-lasting immunity.
    - Routine Immunization at 6 weeks of age, 80% coverage (OPV)
    - SIAs every 6 months targeting 0-5 years, 90% coverage
    - Correlated missedness via per-agent reachable flag
    - Monsoon seasonal forcing (Jul-Oct peak, matching Pakistan polio seasonality)
    - Gravity-model spatial coupling between districts
    - Endemic corridor importation (Quetta, Peshawar, Bannu, Karachi)

Usage:
    /opt/anaconda3/bin/python3 scripts/polio_seir_10patch.py
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
from custom_components import (VaccinatedCompartment, PerPatchVaccinationSEIRV,
                               PatchImportation, VACCINATED)

# ============================================================================
# District Configuration (MMWR Vol. 72 No. 33, Vol. 73 No. 36)
#
# Vaccination uses CORRELATED MISSEDNESS model: each agent has a permanent
# "reachable" flag set at birth. Unreachable agents (zero-dose children in
# hard-to-reach, refusal, or inaccessible communities) are never vaccinated
# by either RI or SIA. This prevents the independent Bernoulli SIA model
# from overestimating cumulative coverage.
#
# RI/SIA coverage values represent MMWR-reported rates AMONG REACHABLE
# populations. The unreachable_frac determines the structural immunity gap.
#
# For endemic equilibrium: f_escape = unreachable + reachable × (1-RI) × (1-SIA)^k
# Must exceed 1/R0 = 0.167 for sustained transmission.
#
# 13 districts: original 10 + 3 southern KP tribal area districts
# (N Waziristan, S Waziristan, DI Khan) which account for >80% of KP cases.
# ============================================================================

DISTRICTS = pd.DataFrame([
    # name,          province,      pop,     lat,     lon,  ri_cov, sia_cov, unreach, init_immun, seed_I
    ("Quetta",       "Balochistan", 1100000, 30.1798, 66.9750, 0.55, 0.60, 0.35, 0.50, 3),
    ("Peshawar",     "KP",         2000000, 34.0150, 71.5249, 0.78, 0.72, 0.15, 0.80, 2),
    ("Karachi",      "Sindh",       3500000, 24.8607, 67.0011, 0.80, 0.74, 0.16, 0.82, 2),
    ("Lahore",       "Punjab",      2500000, 31.5204, 74.3587, 0.95, 0.90, 0.02, 0.97, 0),
    ("Islamabad",    "ICT",          500000, 33.6844, 73.0479, 0.92, 0.87, 0.03, 0.93, 0),
    ("Multan",       "Punjab",      1800000, 30.1575, 71.5249, 0.93, 0.88, 0.02, 0.96, 0),
    ("Faisalabad",   "Punjab",      2200000, 31.4504, 73.1350, 0.94, 0.89, 0.02, 0.96, 0),
    ("Rawalpindi",   "Punjab",      1800000, 33.5651, 73.0169, 0.93, 0.88, 0.02, 0.96, 0),
    ("Hyderabad",    "Sindh",       1700000, 25.3960, 68.3578, 0.78, 0.72, 0.16, 0.82, 0),
    ("Bannu",        "KP/South",     900000, 32.9860, 70.6027, 0.50, 0.55, 0.28, 0.55, 3),
    ("N_Waziristan", "KP/FATA",      550000, 32.3045, 69.8597, 0.30, 0.40, 0.40, 0.40, 3),
    ("S_Waziristan", "KP/FATA",      600000, 32.0711, 69.5310, 0.30, 0.40, 0.40, 0.40, 3),
    ("DI_Khan",      "KP/South",    1500000, 31.8320, 70.9015, 0.45, 0.50, 0.30, 0.52, 2),
], columns=["name", "province", "population", "lat", "lon",
            "ri_coverage", "sia_coverage", "unreachable_frac",
            "initial_immunity", "seed_I"])

# ============================================================================
# Parameters
# ============================================================================

PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": 20 * 365,       # 20-year simulation (first 10 = burn-in)
    "beta": 0.2143,           # R0 = 6, infectious period = 28d → beta = R0/D = 6/28
    "exp_shape": 3,           # Exposed ~3 days (gamma shape=3, scale=1)
    "exp_scale": 1.0,
    "inf_mean": 28,           # Infectious ~28 days
    "inf_sigma": 3,
    "cbr": 29.0,              # Crude birth rate per 1000/year
    "gravity_k": 0.005,
    "gravity_b": 0.5,
    "gravity_c": 1.5,
    "capacity_safety_factor": 2.0,
})

BURNIN_YEARS = 10  # Discard first 10 years of transient dynamics

# ============================================================================
# Monsoon Seasonal Forcing (Pakistan polio: Jul-Oct peak)
# ============================================================================

def build_monsoon_seasonality(nticks, nnodes):
    """Build 365-day seasonal profile with monsoon peak (Jul-Oct).

    Uses cosine-based forcing centered on day 245 (early September)
    with amplitude creating ~2x peak during monsoon vs baseline.
    """
    days = np.arange(365)
    # Peak around day 245 (Sep 2), monsoon season Jul-Oct (days 182-304)
    # Cosine forcing: 1 + amplitude * cos(2*pi*(day - peak_day)/365)
    peak_day = 245
    amplitude = 0.30  # ±30% modulation
    beta_season = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)

    # Tile across simulation
    season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
    return ValuesMap.from_timeseries(season_tiled, nnodes), beta_season


# ============================================================================
# Age Pyramid for Pakistan (young population, high birth rate)
# ============================================================================

def build_pakistan_age_pyramid():
    """Approximate Pakistan age distribution — heavily weighted toward youth.

    Follows the LASER skill convention: integer counts per year of age
    with exponential decay reflecting ~65-year life expectancy.
    """
    mean_life = 65.0
    ages = np.arange(100)  # 0 to 99 years
    # Stable age distribution: count at each age ~ exp(-age/mean_life)
    stable_age_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_age_dist = np.maximum(stable_age_dist, 1)  # Ensure non-zero
    return AliasedDistribution(stable_age_dist), stable_age_dist



# ============================================================================
# Build and Run Model
# ============================================================================

def build_scenario():
    """Build GeoDataFrame scenario from district configuration."""
    geometry = [Point(lon, lat) for lat, lon in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(DISTRICTS, geometry=geometry, crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))

    # Initial compartments based on per-district immunity
    # At t=0, all pre-existing immunity is counted as natural (R).
    # V starts at 0 — vaccination builds up the V compartment during simulation.
    scenario["I"] = scenario["seed_I"].astype(np.uint32)
    scenario["E"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["R"] = np.round(scenario["initial_immunity"] * scenario["population"]).astype(np.uint32)
    scenario["V"] = np.zeros(len(scenario), dtype=np.uint32)
    scenario["S"] = (scenario["population"] - scenario["E"] - scenario["I"] - scenario["R"]).astype(np.uint32)

    return scenario


def initialize_ages(model, pyramid, tick=0):
    """Set initial agent DOBs from age pyramid BEFORE model.run().

    Without this, all agents appear as newborns (dob=0), creating a massive
    artificial vaccination wave at day 270 when RI kicks in.
    """
    count = model.people.count
    # Sample ages in years, convert to days, set dob as negative values
    # Cap at 99 years to stay within survival estimator range
    ages_years = pyramid.sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, 99)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def run_model():
    """Build, configure, and run the Pakistan polio SEIR model."""
    scenario = build_scenario()
    nnodes = len(scenario)

    # Duration distributions
    expdurdist = dists.gamma(shape=PARAMS.exp_shape, scale=PARAMS.exp_scale)
    infdurdist = dists.normal(loc=PARAMS.inf_mean, scale=PARAMS.inf_sigma)

    # Seasonal forcing
    seasonality, beta_season_365 = build_monsoon_seasonality(PARAMS.nticks, nnodes)

    # Birth rates — uniform CBR across all districts
    # LASER expects raw CBR per 1000/year (BirthsByCBR divides by 1000 internally)
    birthrate_array = np.full((PARAMS.nticks, nnodes), PARAMS.cbr, dtype=np.float32)

    # Age pyramid and mortality rates
    pyramid, _stable_age_dist = build_pakistan_age_pyramid()
    # CDR ≈ 7/1000/year for Pakistan (life expectancy ~65, young population)
    # LASER expects raw CDR per 1000/year (MortalityByCDR divides by 1000 internally)
    cdr = 7.0
    mortality_array = np.full((PARAMS.nticks, nnodes), cdr, dtype=np.float32)

    # Build model with V compartment for vaccine tracking
    model = Model(scenario, PARAMS, birthrates=birthrate_array,
                  additional_states=["V"])

    # Gravity network (manual setup for custom normalization)
    from laser.core.migration import distance as haversine_distance
    lats = np.array(scenario.lat)
    lons = np.array(scenario.lon)
    dist_matrix = np.zeros((nnodes, nnodes))
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                dist_matrix[i, j] = haversine_distance(lats[i], lons[i], lats[j], lons[j])

    pops = np.array(scenario.population, dtype=np.float64)
    network = gravity(pops, dist_matrix, 1, 0, PARAMS.gravity_b, PARAMS.gravity_c)
    # Normalize
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * PARAMS.gravity_k
    network = row_normalizer(network, 0.2)  # Cap at 20% export per node
    model.network = network

    # Per-patch unreachable fractions (correlated missedness)
    unreachable_frac = np.array(scenario.unreachable_frac, dtype=np.float32)

    # OPV vaccine waning: gamma distribution, mean ~3 years (1095 days)
    # OPV mucosal immunity is shorter-lived than natural infection immunity
    vax_waning_dist = dists.gamma(shape=3, scale=365.0)  # mean = 3×365 = 1095 days

    # Endemic corridor patch indices
    name_to_idx = {n: i for i, n in enumerate(scenario.name)}
    endemic_names = ["Quetta", "Peshawar", "Karachi", "Hyderabad", "Bannu",
                     "N_Waziristan", "S_Waziristan", "DI_Khan"]
    endemic_patch_indices = [name_to_idx[n] for n in endemic_names if n in name_to_idx]

    # Assemble components in correct order.
    # SEIR.Susceptible must be created first (it adds people.nodeid, people.state).
    # VaccinatedCompartment and PerPatchVaccinationSEIRV depend on these properties.
    susceptible = SEIR.Susceptible(model)                           # 1. Count S (adds nodeid, state)
    exposed = SEIR.Exposed(model, expdurdist, infdurdist)          # 2. E→I transitions
    infectious = SEIR.Infectious(model, infdurdist)                 # 3. I→R transitions
    recovered = SEIR.Recovered(model)                               # 4. Count R (natural immunity)

    # Now create vaccination components (need people.nodeid from Susceptible)
    vax_compartment = VaccinatedCompartment(model, vax_waning_dist, wandurmin=365)  # 5. V waning
    vaccination = PerPatchVaccinationSEIRV(                         # 8. RI + SIA (S→V)
        model, vax_compartment, unreachable_frac,
        ri_coverage=0.80,       # 80% RI coverage (OPV at 6 weeks)
        sia_coverage=0.90,      # 90% SIA coverage
        ri_age=42,              # 6 weeks = 42 days
        sia_period=180,         # Every 6 months
        sia_max_age=1825,       # 0-5 years
    )

    model.components = [
        susceptible,                                                # 1. Count S
        exposed,                                                    # 2. E→I transitions
        infectious,                                                 # 3. I→R transitions
        recovered,                                                  # 4. Count R (natural)
        vax_compartment,                                            # 5. V waning (V→S)
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),  # 6. Births
        MortalityByCDR(model, mortalityrates=mortality_array, mappings=[  # 7. Deaths
            (SEIR.State.SUSCEPTIBLE.value, "S"),
            (SEIR.State.EXPOSED.value, "E"),
            (SEIR.State.INFECTIOUS.value, "I"),
            (SEIR.State.RECOVERED.value, "R"),
            (VACCINATED, "V"),
        ]),
        vaccination,                                                # 8. RI + SIA (S→V)
        PatchImportation(model, infdurdist, endemic_patch_indices,  # 9. Endemic importation
                         period=30, count=3),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),  # 10. Transmission
    ]

    # Initialize ages from pyramid AFTER components are set (BirthsByCBR creates dob)
    # but BEFORE model.run() — prevents all agents appearing as newborns
    initialize_ages(model, pyramid)

    nyears_total = PARAMS.nticks // 365
    print(f"Running {nyears_total}-year Pakistan polio SEIRV simulation ({BURNIN_YEARS}yr burn-in)...")
    print(f"  Districts: {len(scenario)}")
    print(f"  Total population: {scenario.population.sum():,}")
    print(f"  Vaccination: RI at 6wk (80%), SIA biannual 0-5yr (90%)")
    print(f"  OPV waning: ~3yr mean (gamma), Natural immunity: permanent")
    print(f"  Ticks: {PARAMS.nticks}")
    model.run("Pakistan Polio SEIRV")

    return model, scenario, beta_season_365


# ============================================================================
# 4-Panel Diagnostic Plot
# ============================================================================

def _get_node_arrays(model, ticks, node_idx):
    """Helper: return S, E, I, R, V, N arrays for a node at given ticks."""
    S = model.nodes.S[ticks, node_idx].astype(np.float64)
    E = model.nodes.E[ticks, node_idx].astype(np.float64)
    I = model.nodes.I[ticks, node_idx].astype(np.float64)
    R = model.nodes.R[ticks, node_idx].astype(np.float64)
    V = model.nodes.V[ticks, node_idx].astype(np.float64)
    N = np.maximum(S + E + I + R + V, 1.0)
    return S, E, I, R, V, N


def plot_diagnostics(model, scenario, beta_season_365):
    """Generate 6-panel overview figure."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    provinces = list(scenario.province)
    R0 = 6.0

    prov_colors = {
        "Balochistan": "#e74c3c", "KP": "#e67e22", "KP/South": "#f39c12",
        "Sindh": "#3498db", "Punjab": "#2ecc71", "ICT": "#9b59b6"
    }
    district_colors = [prov_colors.get(p, "gray") for p in provinces]

    incidence = model.nodes.newly_infected
    nyears = nticks // 365
    annual_inc = np.zeros((nyears, nnodes))
    for y in range(nyears):
        annual_inc[y] = incidence[y * 365:(y + 1) * 365].sum(axis=0)
    pops = np.array(scenario.population, dtype=np.float64)

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))

    # --- (A) Seasonal forcing + district RI coverage context ---
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, beta_season_365, "b-", lw=2)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
    ax.axvspan(182, 304, alpha=0.12, color="blue", label="Monsoon (Jul-Oct)")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Seasonal Multiplier")
    ax.set_title("(A) Monsoon Seasonal Forcing")
    ax.legend(fontsize=8)

    # --- (B) Annual incidence rate per million (all years, all districts) ---
    ax = axes[0, 1]
    years_x = np.arange(nyears) + 0.5
    for i in range(nnodes):
        rate = annual_inc[:, i] / pops[i] * 1e6
        if rate.max() > 1:
            ax.plot(years_x, rate, color=district_colors[i], alpha=0.8,
                    lw=1.2, label=names[i])
    ax.axvline(BURNIN_YEARS, color="black", ls=":", alpha=0.4, label="Burn-in")
    ax.axhspan(50, 500, alpha=0.08, color="green", label="Target range")
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual Incidence per Million")
    ax.set_title("(B) Incidence Rate per Million by District")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (C) R_eff over time for ALL districts ---
    ax = axes[0, 2]
    sample30 = np.arange(0, nticks, 30)
    time30 = sample30 / 365.0
    for i in range(nnodes):
        S, E, I, R, V, N = _get_node_arrays(model, sample30, i)
        ax.plot(time30, R0 * S / N, color=district_colors[i], alpha=0.7,
                lw=1.2, label=names[i])
    ax.axhline(1.0, color="black", ls="-", lw=2, alpha=0.4, label="R_eff = 1")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Year")
    ax.set_ylabel("R_eff = R0 × S/N")
    ax.set_title("(C) Effective Reproduction Number")
    ax.set_ylim(0, 4)
    ax.legend(fontsize=5, ncol=2, loc="upper right")

    # --- (D) S/N for ALL districts (not just 4) ---
    ax = axes[1, 0]
    for i in range(nnodes):
        S, E, I, R, V, N = _get_node_arrays(model, sample30, i)
        ls = "-" if provinces[i] in ("Balochistan", "KP", "KP/South", "Sindh") else "--"
        ax.plot(time30, S / N, color=district_colors[i], ls=ls, alpha=0.8,
                lw=1.2, label=names[i])
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.5, label="S* = 1/R0")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.3)
    ax.set_xlabel("Year")
    ax.set_ylabel("Susceptible Fraction (S/N)")
    ax.set_title("(D) Susceptible Fraction — All Districts")
    ax.legend(fontsize=5, ncol=2, loc="upper right")
    ax.set_ylim(0, 0.55)

    # --- (E) Infectious count (I) over time — shows outbreak pulses ---
    ax = axes[1, 1]
    sample7 = np.arange(0, nticks, 7)
    time7 = sample7 / 365.0
    for i in range(nnodes):
        I_t = model.nodes.I[sample7, i].astype(np.float64)
        if I_t.max() > 5:
            ax.plot(time7, I_t, color=district_colors[i], alpha=0.7,
                    lw=0.8, label=names[i])
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Year")
    ax.set_ylabel("Infectious Count (I)")
    ax.set_title("(E) Active Infections Over Time")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (F) Post-burn-in weekly incidence (zoomed, seasonality visible) ---
    ax = axes[1, 2]
    burnin_ticks = BURNIN_YEARS * 365
    post_inc = incidence[burnin_ticks:]
    post_weeks = post_inc.shape[0] // 7
    post_weekly = post_inc[:post_weeks * 7].reshape(post_weeks, 7, nnodes).sum(axis=1)
    weeks_x = np.arange(post_weeks) / 52.0 + BURNIN_YEARS
    for i in range(nnodes):
        col_max = post_weekly[:, i].max()
        if col_max > 0:
            ax.plot(weeks_x, post_weekly[:, i], color=district_colors[i],
                    alpha=0.8, lw=0.8, label=names[i])
    # Add monsoon shading for each post-burn-in year
    for yr in range(BURNIN_YEARS, nyears):
        ax.axvspan(yr + 182 / 365, yr + 304 / 365, alpha=0.04, color="blue")
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly Infections")
    ax.set_title("(F) Post-Burn-in Weekly Incidence (blue = monsoon)")
    ax.legend(fontsize=6, ncol=2)

    plt.tight_layout()
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "polio_seir_diagnostics.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Diagnostic plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# Summary Statistics
# ============================================================================

def print_summary(model, scenario):
    """Print key summary statistics for verification (post-burn-in only)."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)

    incidence = model.nodes.newly_infected
    nyears = nticks // 365
    analysis_years = nyears - BURNIN_YEARS
    burnin_ticks = BURNIN_YEARS * 365

    print("\n" + "=" * 70)
    print(f"SIMULATION SUMMARY (post-burn-in: years {BURNIN_YEARS}–{nyears})")
    print("=" * 70)

    # Annual incidence per district (post-burn-in)
    print(f"\n{'District':<14} {'Province':<14} {'Ann.Inf(med)':<14} {'Ann.Inf(post)':<14} {'S/N':<8} {'V/N':<8} {'R/N':<8}")
    print("-" * 84)

    total_infections = 0
    for i in range(nnodes):
        annual = np.array([incidence[y * 365:(y + 1) * 365, i].sum() for y in range(nyears)])
        post_burnin = annual[BURNIN_YEARS:]
        median_ann = np.median(post_burnin)
        total_post = post_burnin.sum()
        total_infections += total_post

        # Final compartment fractions
        final_tick = nticks - 1
        S_f = float(model.nodes.S[final_tick, i])
        R_f = float(model.nodes.R[final_tick, i])
        V_f = float(model.nodes.V[final_tick, i])
        N_f = float(S_f + model.nodes.E[final_tick, i] +
                     model.nodes.I[final_tick, i] + R_f + V_f)
        sn = S_f / max(N_f, 1)
        vn = V_f / max(N_f, 1)
        rn = R_f / max(N_f, 1)

        print(f"{names[i]:<14} {scenario.province.iloc[i]:<14} {median_ann:<14.0f} {total_post:<14.0f} {sn:<8.3f} {vn:<8.3f} {rn:<8.3f}")

    total_pop = scenario.population.sum()
    print(f"\n  Total population: {total_pop:,}")
    print(f"  Total infections (years {BURNIN_YEARS}–{nyears}): {total_infections:,.0f}")
    print(f"  Annual rate per million: {total_infections / analysis_years / total_pop * 1e6:.0f}")

    # Compartment integrity check (SEIRV)
    print("\nCompartment integrity check (SEIRV):")
    for tick in [0, nticks // 2, nticks - 1]:
        S = model.nodes.S[tick].sum()
        E = model.nodes.E[tick].sum()
        I = model.nodes.I[tick].sum()
        R = model.nodes.R[tick].sum()
        V = model.nodes.V[tick].sum()
        print(f"  Tick {tick:>5}: S={S:>8} E={E:>5} I={I:>5} R={R:>8} V={V:>8} Total={S+E+I+R+V:>8}")

    # Check for negative compartments
    any_negative = False
    for comp, arr in [("S", model.nodes.S), ("E", model.nodes.E),
                      ("I", model.nodes.I), ("R", model.nodes.R),
                      ("V", model.nodes.V)]:
        if np.any(arr[:nticks] < 0):
            print(f"  WARNING: Negative {comp} values detected!")
            any_negative = True
    if not any_negative:
        print("  All compartments non-negative: OK")

    # Population trajectory (yearly)
    print(f"\nPopulation trajectory (yearly):")
    for y in range(0, nyears + 1, 2):
        t = min(y * 365, nticks - 1)
        pop_t = int(model.nodes.S[t].sum() + model.nodes.E[t].sum() +
                     model.nodes.I[t].sum() + model.nodes.R[t].sum() +
                     model.nodes.V[t].sum())
        print(f"  Year {y:>2}: {pop_t:>10,}")
    print(f"  Agent capacity: {model.people.capacity:,}, Active: {model.people.count:,}")


# ============================================================================
# Main
# ============================================================================

def plot_dynamics(model, scenario):
    """Generate 6-panel deep-dive into disease dynamics."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    provinces = list(scenario.province)
    R0 = 6.0
    burnin_ticks = BURNIN_YEARS * 365

    prov_colors = {
        "Balochistan": "#e74c3c", "KP": "#e67e22", "KP/South": "#f39c12",
        "Sindh": "#3498db", "Punjab": "#2ecc71", "ICT": "#9b59b6"
    }
    district_colors = [prov_colors.get(p, "gray") for p in provinces]

    incidence = model.nodes.newly_infected
    nyears = nticks // 365
    pops = np.array(scenario.population, dtype=np.float64)

    fig, axes = plt.subplots(3, 2, figsize=(18, 16))

    # --- (A) Heatmap: annual incidence per million (log scale) ---
    ax = axes[0, 0]
    annual_inc = np.zeros((nyears, nnodes))
    for y in range(nyears):
        annual_inc[y] = incidence[y * 365:(y + 1) * 365].sum(axis=0)
    rate_per_m = annual_inc / pops[None, :] * 1e6
    log_rate = np.log10(rate_per_m + 1)
    im = ax.imshow(log_rate.T, aspect="auto", cmap="YlOrRd", origin="lower",
                   extent=[0, nyears, -0.5, nnodes - 0.5])
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Year")
    ax.axvline(BURNIN_YEARS, color="white", ls="--", lw=1.5)
    ax.set_title("(A) Annual Incidence Rate per Million (log10)")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log10(per million + 1)")

    # --- (B) Susceptible (S) and Infectious (I) for Peshawar (dual axis) ---
    ax = axes[0, 1]
    idx = names.index("Peshawar")
    sample7 = np.arange(0, nticks, 7)
    time7 = sample7 / 365.0
    S, E, I, R, V, N = _get_node_arrays(model, sample7, idx)

    ax.plot(time7, S / N, color="#3498db", lw=1, label="S/N (left)")
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.4, label="S* = 1/R0")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.3)
    ax.set_ylabel("Susceptible Fraction (S/N)", color="#3498db")
    ax.set_ylim(0, 0.35)
    ax.tick_params(axis="y", labelcolor="#3498db")

    ax2 = ax.twinx()
    ax2.plot(time7, I, color="#e74c3c", lw=0.7, alpha=0.8, label="I count (right)")
    ax2.set_ylabel("Infectious Count", color="#e74c3c")
    ax2.tick_params(axis="y", labelcolor="#e74c3c")

    ax.set_xlabel("Year")
    ax.set_title("(B) Peshawar — Susceptibles vs Infectious")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="upper right")

    # --- (C) Same dual-axis for Quetta ---
    ax = axes[1, 0]
    idx = names.index("Quetta")
    S, E, I, R, V, N = _get_node_arrays(model, sample7, idx)

    ax.plot(time7, S / N, color="#3498db", lw=1, label="S/N (left)")
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.4, label="S* = 1/R0")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.3)
    ax.set_ylabel("Susceptible Fraction (S/N)", color="#3498db")
    ax.set_ylim(0, 0.55)
    ax.tick_params(axis="y", labelcolor="#3498db")

    ax2 = ax.twinx()
    ax2.plot(time7, I, color="#e74c3c", lw=0.7, alpha=0.8, label="I count (right)")
    ax2.set_ylabel("Infectious Count", color="#e74c3c")
    ax2.tick_params(axis="y", labelcolor="#e74c3c")

    ax.set_xlabel("Year")
    ax.set_title("(C) Quetta — Susceptibles vs Infectious")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="upper right")

    # --- (D) Post-burn-in monthly incidence per million, all districts ---
    ax = axes[1, 1]
    post_inc = incidence[burnin_ticks:]
    nmonths = post_inc.shape[0] // 30
    monthly = post_inc[:nmonths * 30].reshape(nmonths, 30, nnodes).sum(axis=1)
    months_x = np.arange(nmonths) / 12.0 + BURNIN_YEARS
    for i in range(nnodes):
        monthly_rate = monthly[:, i] / pops[i] * 1e6
        if monthly_rate.max() > 0.5:
            ax.plot(months_x, monthly_rate, color=district_colors[i],
                    alpha=0.8, lw=0.9, label=names[i])
    # Monsoon shading
    for yr in range(BURNIN_YEARS, nyears):
        ax.axvspan(yr + 182 / 365, yr + 304 / 365, alpha=0.05, color="blue")
    ax.set_xlabel("Year")
    ax.set_ylabel("Monthly Incidence per Million")
    ax.set_title("(D) Post-Burn-in Monthly Incidence Rate (blue = monsoon)")
    ax.legend(fontsize=7, ncol=2)

    # --- (E) Cumulative infections per million by district ---
    ax = axes[2, 0]
    for i in range(nnodes):
        cum_rate = np.cumsum(annual_inc[:, i]) / pops[i] * 1e6
        ax.plot(np.arange(nyears) + 0.5, cum_rate, color=district_colors[i],
                lw=1.5, alpha=0.8, label=names[i])
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.3)
    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative Infections per Million")
    ax.set_title("(E) Cumulative Incidence per Million")
    ax.legend(fontsize=6, ncol=2)

    # --- (F) Force of infection (FOI) over time for endemic districts ---
    ax = axes[2, 1]
    foi = model.nodes.forces  # shape: (nticks, nnodes)
    # Smooth with 30-day rolling mean
    window = 30
    endemic_show = ["Quetta", "Peshawar", "Karachi", "Bannu"]
    for dname in endemic_show:
        idx = names.index(dname)
        foi_d = foi[:, idx].astype(np.float64)
        # Rolling mean
        kernel = np.ones(window) / window
        foi_smooth = np.convolve(foi_d, kernel, mode="valid")
        time_smooth = np.arange(len(foi_smooth)) / 365.0
        ax.plot(time_smooth, foi_smooth, color=district_colors[idx],
                lw=0.8, alpha=0.8, label=dname)
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.3)
    ax.set_xlabel("Year")
    ax.set_ylabel("Force of Infection (30-day avg)")
    ax.set_title("(F) Force of Infection — Endemic Districts")
    ax.legend(fontsize=8)

    plt.tight_layout()
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outpath = outdir / "polio_dynamics_detail.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Dynamics detail plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# MMWR/GPEI Comparison: Model vs Observed WPV1 Paralytic Cases
# ============================================================================

# GPEI/MMWR observed WPV1 paralytic cases by province (2020-2025)
# Sources:
#   CDC MMWR Vol. 69 No. 46 (2020), Vol. 71 No. 42 (2022),
#   Vol. 72 No. 33 (2023), Vol. 73 No. 36 (2024)
#   Pakistan NEOC endpolio.com.pk (2024 year-end)
#   WHO 43rd Polio IHR Emergency Committee (Nov 2025)
#   Arab News (2025 year-end tally)
OBSERVED_WPV1 = pd.DataFrame({
    "year": [2020, 2021, 2022, 2023, 2024, 2025],
    "Balochistan": [26, 1, 0, 0, 27, 0],
    "KP":          [22, 0, 20, 4, 22, 20],
    "Sindh":       [22, 0, 0, 2, 23, 9],
    "Punjab":      [14, 0, 0, 0, 1, 1],
    "ICT":         [0, 0, 0, 0, 1, 0],
    "Total":       [84, 1, 20, 6, 74, 31],
})

# WPV1 paralysis-to-infection ratio (CDC: ~1:200 for Type 1)
PARALYSIS_RATIO = 1.0 / 200.0

# Province populations (approximate 2024, millions)
PROVINCE_POP = {
    "Balochistan": 14.0e6,
    "KP": 40.0e6,
    "Sindh": 50.0e6,
    "Punjab": 120.0e6,
    "ICT": 2.0e6,
}

# Mapping: model districts → provinces (for MMWR comparison)
DISTRICT_TO_PROVINCE = {
    "Quetta": "Balochistan",
    "Peshawar": "KP", "Bannu": "KP",
    "N_Waziristan": "KP", "S_Waziristan": "KP", "DI_Khan": "KP",
    "Karachi": "Sindh", "Hyderabad": "Sindh",
    "Lahore": "Punjab", "Multan": "Punjab",
    "Faisalabad": "Punjab", "Rawalpindi": "Punjab",
    "Islamabad": "ICT",
}


def compare_with_mmwr(model, scenario):
    """Compare model output with MMWR/GPEI observed WPV1 paralytic cases.

    Translates model infections → expected paralytic cases using the
    WPV1 paralysis-to-infection ratio (~1:200).
    """
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    provinces = list(scenario.province)
    nyears = nticks // 365
    analysis_years = nyears - BURNIN_YEARS

    incidence = model.nodes.newly_infected

    # --- Aggregate model infections by province ---
    model_by_province = {}
    model_annual_by_district = {}
    for i in range(nnodes):
        annual = np.array([incidence[y * 365:(y + 1) * 365, i].sum()
                           for y in range(BURNIN_YEARS, nyears)])
        model_annual_by_district[names[i]] = annual
        prov = DISTRICT_TO_PROVINCE[names[i]]
        if prov not in model_by_province:
            model_by_province[prov] = np.zeros(analysis_years)
        model_by_province[prov] += annual

    # --- Summary table ---
    print("\n" + "=" * 80)
    print("MODEL vs MMWR/GPEI COMPARISON")
    print("=" * 80)
    print(f"\nWPV1 paralysis-to-infection ratio: 1:{int(1/PARALYSIS_RATIO)}")
    print(f"Model analysis period: years {BURNIN_YEARS}–{nyears} ({analysis_years} years)")
    print(f"Observed data: 2020–2025 (6 years)\n")

    # Province-level comparison
    print(f"{'Province':<14} {'Model Inf/yr':<14} {'Model Para/yr':<14} "
          f"{'Obs Para/yr':<14} {'Obs 2024':<10} {'Ratio':<8}")
    print("-" * 74)

    obs_avg = OBSERVED_WPV1[["Balochistan", "KP", "Sindh", "Punjab", "ICT"]].mean()
    total_model_inf = 0
    total_model_para = 0
    total_obs_avg = 0

    for prov in ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]:
        if prov in model_by_province:
            model_inf_yr = model_by_province[prov].mean()
            model_para_yr = model_inf_yr * PARALYSIS_RATIO
            obs_avg_yr = obs_avg[prov]
            obs_2024 = OBSERVED_WPV1[OBSERVED_WPV1.year == 2024][prov].values[0]
            ratio = model_para_yr / max(obs_avg_yr, 0.1)
            total_model_inf += model_inf_yr
            total_model_para += model_para_yr
            total_obs_avg += obs_avg_yr
            print(f"{prov:<14} {model_inf_yr:<14.0f} {model_para_yr:<14.1f} "
                  f"{obs_avg_yr:<14.1f} {obs_2024:<10} {ratio:<8.2f}")
        else:
            print(f"{prov:<14} {'N/A':<14} {'N/A':<14} {obs_avg[prov]:<14.1f}")

    total_obs_2024 = OBSERVED_WPV1[OBSERVED_WPV1.year == 2024]["Total"].values[0]
    total_ratio = total_model_para / max(total_obs_avg, 0.1)
    print("-" * 74)
    print(f"{'TOTAL':<14} {total_model_inf:<14.0f} {total_model_para:<14.1f} "
          f"{total_obs_avg:<14.1f} {total_obs_2024:<10} {total_ratio:<8.2f}")

    # --- Per-capita rate comparison ---
    print(f"\n{'Province':<14} {'Model rate*':<14} {'Obs rate*':<14} {'Note'}")
    print("-" * 60)
    print("  * Paralytic cases per million per year")
    for prov in ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]:
        if prov in model_by_province:
            # Model: use initial district populations (model represents subset)
            model_prov_pop = sum(scenario.population.iloc[i]
                                for i in range(nnodes)
                                if DISTRICT_TO_PROVINCE[names[i]] == prov)
            model_rate = model_by_province[prov].mean() * PARALYSIS_RATIO / model_prov_pop * 1e6
            obs_rate = obs_avg[prov] / PROVINCE_POP[prov] * 1e6
            print(f"{prov:<14} {model_rate:<14.2f} {obs_rate:<14.2f}")

    # --- District-level detail ---
    print(f"\n{'District':<14} {'Infections/yr':<14} {'Paralytic/yr':<14} {'Interquartile range'}")
    print("-" * 56)
    for i in range(nnodes):
        arr = model_annual_by_district[names[i]]
        mean_inf = arr.mean()
        para = mean_inf * PARALYSIS_RATIO
        q25 = np.percentile(arr, 25) * PARALYSIS_RATIO
        q75 = np.percentile(arr, 75) * PARALYSIS_RATIO
        if mean_inf > 0.5:
            print(f"{names[i]:<14} {mean_inf:<14.0f} {para:<14.1f} [{q25:.1f} – {q75:.1f}]")
        else:
            print(f"{names[i]:<14} {'<1':<14} {'<0.01':<14} —")

    # --- Interpretation ---
    print(f"\nKey findings:")
    print(f"  Model predicts ~{total_model_para:.0f} paralytic cases/year nationally")
    print(f"  MMWR/GPEI average: ~{total_obs_avg:.0f} paralytic cases/year (2020-2025)")
    print(f"  Model/Observed ratio: {total_ratio:.2f}")
    if total_ratio < 0.5:
        print(f"  → Model UNDERESTIMATES by ~{1/total_ratio:.1f}x")
    elif total_ratio > 2.0:
        print(f"  → Model OVERESTIMATES by ~{total_ratio:.1f}x")
    else:
        print(f"  → Model is within 2x of observed data")
    print(f"\n  Note: Observed data shows extreme inter-annual variability")
    print(f"        (1 case in 2021 vs 84 in 2020), reflecting stochastic")
    print(f"        dynamics and campaign effectiveness variation.")
    print(f"        Model represents 10 districts ({scenario.population.sum()/1e6:.0f}M),")
    print(f"        not all ~160 districts in Pakistan (~240M).")


def plot_mmwr_comparison(model, scenario):
    """Generate comparison figure: model predictions vs MMWR observed data."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    nyears = nticks // 365
    analysis_years = nyears - BURNIN_YEARS

    incidence = model.nodes.newly_infected

    # Aggregate model by province
    model_by_province = {}
    for i in range(nnodes):
        annual = np.array([incidence[y * 365:(y + 1) * 365, i].sum()
                           for y in range(BURNIN_YEARS, nyears)])
        prov = DISTRICT_TO_PROVINCE[names[i]]
        if prov not in model_by_province:
            model_by_province[prov] = np.zeros(analysis_years)
        model_by_province[prov] += annual

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    prov_colors = {
        "Balochistan": "#e74c3c", "KP": "#e67e22",
        "Sindh": "#3498db", "Punjab": "#2ecc71", "ICT": "#9b59b6"
    }

    # --- (A) Annual paralytic cases: model distribution vs observed ---
    ax = axes[0, 0]
    obs_years = OBSERVED_WPV1["year"].values
    obs_totals = OBSERVED_WPV1["Total"].values
    ax.bar(obs_years, obs_totals, width=0.35, color="black", alpha=0.7,
           label="Observed (GPEI/MMWR)", align="center")

    # Model: convert infections to paralytic, show median and IQR
    model_total = np.zeros(analysis_years)
    for prov, arr in model_by_province.items():
        model_total += arr
    model_para = model_total * PARALYSIS_RATIO
    model_years = np.arange(BURNIN_YEARS, nyears)

    # Show model as horizontal band (median ± IQR) since model years don't map to calendar years
    med = np.median(model_para)
    q25, q75 = np.percentile(model_para, [25, 75])
    ax.axhspan(q25, q75, alpha=0.2, color="steelblue", label=f"Model IQR [{q25:.0f}–{q75:.0f}]")
    ax.axhline(med, color="steelblue", ls="--", lw=2,
               label=f"Model median: {med:.0f} cases/yr")
    ax.set_xlabel("Year")
    ax.set_ylabel("WPV1 Paralytic Cases")
    ax.set_title("(A) Annual Paralytic Cases: Model vs Observed")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, max(obs_totals.max(), q75) * 1.3)

    # --- (B) Province-level comparison (bar chart) ---
    ax = axes[0, 1]
    prov_order = ["Balochistan", "KP", "Sindh", "Punjab", "ICT"]
    x = np.arange(len(prov_order))
    width = 0.35

    obs_avg_vals = [OBSERVED_WPV1[p].mean() for p in prov_order]
    model_avg_vals = []
    for p in prov_order:
        if p in model_by_province:
            model_avg_vals.append(model_by_province[p].mean() * PARALYSIS_RATIO)
        else:
            model_avg_vals.append(0)

    bars1 = ax.bar(x - width / 2, obs_avg_vals, width, color="black", alpha=0.7,
                   label="Observed avg (2020-2025)")
    bars2 = ax.bar(x + width / 2, model_avg_vals, width,
                   color=[prov_colors[p] for p in prov_order], alpha=0.7,
                   label="Model (÷200)")
    ax.set_xticks(x)
    ax.set_xticklabels(prov_order, rotation=30, ha="right")
    ax.set_ylabel("Paralytic Cases / Year")
    ax.set_title("(B) Province-Level: Model vs Observed (avg/yr)")
    ax.legend(fontsize=8)

    # --- (C) Per-capita paralytic rate comparison ---
    ax = axes[1, 0]
    obs_rates = [OBSERVED_WPV1[p].mean() / PROVINCE_POP[p] * 1e6 for p in prov_order]
    model_rates = []
    for p in prov_order:
        if p in model_by_province:
            model_prov_pop = sum(scenario.population.iloc[i]
                                for i in range(nnodes)
                                if DISTRICT_TO_PROVINCE[names[i]] == p)
            model_rates.append(model_by_province[p].mean() * PARALYSIS_RATIO
                               / model_prov_pop * 1e6)
        else:
            model_rates.append(0)

    bars1 = ax.bar(x - width / 2, obs_rates, width, color="black", alpha=0.7,
                   label="Observed rate")
    bars2 = ax.bar(x + width / 2, model_rates, width,
                   color=[prov_colors[p] for p in prov_order], alpha=0.7,
                   label="Model rate")
    ax.set_xticks(x)
    ax.set_xticklabels(prov_order, rotation=30, ha="right")
    ax.set_ylabel("Paralytic Cases per Million per Year")
    ax.set_title("(C) Per-Capita Paralytic Rate by Province")
    ax.legend(fontsize=8)

    # --- (D) Model annual time series by province (paralytic scale) ---
    ax = axes[1, 1]
    model_years_x = np.arange(analysis_years) + BURNIN_YEARS + 0.5
    for prov in prov_order:
        if prov in model_by_province:
            para_ts = model_by_province[prov] * PARALYSIS_RATIO
            if para_ts.max() > 0.05:
                ax.plot(model_years_x, para_ts, color=prov_colors[prov],
                        lw=1.5, alpha=0.8, label=prov)
    ax.set_xlabel("Simulation Year")
    ax.set_ylabel("Expected Paralytic Cases (infections ÷ 200)")
    ax.set_title("(D) Model Paralytic Case Time Series by Province")
    ax.legend(fontsize=8)

    plt.tight_layout()
    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outpath = outdir / "polio_mmwr_comparison.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"MMWR comparison plot saved to {outpath}")
    plt.close(fig)


if __name__ == "__main__":
    model, scenario, beta_season_365 = run_model()
    print_summary(model, scenario)
    compare_with_mmwr(model, scenario)
    plot_diagnostics(model, scenario, beta_season_365)
    plot_dynamics(model, scenario)
    plot_mmwr_comparison(model, scenario)
