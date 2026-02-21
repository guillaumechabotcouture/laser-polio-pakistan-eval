#!/usr/bin/env python3
"""
10-Patch Polio SEIR Model — Gravity Migration + Monsoon Seasonal Forcing

Extends the basic 10-patch polio SEIR model with:
  1. Gravity-model migration network (b=0.5, c=1.5) between patches
     arranged in a line 50km apart
  2. Monsoon seasonal forcing: 1.3x peak Jul-Oct, 0.7x trough Dec-Mar
  3. Row-normalization capping FOI export at 15% per patch

Parameters:
    R0 ~ 6 (beta = R0/D_inf = 6/28 ≈ 0.2143)
    Latent period: ~3 days (gamma, shape=3, scale=1)
    Infectious period: ~28 days (normal, mean=28, sigma=3)
    Crude birth rate: 29/1000/year
    Crude death rate: 7/1000/year

Usage:
    /opt/anaconda3/bin/python3 scripts/polio_gravity_seasonal.py
"""

import sys
import numpy as np
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
from custom_components import PatchImportation

# ============================================================================
# Configuration
# ============================================================================

N_PATCHES = 10
POP_PER_PATCH = 100_000
SPACING_KM = 50.0           # Inter-patch spacing along the line
N_SEED_INFECTIOUS = 5        # Initial infectious per patch
INITIAL_IMMUNE_FRAC = 0.95   # 95% initially recovered (immune)

PARAMS = PropertySet({
    "prng_seed": 42,
    "nticks": 10 * 365,            # 10-year simulation
    "beta": 6.0 / 28.0,            # R0=6, D_inf=28d → beta ≈ 0.2143
    "exp_shape": 3,                 # Latent period ~3d (gamma shape=3, scale=1)
    "exp_scale": 1.0,
    "inf_mean": 28,                 # Infectious period ~28d
    "inf_sigma": 3,
    "cbr": 29.0,                    # Crude birth rate per 1000/year
    "cdr": 7.0,                     # Crude death rate per 1000/year
    "gravity_k": 0.005,             # Spatial coupling strength
    "gravity_b": 0.5,               # Destination population exponent
    "gravity_c": 1.5,               # Distance decay exponent
    "max_export_frac": 0.15,        # Row-normalization cap
    "capacity_safety_factor": 2.5,
})

PATCH_NAMES = [f"Patch_{i+1:02d}" for i in range(N_PATCHES)]


# ============================================================================
# 1. Distance Matrix — Patches in a line, 50km apart
# ============================================================================

def build_linear_distance_matrix(n, spacing_km=50.0):
    """Build pairwise distance matrix for n patches along a line.

    Patch i is at position i * spacing_km.
    Distance(i, j) = |i - j| * spacing_km.
    """
    positions = np.arange(n, dtype=np.float64) * spacing_km
    dist_matrix = np.abs(positions[:, None] - positions[None, :])
    np.fill_diagonal(dist_matrix, 0.0)
    return dist_matrix


# ============================================================================
# 2. Gravity-Model Migration Network
# ============================================================================

def build_gravity_network(populations, dist_matrix, k, b, c, max_export_frac):
    """Build gravity-model migration network with row normalization.

    Gravity law: M_{i,j} = k * p_j^b / d_{ij}^c  (source exponent a=0)

    Steps:
      1. Compute raw gravity weights via laser.core.migration.gravity()
      2. Scale so k represents average export fraction
      3. Row-normalize: cap each patch's total FOI export at max_export_frac

    Args:
        populations: 1D array of patch populations
        dist_matrix: 2D pairwise distance matrix (km), zeros on diagonal
        k: Overall coupling constant (average export fraction before capping)
        b: Destination population exponent (0.5)
        c: Distance decay exponent (1.5)
        max_export_frac: Max fraction of FOI any patch can export (0.15)

    Returns:
        network: Row-normalized nnodes x nnodes matrix
    """
    pops = np.asarray(populations, dtype=np.float64)

    # gravity(populations, distances, k, a, b, c)
    # a=0: source population does not affect outward flow
    network = gravity(pops, dist_matrix, 1.0, 0, b, c)

    # Scale so k represents average export fraction
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * k

    # Row-normalize: cap each row sum at max_export_frac
    network = row_normalizer(network, max_export_frac)

    return network


# ============================================================================
# 3. Monsoon Seasonal Forcing
# ============================================================================

def build_monsoon_seasonality(nticks, nnodes):
    """Build 365-day seasonal profile for Pakistan monsoon polio dynamics.

    Specification:
      - Peak Jul-Oct (days 182-304): 1.3x baseline
      - Trough Dec-Mar (days 335-89): 0.7x baseline
      - Smooth cosine transitions between seasons

    Implementation uses cosine forcing centered on day 245 (Sep 2):
      season(t) = 1.0 + 0.3 * cos(2π(t - 245)/365)
        → maximum at day 245: 1.0 + 0.3 = 1.3  (mid-monsoon)
        → minimum at day 62:  1.0 - 0.3 = 0.7  (mid-dry season)
        → mean = 1.0 (cosine integrates to zero over full period)

    Returns:
        seasonality: ValuesMap for SEIR.Transmission
        season_365: 365-day profile array (for plotting)
    """
    days = np.arange(365)
    peak_day = 245   # Sep 2 — center of Jul-Oct monsoon window
    amplitude = 0.3  # ±30% → 1.3x peak, 0.7x trough

    season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)

    # Verify profile properties
    assert abs(season_365.mean() - 1.0) < 1e-3, \
        f"Seasonal profile must average to ~1.0, got {season_365.mean():.4f}"
    assert abs(season_365.max() - 1.3) < 0.01, \
        f"Peak should be ~1.3, got {season_365.max():.4f}"
    assert abs(season_365.min() - 0.7) < 0.01, \
        f"Trough should be ~0.7, got {season_365.min():.4f}"

    # Tile across simulation duration
    season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
    seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

    return seasonality, season_365


# ============================================================================
# Build and Run
# ============================================================================

def build_age_pyramid():
    """Approximate Pakistan age distribution (~65yr life expectancy)."""
    mean_life = 65.0
    ages = np.arange(100)
    stable_age_dist = np.array(1000 * np.exp(-ages / mean_life), dtype=np.int64)
    stable_age_dist = np.maximum(stable_age_dist, 1)
    return AliasedDistribution(stable_age_dist)


def build_scenario():
    """Create GeoDataFrame with N_PATCHES equal-population patches in a line.

    Geometry uses synthetic coordinates along longitude 67-72°E at lat 30°N.
    Actual inter-patch distances come from the explicit 50km-spaced matrix.
    """
    # Synthetic coordinates — only used for GeoDataFrame geometry requirement
    base_lon = 67.0
    lats = np.full(N_PATCHES, 30.0)
    lons = base_lon + np.arange(N_PATCHES) * 0.5

    data = {
        "nodeid": range(N_PATCHES),
        "name": PATCH_NAMES,
        "population": [POP_PER_PATCH] * N_PATCHES,
        "lat": lats,
        "lon": lons,
    }
    geometry = [Point(lon, lat) for lat, lon in zip(lats, lons)]
    scenario = gpd.GeoDataFrame(data, geometry=geometry, crs="EPSG:4326")

    # Initial compartments: 95% recovered, 5 infectious, rest susceptible
    scenario["I"] = np.full(N_PATCHES, N_SEED_INFECTIOUS, dtype=np.uint32)
    scenario["E"] = np.zeros(N_PATCHES, dtype=np.uint32)
    scenario["R"] = np.full(
        N_PATCHES, int(INITIAL_IMMUNE_FRAC * POP_PER_PATCH), dtype=np.uint32
    )
    scenario["S"] = (
        scenario["population"] - scenario["E"]
        - scenario["I"] - scenario["R"]
    ).astype(np.uint32)

    return scenario


def initialize_ages(model, pyramid):
    """Assign realistic ages from pyramid so agents aren't all born at tick 0."""
    count = model.people.count
    ages_years = pyramid.sample(count=count, dtype=np.int32)
    ages_years = np.minimum(ages_years, 99)
    ages_days = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def run_model():
    """Build, configure, and run the extended 10-patch polio SEIR model."""
    scenario = build_scenario()
    nnodes = len(scenario)
    nticks = PARAMS.nticks

    # Verify initial compartments
    for _, row in scenario.iterrows():
        assert row["S"] + row["E"] + row["I"] + row["R"] == row["population"], \
            f"Compartment mismatch for {row['name']}"

    # Duration distributions
    expdurdist = dists.gamma(shape=PARAMS.exp_shape, scale=PARAMS.exp_scale)
    infdurdist = dists.normal(loc=PARAMS.inf_mean, scale=PARAMS.inf_sigma)

    # --- Seasonal forcing ---
    seasonality, season_365 = build_monsoon_seasonality(nticks, nnodes)

    # --- Birth/death rates (per-1000/year — LASER divides by 1000 internally) ---
    assert PARAMS.cbr >= 1 and PARAMS.cbr <= 60, \
        f"CBR must be per-1000/year, got {PARAMS.cbr}"
    birthrate_array = np.full((nticks, nnodes), PARAMS.cbr, dtype=np.float32)
    deathrate_array = np.full((nticks, nnodes), PARAMS.cdr, dtype=np.float32)

    # Age pyramid
    pyramid = build_age_pyramid()

    # --- Build LASER model ---
    model = Model(scenario, PARAMS, birthrates=birthrate_array)

    # --- Gravity migration network (50km-spaced linear patches) ---
    dist_matrix = build_linear_distance_matrix(nnodes, spacing_km=SPACING_KM)
    pops = np.array(scenario.population, dtype=np.float64)
    network = build_gravity_network(
        pops, dist_matrix,
        k=PARAMS.gravity_k,
        b=PARAMS.gravity_b,
        c=PARAMS.gravity_c,
        max_export_frac=PARAMS.max_export_frac,
    )
    model.network = network

    # Print network summary
    row_sums = network.sum(axis=1)
    print("Gravity network:")
    print(f"  Layout: {nnodes} patches in a line, {SPACING_KM}km apart")
    print(f"  Distance decay c={PARAMS.gravity_c}, "
          f"dest pop exponent b={PARAMS.gravity_b}")
    print(f"  Row sums (FOI export fraction): "
          f"min={row_sums.min():.4f}, max={row_sums.max():.4f}, "
          f"mean={row_sums.mean():.4f}")
    print(f"  15% export cap: "
          f"{'SATISFIED' if row_sums.max() <= 0.1501 else 'VIOLATED'}")

    # Print distance and network matrices
    print(f"\n  Distance matrix (km):")
    for i in range(nnodes):
        dists_row = " ".join(f"{dist_matrix[i,j]:6.0f}" for j in range(nnodes))
        print(f"    {PATCH_NAMES[i]:>10s}  {dists_row}")

    print(f"\n  Network matrix (migration weights):")
    header = "              " + " ".join(f"{n:>10s}" for n in PATCH_NAMES)
    print(f"    {header}")
    for i in range(nnodes):
        vals = " ".join(f"{network[i,j]:10.6f}" for j in range(nnodes))
        print(f"    {PATCH_NAMES[i]:>10s}  {vals}  sum={row_sums[i]:.4f}")

    # Endemic importation in patches 0, 4, 9 to prevent stochastic fadeout
    endemic_patches = [0, 4, 9]

    # --- Assemble components (order matters for S+E+I+R=N invariant) ---
    model.components = [
        SEIR.Susceptible(model),
        SEIR.Exposed(model, expdurdist, infdurdist),
        SEIR.Infectious(model, infdurdist),
        SEIR.Recovered(model),
        BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
        MortalityByCDR(model, mortalityrates=deathrate_array),
        PatchImportation(model, infdurdist, endemic_patches,
                         period=30, count=2),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    ]

    # Initialize ages before run
    initialize_ages(model, pyramid)

    total_pop = scenario.population.sum()
    R0 = PARAMS.beta * PARAMS.inf_mean
    print(f"\nRunning {nticks // 365}-year spatial SEIR simulation...")
    print(f"  Patches: {nnodes}, Population/patch: {POP_PER_PATCH:,}")
    print(f"  Total population: {total_pop:,}")
    print(f"  R0 = beta x D_inf = {PARAMS.beta:.4f} x {PARAMS.inf_mean} = {R0:.1f}")
    print(f"  Seasonal forcing: 1.3x peak (Jul-Oct), 0.7x trough (Dec-Mar)")
    print(f"  Gravity: k={PARAMS.gravity_k}, b={PARAMS.gravity_b}, "
          f"c={PARAMS.gravity_c}, 15% export cap")
    print(f"  Importation: 2 infections in 3 patches every 30d")

    model.run("Polio SEIR — Gravity + Monsoon")

    return model, scenario, season_365, network, dist_matrix


# ============================================================================
# Summary Statistics
# ============================================================================

def print_summary(model, scenario, network):
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

    # Network check
    row_sums = network.sum(axis=1)
    print(f"\nGravity network: b={PARAMS.gravity_b}, c={PARAMS.gravity_c}")
    print(f"  Export fraction: min={row_sums.min():.4f}, "
          f"max={row_sums.max():.4f}")
    print(f"  15% cap enforced: {row_sums.max() <= 0.1501}")

    # Per-patch summary
    print(f"\n{'Patch':<12} {'Avg Inf/yr':>12} {'Total Inf':>12} {'S/N(final)':>10}")
    print("-" * 48)

    total_infections = 0
    for i in range(nnodes):
        annual = np.array([
            incidence[y * 365:(y + 1) * 365, i].sum()
            for y in range(nyears)
        ])
        mean_ann = annual.mean()
        total = annual.sum()
        total_infections += total

        t = nticks - 1
        S_f = float(model.nodes.S[t, i])
        N_f = float(model.nodes.S[t, i] + model.nodes.E[t, i]
                     + model.nodes.I[t, i] + model.nodes.R[t, i])
        sn = S_f / max(N_f, 1)
        print(f"{names[i]:<12} {mean_ann:>12.0f} {total:>12.0f} {sn:>10.3f}")

    total_pop = scenario.population.sum()
    print(f"\n  Total infections over {nyears} years: {total_infections:,.0f}")
    print(f"  Avg annual infections: {total_infections / nyears:,.0f}")
    print(f"  Annual rate per million: "
          f"{total_infections / nyears / total_pop * 1e6:.0f}")

    # Population trajectory
    print(f"\nPopulation trajectory:")
    for y in range(0, nyears + 1, 2):
        t = min(y * 365, nticks - 1)
        pop_t = int(model.nodes.S[t].sum() + model.nodes.E[t].sum()
                     + model.nodes.I[t].sum() + model.nodes.R[t].sum())
        print(f"  Year {y:>2}: {pop_t:>10,}")

    # Compartment integrity
    print(f"\nCompartment integrity:")
    for tick in [0, nticks // 2, nticks - 1]:
        S = model.nodes.S[tick].sum()
        E = model.nodes.E[tick].sum()
        I = model.nodes.I[tick].sum()
        R = model.nodes.R[tick].sum()
        print(f"  Tick {tick:>5}: S={S:>8} E={E:>5} I={I:>5} "
              f"R={R:>8} N={S+E+I+R:>8}")

    # R_eff in second half
    sample = np.arange(nticks // 2, nticks, 30)
    S_all = model.nodes.S[sample].sum(axis=1).astype(np.float64)
    N_all = (model.nodes.S[sample].sum(axis=1)
             + model.nodes.E[sample].sum(axis=1)
             + model.nodes.I[sample].sum(axis=1)
             + model.nodes.R[sample].sum(axis=1)).astype(np.float64)
    reff = R0 * S_all / np.maximum(N_all, 1)
    print(f"\n  R_eff (years 5-10): mean={reff.mean():.2f}, "
          f"range=[{reff.min():.2f}, {reff.max():.2f}]")


# ============================================================================
# Diagnostic Plots
# ============================================================================

def plot_diagnostics(model, scenario, season_365, network, dist_matrix):
    """Generate 6-panel diagnostic figure."""
    nticks = PARAMS.nticks
    nnodes = len(scenario)
    names = list(scenario.name)
    R0 = PARAMS.beta * PARAMS.inf_mean
    incidence = model.nodes.newly_infected
    nyears = nticks // 365

    cmap = plt.cm.tab10
    colors = [cmap(i) for i in range(nnodes)]

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))

    # --- (A) Seasonal forcing profile ---
    ax = axes[0, 0]
    days = np.arange(365)
    ax.plot(days, season_365, "b-", lw=2)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.5)
    ax.axvspan(182, 304, alpha=0.12, color="blue", label="Monsoon (Jul-Oct)")
    ax.axvspan(335, 365, alpha=0.12, color="orange", label="Dry (Dec-Mar)")
    ax.axvspan(0, 89, alpha=0.12, color="orange")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Seasonal Multiplier")
    ax.set_title("(A) Monsoon Seasonal Forcing (1.3x / 0.7x)")
    ax.legend(fontsize=8)
    ax.set_ylim(0.6, 1.4)

    # --- (B) Gravity network heatmap ---
    ax = axes[0, 1]
    im = ax.imshow(network, cmap="YlOrRd", aspect="auto",
                   origin="lower", interpolation="nearest")
    ax.set_xticks(range(nnodes))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title(f"(B) Gravity Network (b={PARAMS.gravity_b}, "
                 f"c={PARAMS.gravity_c}, cap=15%)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Migration weight")

    # --- (C) Weekly incidence by patch ---
    ax = axes[0, 2]
    nweeks = nticks // 7
    for i in range(nnodes):
        weekly = incidence[:nweeks * 7, i].reshape(nweeks, 7).sum(axis=1)
        weeks_x = np.arange(nweeks) / 52.0
        ax.plot(weeks_x, weekly, color=colors[i], alpha=0.7, lw=0.8,
                label=names[i])
    # Monsoon shading
    for yr in range(nyears):
        ax.axvspan(yr + 182 / 365, yr + 304 / 365, alpha=0.04, color="blue")
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly Infections")
    ax.set_title("(C) Weekly Incidence (blue = monsoon)")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (D) Susceptible fraction over time ---
    ax = axes[1, 0]
    sample30 = np.arange(0, nticks, 30)
    time30 = sample30 / 365.0
    for i in range(nnodes):
        S = model.nodes.S[sample30, i].astype(np.float64)
        E = model.nodes.E[sample30, i].astype(np.float64)
        I = model.nodes.I[sample30, i].astype(np.float64)
        R = model.nodes.R[sample30, i].astype(np.float64)
        N = np.maximum(S + E + I + R, 1.0)
        ax.plot(time30, S / N, color=colors[i], alpha=0.7, lw=0.8,
                label=names[i])
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.5,
               label=f"S* = 1/R0 = {1/R0:.3f}")
    ax.set_xlabel("Year")
    ax.set_ylabel("Susceptible Fraction (S/N)")
    ax.set_title("(D) Susceptible Fraction")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # --- (E) R_eff over time ---
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
    ax.set_title("(E) Effective Reproduction Number")
    ax.set_ylim(0, 4)
    ax.legend(fontsize=8)

    # --- (F) Network coupling vs distance ---
    ax = axes[1, 2]
    # Plot migration weight vs distance for all off-diagonal pairs
    dists_flat = []
    weights_flat = []
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                dists_flat.append(dist_matrix[i, j])
                weights_flat.append(network[i, j])
    ax.scatter(dists_flat, weights_flat, alpha=0.6, s=30, c="steelblue")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Migration Weight")
    ax.set_title(f"(F) Gravity Coupling vs Distance (c={PARAMS.gravity_c})")
    ax.set_yscale("log")

    # Overlay theoretical decay curve
    d_range = np.linspace(SPACING_KM, max(dists_flat), 100)
    # Theoretical: weight ∝ pop^b / d^c, all pops equal so just 1/d^c
    theoretical = (POP_PER_PATCH ** PARAMS.gravity_b) / (d_range ** PARAMS.gravity_c)
    # Scale to match data
    scale = np.median(weights_flat) / np.median(
        (POP_PER_PATCH ** PARAMS.gravity_b)
        / (np.array(dists_flat) ** PARAMS.gravity_c)
    )
    ax.plot(d_range, scale * theoretical, "r--", lw=1.5,
            label=f"d^{{-{PARAMS.gravity_c}}}")
    ax.legend(fontsize=8)

    plt.suptitle(
        f"10-Patch Polio SEIR — Gravity Migration + Monsoon Forcing\n"
        f"(R0={R0:.0f}, {POP_PER_PATCH/1000:.0f}K/patch, "
        f"{SPACING_KM:.0f}km spacing, "
        f"b={PARAMS.gravity_b}, c={PARAMS.gravity_c}, 15% cap)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()

    outdir = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "polio_gravity_seasonal_diagnostics.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nDiagnostic plot saved to {outpath}")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    model, scenario, season_365, network, dist_matrix = run_model()
    print_summary(model, scenario, network)
    plot_diagnostics(model, scenario, season_365, network, dist_matrix)
