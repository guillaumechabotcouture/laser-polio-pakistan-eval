The filesystem is read-only, so here is the complete, runnable script. It follows the exact patterns from your existing `scripts/polio_seir_basic_10patch.py` and `scripts/custom_components.py` (which has `PatchImportation`, not `Importation`).

Save this as `scripts/polio_spatial_seir_basic.py`:

```python
#!/usr/bin/env python3
"""
Spatial SEIR Model — Polio Transmission, 10 Pakistan Districts
==============================================================

10 patches of 100,000 each, 10-year simulation.

Epidemiological parameters (WPV1):
    R0 ≈ 6
    Latent period (E→I):     3 days  (gamma, shape=3, scale=1.0 → mean=3d)
    Infectious period (I→R): 28 days (normal, mean=28d, sd=3d)
    beta = R0 / D_infectious = 6 / 28 ≈ 0.2143 d⁻¹

Initial conditions (per patch):
    Population: 100,000
    R:  95,000  (95% immune — prior immunity/vaccination)
    I:       5  (active seed infections)
    E:       0
    S:   4,995  (remaining susceptible pool)

Vital dynamics (Pakistan averages):
    CBR: 29 per 1,000/year   CDR: 7 per 1,000/year

Spatial coupling:
    Gravity model (k=0.005, b=0.5, c=1.5), manually normalized so k equals
    the mean export fraction per node (~0.5%).

Periodic importation (all 10 patches):
    2 infections per patch every 30 days for the full 10 years.
    Required: initial R_eff = R0 × 0.05 = 0.30 < 1, so births must replenish
    susceptibles before sustained transmission is possible.  Without importation,
    stochastic extinction removes the virus before R_eff recovers to 1.

Usage:
    python scripts/polio_spatial_seir_basic.py
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
from laser.core.migration import gravity, row_normalizer, distance as haversine

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from custom_components import PatchImportation


# =============================================================================
# 1. Constants
# =============================================================================

N_PATCHES       = 10
POP_PER_PATCH   = 100_000
N_YEARS         = 10
NTICKS          = N_YEARS * 365

# Epidemiological
R0              = 6.0
LATENT_DAYS     = 3           # mean E→I duration, days
INFECTIOUS_DAYS = 28          # mean I→R duration, days
BETA            = R0 / INFECTIOUS_DAYS   # ≈ 0.2143 d⁻¹

# Initial conditions
IMMUNE_FRAC     = 0.95        # fraction placed in R at t=0
I_SEED          = 5           # infectious seed per patch

# Vital dynamics
CBR             = 29.0        # crude birth rate per 1,000/year
CDR             = 7.0         # crude death rate per 1,000/year

# Gravity network
GRAVITY_K       = 0.005       # mean export fraction per node
GRAVITY_B       = 0.5         # destination-population exponent
GRAVITY_C       = 1.5         # distance-decay exponent

# Importation
IMPORT_PERIOD   = 30          # ticks between events
IMPORT_COUNT    = 2           # infections per patch per event


# =============================================================================
# 2. Geographic scenario — 10 Pakistani districts
# =============================================================================

DISTRICTS = [
    # (name,          lat,     lon)
    ("Karachi",     24.86,   67.01),
    ("Lahore",      31.55,   74.34),
    ("Islamabad",   33.72,   73.06),
    ("Peshawar",    34.01,   71.57),
    ("Quetta",      30.19,   67.01),
    ("Multan",      30.20,   71.47),
    ("Faisalabad",  31.42,   73.08),
    ("Rawalpindi",  33.60,   73.06),
    ("Hyderabad",   25.38,   68.37),
    ("Sialkot",     32.49,   74.54),
]
assert len(DISTRICTS) == N_PATCHES

LATS = np.array([d[1] for d in DISTRICTS])
LONS = np.array([d[2] for d in DISTRICTS])


def build_scenario():
    """GeoDataFrame with one row per patch and initial SEIR compartments."""
    R_init = int(IMMUNE_FRAC * POP_PER_PATCH)          # 95,000
    S_init = POP_PER_PATCH - R_init - I_SEED            # 4,995
    assert S_init >= 0, f"Negative initial S: {S_init}"

    rows = [
        {
            "nodeid":     i,
            "name":       name,
            "population": POP_PER_PATCH,
            "lat":        lat,
            "lon":        lon,
            "geometry":   Point(lon, lat),
            "S":          S_init,
            "E":          0,
            "I":          I_SEED,
            "R":          R_init,
        }
        for i, (name, lat, lon) in enumerate(DISTRICTS)
    ]
    scenario = gpd.GeoDataFrame(rows, crs="EPSG:4326")

    # Required invariant: S + E + I + R == population for every patch
    assert (scenario.S + scenario.E + scenario.I + scenario.R
            == scenario.population).all(), "Compartment sum mismatch"
    assert (scenario.I > 0).all(), "No seed infections"

    return scenario, S_init, R_init


# =============================================================================
# 3. Age pyramid (stable exponential, Pakistan ~65-year life expectancy)
# =============================================================================

def build_age_pyramid():
    """Approximate Pakistan age pyramid (integer weights per year of age)."""
    ages    = np.arange(100)
    weights = np.array(1000 * np.exp(-ages / 65.0), dtype=np.int64)
    return AliasedDistribution(np.maximum(weights, 1))


def initialize_ages(model, pyramid):
    """Assign realistic initial ages before model.run().

    Without this every agent has dob=0 (appears as a newborn), which
    distorts mortality and any age-targeted vaccination added later.
    """
    count      = model.people.count
    ages_years = np.minimum(pyramid.sample(count=count, dtype=np.int32), 99)
    ages_days  = ages_years * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


# =============================================================================
# 4. Gravity network (manual normalization)
# =============================================================================

def build_gravity_network(scenario):
    """Pairwise Haversine distances → gravity migration matrix.

    Returns (nnodes × nnodes) matrix where [i, j] is the fraction of
    patch i's infectious pressure that reaches patch j.
    """
    nnodes      = len(scenario)
    dist_matrix = np.zeros((nnodes, nnodes))
    for i in range(nnodes):
        for j in range(nnodes):
            if i != j:
                dist_matrix[i, j] = haversine(
                    LATS[i], LONS[i], LATS[j], LONS[j]
                )

    pops    = np.array(scenario.population, dtype=np.float64)
    network = gravity(pops, dist_matrix, 1, 0, GRAVITY_B, GRAVITY_C)

    # Normalize: GRAVITY_K == mean export fraction per node
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * GRAVITY_K
    network = row_normalizer(network, 0.2)   # cap at 20 % per row

    # Validation — catches silently broken networks
    assert network.sum() > 0, \
        "Gravity network is all zeros — check GRAVITY_K or coordinates"
    assert network.sum(axis=1).max() < 0.3, \
        f"Network row sum too high: {network.sum(axis=1).max():.3f}"

    return network


# =============================================================================
# 5. Build and run
# =============================================================================

def run_model():
    scenario, S_init, R_init = build_scenario()
    nnodes = len(scenario)

    print("=" * 64)
    print("Spatial Polio SEIR — 10 Districts, 10-year simulation")
    print("=" * 64)
    print(f"  Patches : {nnodes} × {POP_PER_PATCH:,} = {nnodes*POP_PER_PATCH:,} total")
    print(f"  R0      : {R0}  →  beta = {BETA:.4f} d⁻¹")
    print(f"  Latent  : ~{LATENT_DAYS}d (gamma)  |  Infectious: ~{INFECTIOUS_DAYS}d (normal)")
    print(f"  t=0/patch: S={S_init:,}  E=0  I={I_SEED}  R={R_init:,}")
    print(f"  Immune  : {IMMUNE_FRAC:.0%}  |  Initial R_eff = {R0*(1-IMMUNE_FRAC):.2f}")
    print(f"  CBR={CBR}/1000yr  CDR={CDR}/1000yr")

    # --- Duration distributions ---
    # gamma(shape=3, scale=1.0): mean = shape × scale = 3 days
    expdurdist = dists.gamma(shape=3, scale=float(LATENT_DAYS) / 3)
    # normal(mean=28, sd=3): mean = 28 days
    infdurdist = dists.normal(loc=float(INFECTIOUS_DAYS), scale=3.0)

    # --- Vital dynamics arrays (nticks × nnodes, per-1000/year) ---
    birthrate_array = np.full((NTICKS, nnodes), CBR, dtype=np.float32)
    deathrate_array = np.full((NTICKS, nnodes), CDR, dtype=np.float32)
    # Unit assertion: BirthsByCBR expects per-1000/year (common silent error)
    assert (birthrate_array >= 1).all() and (birthrate_array <= 60).all(), \
        "Birthrates must be per-1000/year (typical range 10–50)"

    pyramid = build_age_pyramid()

    # --- Parameters (no gravity keys — network set manually below) ---
    params = PropertySet({
        "prng_seed":              42,
        "nticks":                 NTICKS,
        "beta":                   BETA,
        "capacity_safety_factor": 2.5,   # pre-allocate for 10-year birth growth
    })

    # --- Instantiate model ---
    model = Model(scenario, params, birthrates=birthrate_array)

    # --- Gravity network (manual, for custom normalization) ---
    model.network = build_gravity_network(scenario)
    print(f"\n  Network export fractions: "
          f"min={model.network.sum(axis=1).min():.4f}  "
          f"max={model.network.sum(axis=1).max():.4f}")

    # --- No seasonal forcing (flat profile) ---
    # Polio has weak seasonality; add monsoon forcing later if needed.
    seasonality = ValuesMap.from_scalar(1.0, NTICKS, nnodes)

    # --- Component pipeline ---
    # Order preserves the S + E + I + R = N invariant each tick.
    model.components = [
        SEIR.Susceptible(model),                              # propagate S[t+1]
        SEIR.Exposed(model, expdurdist, infdurdist),         # E→I
        SEIR.Infectious(model, infdurdist),                   # I→R
        SEIR.Recovered(model),                                # propagate R[t+1]
        BirthsByCBR(model,
                    birthrates=birthrate_array,
                    pyramid=pyramid),
        MortalityByCDR(model,
                       mortalityrates=deathrate_array),
        PatchImportation(model, infdurdist,
                         endemic_patches=list(range(nnodes)),
                         period=IMPORT_PERIOD,
                         count=IMPORT_COUNT,
                         end_tick=NTICKS),
        SEIR.Transmission(model, expdurdist,                 # S→E via FOI
                          seasonality=seasonality),
    ]

    # Assign realistic initial ages BEFORE run() (must follow component setup
    # so that BirthsByCBR has already created model.people.dob)
    initialize_ages(model, pyramid)

    print(f"\nRunning {N_YEARS}-year simulation ({NTICKS} ticks)…")
    model.run("Polio SEIR 10-patch")
    print("Done.\n")

    return model, scenario


# =============================================================================
# 6. Summary statistics
# =============================================================================

def print_summary(model, scenario):
    nticks = NTICKS
    nnodes = len(scenario)
    names  = list(scenario.name)
    nyears = nticks // 365
    inc    = model.nodes.newly_infected    # (nticks, nnodes)

    print("=" * 64)
    print("SIMULATION SUMMARY")
    print("=" * 64)
    print(f"\n{'District':<13} {'Mean Inf/yr':<15} {'Total Inf':<12} {'S/N (final)'}")
    print("-" * 54)
    grand_total = 0
    for i in range(nnodes):
        annual      = np.array([inc[y*365:(y+1)*365, i].sum() for y in range(nyears)])
        total       = int(annual.sum())
        grand_total += total
        t_fin       = nticks - 1
        S_f = float(model.nodes.S[t_fin, i])
        N_f = float(model.nodes.S[t_fin, i] + model.nodes.E[t_fin, i]
                    + model.nodes.I[t_fin, i] + model.nodes.R[t_fin, i])
        print(f"{names[i]:<13} {annual.mean():<15.0f} {total:<12,} {S_f/max(N_f,1):.3f}")
    print("-" * 54)
    print(f"{'TOTAL':<13} {'':<15} {grand_total:,}")

    # Population trajectory
    print(f"\nPopulation (every 2 years):")
    for y in range(0, nyears + 1, 2):
        t   = min(y * 365, nticks - 1)
        pop = int(model.nodes.S[t].sum() + model.nodes.E[t].sum()
                  + model.nodes.I[t].sum() + model.nodes.R[t].sum())
        print(f"  Year {y:>2d}: {pop:>10,}")
    print(f"  Agent capacity: {model.people.capacity:,}  "
          f"Active: {model.people.count:,}")

    # Compartment integrity at key ticks
    print(f"\nCompartment integrity (S+E+I+R = N):")
    for tick in [0, nticks // 2, nticks - 1]:
        S = model.nodes.S[tick].sum()
        E = model.nodes.E[tick].sum()
        I = model.nodes.I[tick].sum()
        R = model.nodes.R[tick].sum()
        print(f"  Tick {tick:>5}: S={S:>8,}  E={E:>5,}  "
              f"I={I:>5,}  R={R:>9,}  N={S+E+I+R:>10,}")

    # Compartment non-negativity
    ok = True
    for label, arr in [("S", model.nodes.S), ("E", model.nodes.E),
                        ("I", model.nodes.I), ("R", model.nodes.R)]:
        if np.any(arr[:nticks] < 0):
            print(f"  WARNING: negative {label} values detected")
            ok = False
    if ok:
        print("  All compartments non-negative: OK")

    # National R_eff during second half of simulation
    ticks30 = np.arange(nticks // 2, nticks, 30)
    S_t = model.nodes.S[ticks30].sum(axis=1).astype(float)
    N_t = (model.nodes.S[ticks30].sum(axis=1)
           + model.nodes.E[ticks30].sum(axis=1)
           + model.nodes.I[ticks30].sum(axis=1)
           + model.nodes.R[ticks30].sum(axis=1)).astype(float)
    reff = R0 * S_t / np.maximum(N_t, 1)
    print(f"\n  National R_eff (years 5–10): "
          f"mean={reff.mean():.2f}  "
          f"range=[{reff.min():.2f}, {reff.max():.2f}]")


# =============================================================================
# 7. Diagnostic plots
# =============================================================================

def plot_results(model, scenario):
    nticks  = NTICKS
    nnodes  = len(scenario)
    names   = list(scenario.name)
    nyears  = nticks // 365
    nweeks  = nticks // 7
    inc     = model.nodes.newly_infected

    weekly_national = inc[:nweeks*7].reshape(nweeks, 7, nnodes).sum(axis=(1, 2))
    weekly_by_patch = inc[:nweeks*7].reshape(nweeks, 7, nnodes).sum(axis=1)
    weeks_x         = np.arange(nweeks) / 52.0    # in years

    ticks30 = np.arange(0, nticks, 30)
    time30  = ticks30 / 365.0

    cmap   = plt.cm.tab10
    colors = [cmap(i) for i in range(nnodes)]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(
        f"Spatial Polio SEIR — 10 Pakistan Districts\n"
        f"R₀={R0}, latent={LATENT_DAYS}d, infectious={INFECTIOUS_DAYS}d, "
        f"{N_YEARS}-yr  |  {N_PATCHES} patches × {POP_PER_PATCH//1000}K",
        fontsize=12, fontweight="bold",
    )

    # (A) National weekly incidence
    ax = axes[0, 0]
    ax.plot(weeks_x, weekly_national, color="steelblue", lw=1)
    ax.set_xlabel("Year")
    ax.set_ylabel("New infections / week")
    ax.set_title("(A) National Weekly Incidence")
    ax.set_xlim(0, nyears)

    # (B) Susceptible fraction per district
    ax = axes[0, 1]
    for i in range(nnodes):
        S = model.nodes.S[ticks30, i].astype(float)
        N = (model.nodes.S[ticks30, i] + model.nodes.E[ticks30, i]
             + model.nodes.I[ticks30, i] + model.nodes.R[ticks30, i]).astype(float)
        ax.plot(time30, S / np.maximum(N, 1), color=colors[i],
                lw=0.9, alpha=0.8, label=names[i])
    ax.axhline(1 / R0, color="black", ls=":", alpha=0.5,
               label=f"S* = 1/R₀ = {1/R0:.3f}")
    ax.set_xlabel("Year")
    ax.set_ylabel("S / N")
    ax.set_title("(B) Susceptible Fraction by District")
    ax.legend(fontsize=6, ncol=2, loc="upper left")
    ax.set_xlim(0, nyears)

    # (C) R_eff per district
    ax = axes[1, 0]
    for i in range(nnodes):
        S = model.nodes.S[ticks30, i].astype(float)
        N = (model.nodes.S[ticks30, i] + model.nodes.E[ticks30, i]
             + model.nodes.I[ticks30, i] + model.nodes.R[ticks30, i]).astype(float)
        ax.plot(time30, R0 * S / np.maximum(N, 1),
                color=colors[i], lw=0.9, alpha=0.8, label=names[i])
    ax.axhline(1.0, color="black", lw=2, alpha=0.35, label="R_eff = 1")
    ax.set_xlabel("Year")
    ax.set_ylabel("R_eff = R₀ × S/N")
    ax.set_title("(C) Effective Reproduction Number")
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=6, ncol=2, loc="upper left")
    ax.set_xlim(0, nyears)

    # (D) Spatial heatmap — weekly incidence per district
    ax = axes[1, 1]
    im = ax.imshow(weekly_by_patch.T, aspect="auto", cmap="YlOrRd",
                   origin="lower", extent=[0, nyears, -0.5, nnodes - 0.5])
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Year")
    ax.set_title("(D) Weekly Incidence Heatmap by District")
    plt.colorbar(im, ax=ax, label="Cases / week")

    plt.tight_layout()
    outdir  = Path(__file__).parent.parent / "eval" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "polio_spatial_seir_basic.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Plot saved → {outpath}")
    plt.close(fig)


# =============================================================================
# 8. Entry point
# =============================================================================

if __name__ == "__main__":
    model, scenario = run_model()
    print_summary(model, scenario)
    plot_results(model, scenario)
```

---

**Key design decisions and why:**

| Decision | Rationale |
|---|---|
| `BETA = R0 / INFECTIOUS_DAYS = 6/28 ≈ 0.2143` | LASER's FOI is `beta × I/N`; this recovers R0 exactly at full susceptibility |
| `gamma(shape=3, scale=1.0)` for latent period | Mean = 3d, sd ≈ 1.73d; gamma gives realistic right-skewed residence time |
| `normal(mean=28, scale=3)` for infectious period | Matches polio shedding duration; `infdurmin=1` in `InfectiousIR` floors at 1d |
| `PatchImportation` (all 10 patches, count=2, period=30) | Initial R_eff = 6 × 0.05 = **0.30 < 1**; without importation, stochastic extinction removes the virus before births replenish the susceptible pool |
| Manual gravity network (not params auto-setup) | Lets you normalize so `GRAVITY_K` has a clear physical meaning (mean export fraction) |
| `initialize_ages()` before `run()` | Prevents all 1M agents appearing as newborns, which would cause a spurious mortality pulse |
| `capacity_safety_factor=2.5` | 10 years × (29−7)/1000 ≈ 22% net growth per decade; 2.5× gives safe headroom |
| `mortalityrates=` (not `deathrates=`) | `MortalityByCDR` uses this keyword; the wrong name silently ignores deaths |
