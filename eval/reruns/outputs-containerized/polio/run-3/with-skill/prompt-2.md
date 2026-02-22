The workspace filesystem is read-only, so I'm presenting the complete code here. Save it as `polio_gravity_seasonal_model.py`.

---

Three design decisions worth noting before the code:

**Gravity formula:** `M_{i,j} = k · pᵢ · pⱼ^0.5 / dᵢⱼ^1.5` (source exponent `a=1`, destination `b=0.5`, decay `c=1.5`). The network is built with `k=1`, re-scaled so the mean row-sum equals a target coupling fraction of 0.01, then capped at 15 % by `row_normalizer`.

**Seasonal profile:** `season(t) = 1 + 0.3·cos(2π(t−242)/365)`. Day 242 is late August (centre of Jul–Oct); the cosine trough falls at day ~60 (late February, within Dec–Mar). The cosine's integral over a full period is zero, so the mean is exactly 1.0 without approximation—no bias to R_eff.

**Distance matrix:** computed directly as `|i−j|×50 km` (exact for a linear arrangement; no haversine approximation needed).

```python
#!/usr/bin/env python3
"""
Pakistan Polio Corridor — 6-Patch Spatial SEIR Model

Extends the basic LASER polio model with three features:
    1. Gravity-model migration network (b=0.5, c=1.5, source a=1)
    2. Monsoon seasonal forcing (Jul–Oct peak 1.3×, Dec–Mar trough 0.7×)
    3. Row-normalised network capped at 15 % outward FOI export per patch

Patch geometry:
    6 districts in a linear corridor, 50 km apart.
    Distance matrix computed directly as |i − j| × 50 km (exact for a line).

Poliovirus epidemiology:
    R0 ≈ 8  (endemic Pakistan, low OPV coverage)
    Latent period  ≈ 7 d   — gamma(shape=3, scale=7/3)
    Infectious period ≈ 14 d — gamma(shape=4, scale=3.5)
    beta = R0 / inf_mean = 8 / 14 ≈ 0.571

Vital dynamics:
    CBR = 30 / 1000 / yr   CDR = 8 / 1000 / yr   (Pakistan national rates)

Duration:  10 years  (3 650 ticks)

Usage:
    python polio_gravity_seasonal_model.py
"""

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


# ============================================================================
# PatchImportation — periodic seeding to maintain endemic transmission
# ============================================================================

class PatchImportation:
    """Seeds infections into specified patches at regular intervals.

    Every `period` ticks, selects up to `count` susceptible agents in each
    target patch and transitions them to INFECTIOUS with sampled infectious
    durations.  Stops seeding after `end_tick`.

    Parameters
    ----------
    model           : LASER Model instance
    infdurdist      : Distribution for infectious duration sampling
    endemic_patches : Patch indices to seed  (default: [0])
    period          : Ticks between importation events  (default: 30)
    count           : Infections per patch per event  (default: 2)
    end_tick        : Stop after this tick  (default: full simulation)
    """

    def __init__(self, model, infdurdist, endemic_patches=None,
                 period=30, count=2, end_tick=None):
        self.model           = model
        self.infdurdist      = infdurdist
        self.endemic_patches = np.asarray(
            endemic_patches if endemic_patches is not None else [0],
            dtype=np.int32,
        )
        self.period   = period
        self.count    = count
        self.end_tick = end_tick if end_tick is not None else model.params.nticks

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0 or tick >= self.end_tick:
            return

        people = self.model.people
        nodes  = self.model.nodes
        n_live = people.count

        for patch in self.endemic_patches:
            eligible = np.nonzero(
                (people.state[:n_live]  == SEIR.State.SUSCEPTIBLE.value) &
                (people.nodeid[:n_live] == patch)
            )[0]
            if len(eligible) == 0:
                continue

            n_infect = min(self.count, len(eligible))
            chosen   = np.random.choice(eligible, size=n_infect, replace=False)

            people.state[chosen] = SEIR.State.INFECTIOUS.value

            samples = dists.sample_floats(
                self.infdurdist, np.zeros(n_infect, np.float32)
            )
            people.itimer[chosen] = np.maximum(
                np.round(samples), 1
            ).astype(people.itimer.dtype)

            nodes.S[tick + 1, patch] -= n_infect
            nodes.S[tick + 1, patch]  = max(nodes.S[tick + 1, patch], 0)
            nodes.I[tick + 1, patch] += n_infect


# ============================================================================
# Patch Configuration — 6 districts in a 50-km-spaced corridor
# ============================================================================

NNODES         = 6
INTER_PATCH_KM = 50.0      # inter-patch spacing (km)

patches = [
    {"nodeid": 0, "name": "Dera_Ismail",  "population": 120_000},
    {"nodeid": 1, "name": "Jampur",        "population": 250_000},
    {"nodeid": 2, "name": "Rajanpur",      "population": 180_000},
    {"nodeid": 3, "name": "Fazilpur",      "population": 310_000},
    {"nodeid": 4, "name": "Muzaffargarh",  "population": 420_000},
    {"nodeid": 5, "name": "Alipur",        "population":  90_000},
]

# Lon spacing for exactly 50 km at 30 °N  (1 ° lon ≈ 96.1 km at 30 °N)
BASE_LAT, BASE_LON = 30.0, 67.0
LON_PER_50KM = 50.0 / (111.0 * np.cos(np.radians(BASE_LAT)))   # ≈ 0.519 °

for k, p in enumerate(patches):
    p["lat"] = BASE_LAT
    p["lon"] = BASE_LON + k * LON_PER_50KM


# ============================================================================
# Feature 1 — Pairwise Distance Matrix  (exact linear geometry)
# ============================================================================
#
# d[i, j] = |i − j| × 50 km, which is exact for a linear arrangement.
# Diagonal set to ∞ so that gravity() yields M_{i,i} = 0 (no self-coupling).

dist_matrix = np.array(
    [[abs(i - j) * INTER_PATCH_KM for j in range(NNODES)]
     for i in range(NNODES)],
    dtype=np.float64,
)
np.fill_diagonal(dist_matrix, np.inf)


# ============================================================================
# Scenario GeoDataFrame
# ============================================================================

geometry = [Point(p["lon"], p["lat"]) for p in patches]
scenario = gpd.GeoDataFrame(patches, geometry=geometry, crs="EPSG:4326")

# Initial conditions for endemic low-coverage setting:
#   S ≈ 15 %,  E = 0,  I ≈ 0.2 %,  R ≈ 84.8 %
scenario["E"] = np.zeros(NNODES, dtype=np.uint32)
scenario["I"] = np.maximum(
    np.round(0.002 * scenario["population"]).astype(np.uint32), 1
)
scenario["R"] = np.round(0.848 * scenario["population"]).astype(np.uint32)
scenario["S"] = (
    scenario["population"] - scenario["E"] - scenario["I"] - scenario["R"]
).astype(np.uint32)

assert (scenario.S + scenario.E + scenario.I + scenario.R
        == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch needs initial infections"

names = list(scenario.name)

print("Patch configuration (6 districts, 50 km apart):")
for _, row in scenario.iterrows():
    print(f"  {row['name']:<18s}: pop={row['population']:>7,}  "
          f"S={row['S']:>6,}  I={int(row['I']):>4}  R={row['R']:>7,}")


# ============================================================================
# Simulation Parameters
# ============================================================================

NTICKS   = 10 * 365        # 10-year simulation
CBR      = 30.0            # crude birth rate  (per 1 000 / yr)
CDR      =  8.0            # crude death rate  (per 1 000 / yr)
R0       =  8.0            # poliovirus R0 in low-immunity setting
INF_MEAN = 14.0            # mean infectious period (days)
BETA     = R0 / INF_MEAN   # ≈ 0.571

params = PropertySet({
    "prng_seed":              42,
    "nticks":                 NTICKS,
    "beta":                   BETA,
    "capacity_safety_factor": 3.0,
})


# ============================================================================
# Birth / Death Rates  (per 1 000 / yr — LASER divides by 1 000 internally)
# ============================================================================

birthrate_array = np.full((NTICKS, NNODES), CBR, dtype=np.float32)
mortality_array = np.full((NTICKS, NNODES), CDR, dtype=np.float32)

# Unit guard: must be per-1000/yr, NOT daily per-capita
assert np.all(birthrate_array >= 1) and np.all(birthrate_array <= 60), \
    f"Birthrates out of range: {birthrate_array.min():.4f}–{birthrate_array.max():.4f}"
assert np.all(mortality_array >= 1) and np.all(mortality_array <= 60), \
    f"Death rates out of range: {mortality_array.min():.4f}–{mortality_array.max():.4f}"

# Age pyramid: exponential with ~65-yr life expectancy
ages    = np.arange(100)
pyramid = AliasedDistribution(
    np.maximum(np.round(1000 * np.exp(-ages / 65.0)).astype(np.int64), 1)
)


# ============================================================================
# Feature 2 — Monsoon Seasonal Forcing
# ============================================================================
#
# Profile:  season(t) = 1 + A · cos( 2π(t − peak_day) / 365 )
#
# Window alignment:
#   Jul 1  = day 181,  Oct 31 = day 303  →  centre ≈ day 242  (late August)
#   Dec 1  = day 334,  Mar 31 = day  89  →  trough ≈ day  60  (late February)
#
# With peak_day = 242, A = 0.3:
#   t = 242  →  1 + 0.3·cos(0)   = 1.30   (peak in Jul–Oct)    ✓
#   t =  60  →  1 + 0.3·cos(π)  ≈ 0.70   (trough in Dec–Mar)  ✓
#   mean     = 1.0 exactly  (cosine integrates to 0 over a full period)
#
# The profile is normalised to mean == 1.0 before use.  An un-normalised
# profile biases R_eff systematically (LASER Layer 1.1 warning).

PEAK_DAY  = 242      # late August — centre of the Jul–Oct monsoon window
AMPLITUDE = 0.3      # yields peak = 1.3× and trough = 0.7×

days       = np.arange(365)
season_365 = 1.0 + AMPLITUDE * np.cos(2 * np.pi * (days - PEAK_DAY) / 365)
season_365 = season_365 / season_365.mean()   # normalise (should already be 1.0)

# ---- Assertions ----
assert abs(season_365.mean() - 1.0) < 0.01, \
    f"Seasonal profile mean={season_365.mean():.4f} — must be ~1.0"
assert abs(season_365.max() - 1.30) < 0.02, \
    f"Seasonal peak={season_365.max():.3f} — expected ~1.30 (Jul–Oct)"
assert abs(season_365.min() - 0.70) < 0.02, \
    f"Seasonal trough={season_365.min():.3f} — expected ~0.70 (Dec–Mar)"

peak_idx   = int(np.argmax(season_365))
trough_idx = int(np.argmin(season_365))
assert 181 <= peak_idx <= 303, \
    f"Seasonal peak at day {peak_idx} — expected within Jul–Oct (181–303)"
assert trough_idx <= 90 or trough_idx >= 334, \
    f"Seasonal trough at day {trough_idx} — expected within Dec–Mar"

print(f"\nMonsoon seasonal forcing:")
print(f"  Peak day:          day {peak_idx}  (≈ late August)")
print(f"  Peak multiplier:   {season_365.max():.3f}×  (Jul–Oct monsoon season)")
print(f"  Trough day:        day {trough_idx}  (≈ late February)")
print(f"  Trough multiplier: {season_365.min():.3f}×  (Dec–Mar dry season)")
print(f"  Annual mean:       {season_365.mean():.4f}  (normalised)")

# Tile across simulation and wrap as ValuesMap for SEIR.Transmission
season_tiled = np.tile(season_365, NTICKS // 365 + 1)[:NTICKS]
seasonality  = ValuesMap.from_timeseries(season_tiled, NNODES)


# ============================================================================
# Build LASER Model
# ============================================================================

model = Model(scenario, params, birthrates=birthrate_array)


# ============================================================================
# Feature 1 — Gravity Migration Network
# ============================================================================
#
# Formula:  M_{i,j} = k · pᵢ^a · pⱼ^b / dᵢⱼ^c
#   a = 1.0  (source population exponent — standard gravity convention)
#   b = 0.5  (destination population exponent — specified)
#   c = 1.5  (distance decay exponent — specified)
#
# Construction:
#   Step 1. Build raw network with k=1  (scale is arbitrary at this stage)
#   Step 2. Re-scale so mean row-sum = GRAVITY_K = 0.01
#           (on average each patch exports 1 % of its FOI)
#   Step 3. row_normalizer(0.15) enforces Feature 3

GRAVITY_A = 1.0     # source population exponent
GRAVITY_B = 0.5     # destination population exponent  (per spec)
GRAVITY_C = 1.5     # distance decay exponent          (per spec)
GRAVITY_K = 0.01    # target mean outward coupling fraction

pops        = np.array(scenario.population, dtype=np.float64)
raw_network = gravity(pops, dist_matrix, 1.0, GRAVITY_A, GRAVITY_B, GRAVITY_C)

avg_row_sum = raw_network.sum(axis=1).mean()
if avg_row_sum == 0:
    raise ValueError(
        "Raw gravity network is all zeros — check distance matrix and populations"
    )
network = raw_network * (GRAVITY_K / avg_row_sum)


# ============================================================================
# Feature 3 — Row-Normalise: cap FOI export at 15 % per patch
# ============================================================================
#
# row_normalizer(network, f) scales each row so that row_sum <= f.
# This prevents small patches adjacent to large population centres from
# losing an unrealistic share of their force of infection to neighbours.
# The 15 % cap is applied after the k-scaling step above.

MAX_EXPORT = 0.15   # maximum outward FOI fraction per patch

network       = row_normalizer(network, MAX_EXPORT)
model.network = network

# ---- Assertions ----
assert model.network.sum() > 0, \
    "Network is all zeros after normalisation — check gravity parameters"

row_sums = model.network.sum(axis=1)
assert row_sums.max() <= MAX_EXPORT + 1e-9, \
    (f"Max row sum {row_sums.max():.6f} exceeds {MAX_EXPORT:.2f} cap "
     f"after row_normalizer — unexpected numerical issue")
assert row_sums.max() < 0.3, \
    "Network safety check: max row sum must be < 0.3"

print(f"\nGravity network  (a={GRAVITY_A}, b={GRAVITY_B}, c={GRAVITY_C}, "
      f"max_export={MAX_EXPORT:.0%}):")
print(f"  {'Patch':<18}  {'Row sum':>10}  {'Max element':>13}")
for i, name in enumerate(names):
    print(f"  {name:<18}  {row_sums[i]:>10.6f}  "
          f"{model.network[i].max():>13.6f}")
print(f"  Max row sum: {row_sums.max():.6f}  (cap enforced at {MAX_EXPORT:.2f}) ✓")

print(f"\nGravity network matrix  (row = source, col = destination):")
header = "  ".join(f"{n:>14}" for n in names)
print(f"  {'':>18} {header}")
for i, src in enumerate(names):
    row_str = "  ".join(f"{model.network[i, j]:>14.6f}" for j in range(NNODES))
    print(f"  {src:<18} {row_str}")


# ============================================================================
# Duration Distributions — Poliovirus type 1
# ============================================================================

# Latent period:    gamma(shape=3, scale=7/3) → mean = 7 d, CV ≈ 0.58
expdurdist = dists.gamma(shape=3, scale=7.0 / 3.0)

# Infectious period: gamma(shape=4, scale=3.5) → mean = 14 d, CV = 0.5
infdurdist = dists.gamma(shape=4, scale=3.5)


# ============================================================================
# Assemble Component Pipeline
# ============================================================================
#
# Component ordering:
#   1–4. Susceptible / Exposed / Infectious / Recovered propagate compartment
#        counts (S[t+1]←S[t], …) and advance state machines (E→I, I→R).
#     5. PatchImportation seeds new infections before Transmission runs.
#     6. Transmission computes FOI — internally uses model.network (gravity)
#        and seasonality (monsoon), then draws S→E transitions.
#   7–8. Vital dynamics run last so the cohort is stable during disease steps.

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    PatchImportation(
        model, infdurdist,
        endemic_patches=[0, 3],   # seed two corridor patches
        period=30, count=2,
    ),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortality_array),
]


# ============================================================================
# Run Simulation
# ============================================================================

total_pop = int(scenario.population.sum())
print(f"\nRunning 10-year Pakistan Polio Corridor SEIR simulation ...")
print(f"  Patches:          {NNODES} districts × {INTER_PATCH_KM:.0f} km spacing")
print(f"  Total population: {total_pop:,}")
print(f"  R0={R0}, beta={BETA:.4f}, latent≈7 d, infectious≈{INF_MEAN:.0f} d")
print(f"  CBR={CBR}/1 000/yr, CDR={CDR}/1 000/yr")
print(f"  Gravity: a={GRAVITY_A}, b={GRAVITY_B}, c={GRAVITY_C}, "
      f"max_export={MAX_EXPORT:.0%}")
print(f"  Seasonality: monsoon peak {season_365.max():.2f}× (Jul–Oct), "
      f"trough {season_365.min():.2f}× (Dec–Mar)")
print(f"  Agent capacity:   {model.people.capacity:,}")

model.run("Pakistan Polio Corridor SEIR")
print("Simulation complete.")


# ============================================================================
# Post-Run Verification
# ============================================================================

print("\n--- Post-Run Verification ---")

# 1. Population trajectory (should grow: CBR=30 > CDR=8)
pop_t0 = int(
    model.nodes.S[0].sum()         + model.nodes.E[0].sum()
    + model.nodes.I[0].sum()       + model.nodes.R[0].sum()
)
pop_tf = int(
    model.nodes.S[NTICKS - 1].sum() + model.nodes.E[NTICKS - 1].sum()
    + model.nodes.I[NTICKS - 1].sum() + model.nodes.R[NTICKS - 1].sum()
)
growth_pct = (pop_tf / pop_t0 - 1) * 100
print(f"Population: {pop_t0:,} → {pop_tf:,}  ({growth_pct:+.1f}%)")
if pop_tf > pop_t0:
    print("  Population growing: OK  (CBR > CDR)")
else:
    print("  WARNING: population not growing — check birth/death rate units")

# 2. Compartment non-negativity
any_neg = False
for arr, cname in [(model.nodes.S, "S"), (model.nodes.E, "E"),
                   (model.nodes.I, "I"), (model.nodes.R, "R")]:
    if np.any(arr[:NTICKS] < 0):
        print(f"  WARNING: negative {cname} values detected!")
        any_neg = True
if not any_neg:
    print("  Compartment non-negativity: OK")

# 3. Epidemic dynamics
total_inf = int(model.nodes.newly_infected[:NTICKS].sum())
print(f"  Total infections (10 yr): {total_inf:,}")
if total_inf > 0:
    print("  Transmission active: OK")
else:
    print("  WARNING: zero infections — check beta or initial conditions")

# 4. Spatial coupling
patch_totals = model.nodes.newly_infected[:NTICKS].sum(axis=0)
n_active     = int(np.sum(patch_totals > 0))
print(f"  Patches with infections: {n_active}/{NNODES}")
if n_active == NNODES:
    print("  Spatial coupling active across all patches: OK")
else:
    print("  WARNING: some patches have zero infections — check network")

# 5. Network row-sums still within bounds after run
rs = model.network.sum(axis=1)
assert rs.max() <= MAX_EXPORT + 1e-9, \
    f"Network row sum {rs.max():.6f} exceeds {MAX_EXPORT:.2f} cap"
print(f"  Network max row sum: {rs.max():.6f}  (≤ {MAX_EXPORT:.2f}) ✓")

# 6. Compartment snapshots
print(f"\nCompartment counts:")
print(f"  {'Tick':<8} {'S':>10} {'E':>8} {'I':>8} {'R':>10} {'Total':>10}")
print(f"  {'-' * 58}")
for tick in [0, NTICKS // 4, NTICKS // 2, 3 * NTICKS // 4, NTICKS - 1]:
    S = int(model.nodes.S[tick].sum())
    E = int(model.nodes.E[tick].sum())
    I = int(model.nodes.I[tick].sum())
    R = int(model.nodes.R[tick].sum())
    print(f"  {tick:<8} {S:>10,} {E:>8,} {I:>8,} {R:>10,} {S+E+I+R:>10,}")

# 7. Per-patch annual infection rates
nyears = NTICKS // 365
print(f"\nAnnual infections by patch:")
print(f"  {'Patch':<18} {'Pop':>8} {'Avg/yr':>10} {'Per million':>12}")
print(f"  {'-' * 52}")
for i, name in enumerate(names):
    annual = np.array([
        int(model.nodes.newly_infected[y * 365:(y + 1) * 365, i].sum())
        for y in range(nyears)
    ])
    avg  = annual.mean()
    rate = avg / pops[i] * 1e6
    print(f"  {name:<18} {int(pops[i]):>8,} {avg:>10,.1f} {rate:>12,.0f}")


# ============================================================================
# Diagnostic Plots
# ============================================================================

COLORS = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c"]

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MONTH_STARTS = [0, 31, 59, 90, 120, 151, 181, 212, 242, 273, 303, 334]

# (A) Monsoon seasonal forcing profile
ax = axes[0, 0]
ax.plot(days, season_365, color="steelblue", lw=2.5)
ax.axhline(1.0, color="gray", ls="--", alpha=0.5, lw=1)
ax.axhline(1.3, color="red",  ls=":",  alpha=0.7, lw=1.5,
           label="Peak 1.3× (Jul–Oct)")
ax.axhline(0.7, color="navy", ls=":",  alpha=0.7, lw=1.5,
           label="Trough 0.7× (Dec–Mar)")
ax.axvspan(181, 303, alpha=0.12, color="red",  label="Jul–Oct monsoon")
ax.axvspan(0,    90, alpha=0.12, color="navy", label="Dec–Mar trough")
ax.axvspan(334, 365, alpha=0.12, color="navy")
ax.set_xticks(MONTH_STARTS)
ax.set_xticklabels(MONTH_LABELS, fontsize=8)
ax.set_xlabel("Month of year")
ax.set_ylabel("Transmission multiplier")
ax.set_title("(A) Monsoon Seasonal Forcing")
ax.legend(fontsize=8)
ax.set_ylim(0.5, 1.55)

# (B) Daily new infections by patch (7-day rolling average)
ax   = axes[0, 1]
kern = np.ones(7) / 7
for i in range(NNODES):
    daily    = model.nodes.newly_infected[:NTICKS, i].astype(np.float64)
    smoothed = np.convolve(daily, kern, mode="valid")
    ax.plot(np.arange(len(smoothed)) / 365, smoothed,
            color=COLORS[i], lw=0.9, label=names[i], alpha=0.85)
ax.set_xlabel("Year")
ax.set_ylabel("New infections (7-day avg)")
ax.set_title("(B) Incidence by Patch")
ax.legend(fontsize=8)

# (C) Susceptible fraction S/N with 1/R0 threshold
ax     = axes[0, 2]
sample = np.arange(0, NTICKS, 7)
for i in range(NNODES):
    S = model.nodes.S[sample, i].astype(np.float64)
    N = (S + model.nodes.E[sample, i] + model.nodes.I[sample, i]
         + model.nodes.R[sample, i]).astype(np.float64)
    ax.plot(sample / 365, S / np.maximum(N, 1),
            color=COLORS[i], lw=1, label=names[i])
ax.axhline(1.0 / R0, color="black", ls=":", lw=1.5, alpha=0.6,
           label=f"1/R0 = {1/R0:.2f}")
ax.set_xlabel("Year")
ax.set_ylabel("S / N")
ax.set_title("(C) Susceptible Fraction")
ax.legend(fontsize=7)

# (D) Population by patch over time
ax = axes[1, 0]
for i in range(NNODES):
    N_i = (model.nodes.S[:NTICKS, i] + model.nodes.E[:NTICKS, i]
           + model.nodes.I[:NTICKS, i] + model.nodes.R[:NTICKS, i]
           ).astype(np.float64)
    ax.plot(np.arange(NTICKS) / 365, N_i,
            color=COLORS[i], lw=1, label=names[i])
ax.set_xlabel("Year")
ax.set_ylabel("Population")
ax.set_title("(D) Population by Patch")
ax.legend(fontsize=8)

# (E) Effective reproduction number R_eff = R0 × S/N
ax = axes[1, 1]
for i in range(NNODES):
    S = model.nodes.S[sample, i].astype(np.float64)
    N = (S + model.nodes.E[sample, i] + model.nodes.I[sample, i]
         + model.nodes.R[sample, i]).astype(np.float64)
    ax.plot(sample / 365, R0 * S / np.maximum(N, 1),
            color=COLORS[i], lw=1, label=names[i])
ax.axhline(1.0, color="black", ls="-", lw=1.5, alpha=0.4, label="R_eff = 1")
ax.set_xlabel("Year")
ax.set_ylabel("R_eff = R0 · S/N")
ax.set_title("(E) Effective Reproduction Number")
ax.legend(fontsize=7)

# (F) Spatial heatmap — weekly incidence per patch
ax        = axes[1, 2]
nweeks    = NTICKS // 7
wkly_mat  = (model.nodes.newly_infected[:nweeks * 7]
             .reshape(nweeks, 7, NNODES).sum(axis=1))
im = ax.imshow(
    wkly_mat.T, aspect="auto", cmap="hot",
    extent=[0, nweeks / 52, -0.5, NNODES - 0.5], origin="lower",
)
ax.set_yticks(range(NNODES))
ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel("Year")
ax.set_title("(F) Weekly Incidence Heatmap")
fig.colorbar(im, ax=ax, shrink=0.8, label="Weekly infections")

plt.suptitle(
    "Pakistan Polio Corridor — 6-Patch Spatial SEIR\n"
    f"Gravity (b={GRAVITY_B}, c={GRAVITY_C}, ≤{MAX_EXPORT:.0%} FOI export)  ·  "
    "Monsoon Seasonality (1.3× Jul–Oct, 0.7× Dec–Mar)",
    fontsize=12, fontweight="bold",
)
plt.tight_layout()

outdir  = Path(__file__).parent / "eval" / "outputs"
outdir.mkdir(parents=True, exist_ok=True)
outpath = outdir / "polio_gravity_seasonal.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nDiagnostic plots saved to {outpath}")
plt.close(fig)
```

---

**Key implementation details:**

**Gravity network** (`gravity` → manual k-scaling → `row_normalizer`):
- Built with source exponent `a=1` (standard), destination `b=0.5`, decay `c=1.5`
- k-scaling normalises the mean row-sum to `GRAVITY_K=0.01` (1 % average outward coupling), then `row_normalizer(0.15)` caps every row at 15 %
- `dist_matrix[i,i] = np.inf` so `gravity()` produces exactly zero self-coupling

**Seasonal profile** (cosine, mean=1.0):
- `season(t) = 1 + 0.3·cos(2π(t−242)/365)` peaks at 1.30× on day 242 (late August, centre of Jul–Oct) and troughs at 0.70× on day ~60 (late February, within Dec–Mar)
- Cosine mean over a full period is analytically zero, so the profile mean is exactly 1.0 — the normalization step is a safety guard
- Six assertions verify the amplitude, window placement, and mean before the run

**Silent-failure guards wired in:** unit assertions on CBR/CDR, seasonal mean check, network sum and row-cap assertions both before and after the run.
