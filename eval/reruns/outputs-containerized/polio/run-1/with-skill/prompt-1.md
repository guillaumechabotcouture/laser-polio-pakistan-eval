The filesystem is read-only, so I'll deliver the complete code here. It closely follows the established pattern from `eastland_respiratory_model.py` and the existing `scripts/polio_seir_basic_10patch.py`.

---

**Key parameter derivation:**

| Parameter | Value | Reasoning |
|---|---|---|
| `beta` | `6/28 ≈ 0.2143` | `R0 = beta × D_I → beta = R0 / 28` |
| `expdurdist` | `gamma(shape=3, scale=1.0)` | mean = 3 days (latent) |
| `infdurdist` | `normal(loc=28, scale=3)` | mean = 28 days (infectious) |
| `S₀` per patch | `4,995` | `100k − 5 − 95k` |
| `R_eff(t=0)` | `≈ 0.30` | `6 × 4995/100k < 1` — births replenish susceptibles |

---

```python
#!/usr/bin/env python3
"""
Pakistan Polio — 10-Patch Spatial SEIR Model

Spatial SEIR simulation for WPV1 transmission dynamics across 10 districts
using the LASER framework (laser-generic v1.0.0).

Specification:
    - 10 patches (districts), 100,000 population each
    - R0 = 6,  beta = R0 / D_I = 6 / 28 ≈ 0.2143
    - Latent period:     gamma(shape=3, scale=1.0)   → mean =  3 days
    - Infectious period: normal(loc=28, scale=3)      → mean = 28 days
    - Initial per patch: S=4,995  E=0  I=5  R=95,000
    - Vital dynamics:    CBR=29, CDR=7  per 1,000/year (Pakistan)
    - Seasonal forcing:  monsoon peak day=245 (+/-30%), ValuesMap
    - Gravity network:   k=0.005, a=0, b=0.5, c=1.5, row_normalizer(0.2)
    - Importation:       2 agents every 30 ticks in 3 endemic patches
    - Duration:          10 years (3,650 ticks)

Usage:
    python polio_pakistan_seir.py

Notes:
    R_eff(t=0) = R0 × S/N = 6 × 0.04995 ≈ 0.30 < 1.
    Disease cannot spread immediately; births (CBR=29) add ~2,900 susceptibles
    per patch per year. R_eff crosses 1 when S/N > 1/R0 ≈ 16.7%,
    around year 3-4. Importation in 3 endemic patches sustains viral presence
    until the susceptible pool is large enough to support epidemic cycles.
"""

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
from laser.core.migration import gravity, row_normalizer, distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================================
# Importation Component — Endemic Seeding in Selected Patches
#
# Infects susceptible agents in specified patches periodically to model
# cross-border importation and persistent transmission foci in endemic areas.
# Only targets susceptible agents (unlike built-in Infect_Random_Agents).
# ============================================================================

class PatchImportation:
    """Periodic importation of poliovirus into susceptible agents in endemic patches.

    Every `period` ticks, selects up to `count` susceptible agents per patch
    from `patch_ids` and transitions them to INFECTIOUS. Infecting only
    susceptibles is epidemiologically precise and prevents wasting importation
    events on already-infectious or recovered agents.

    Parameters:
        model:     LASER Model instance
        infdurdist: Distribution for infectious duration sampling
        patch_ids: List of patch indices to seed (default: all patches)
        period:    Ticks between importation events (default: 30)
        count:     Agents infected per patch per event (default: 2)
        end_tick:  Stop importation after this tick (default: nticks)
    """

    def __init__(self, model, infdurdist, patch_ids=None,
                 period=30, count=2, end_tick=None):
        self.model     = model
        self.infdurdist = infdurdist
        self.patch_ids  = patch_ids if patch_ids is not None else list(
            range(model.nodes.count)
        )
        self.period   = period
        self.count    = count
        self.end_tick = end_tick if end_tick is not None else model.params.nticks

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0 or tick >= self.end_tick:
            return

        people = self.model.people
        nodes  = self.model.nodes
        active = people.count

        for patch in self.patch_ids:
            # Susceptible agents in this patch
            in_patch = np.nonzero(
                (people.state[:active] == SEIR.State.SUSCEPTIBLE.value) &
                (people.nodeid[:active] == patch)
            )[0]
            if len(in_patch) == 0:
                continue

            n_infect = min(self.count, len(in_patch))
            chosen   = np.random.choice(in_patch, size=n_infect, replace=False)

            # Transition to INFECTIOUS
            people.state[chosen] = SEIR.State.INFECTIOUS.value

            # Sample infectious durations
            offsets = np.zeros(n_infect, dtype=np.float32)
            samples = dists.sample_floats(self.infdurdist, offsets)
            samples = np.maximum(np.round(samples), 1).astype(people.itimer.dtype)
            people.itimer[chosen] = samples

            # Update node-level compartment counts
            nodes.S[tick + 1, patch] = max(
                nodes.S[tick + 1, patch] - n_infect, 0
            )
            nodes.I[tick + 1, patch] += n_infect


# ============================================================================
# Patch Configuration — 10 Pakistan Districts
# Approximate coordinates for representative urban/district centres.
# ============================================================================

patches = [
    {"nodeid": 0, "name": "District_01", "population": 100_000, "lat": 34.0,  "lon": 71.5},
    {"nodeid": 1, "name": "District_02", "population": 100_000, "lat": 33.7,  "lon": 73.0},
    {"nodeid": 2, "name": "District_03", "population": 100_000, "lat": 31.5,  "lon": 74.4},
    {"nodeid": 3, "name": "District_04", "population": 100_000, "lat": 30.2,  "lon": 67.0},
    {"nodeid": 4, "name": "District_05", "population": 100_000, "lat": 31.0,  "lon": 70.5},
    {"nodeid": 5, "name": "District_06", "population": 100_000, "lat": 30.5,  "lon": 72.0},
    {"nodeid": 6, "name": "District_07", "population": 100_000, "lat": 31.5,  "lon": 73.1},
    {"nodeid": 7, "name": "District_08", "population": 100_000, "lat": 27.5,  "lon": 65.5},
    {"nodeid": 8, "name": "District_09", "population": 100_000, "lat": 25.4,  "lon": 68.4},
    {"nodeid": 9, "name": "District_10", "population": 100_000, "lat": 24.9,  "lon": 67.0},
]

# ============================================================================
# Simulation Parameters
# ============================================================================

NTICKS = 10 * 365      # 10-year simulation
CBR    = 29.0           # Crude birth rate per 1,000/year (Pakistan ~29)
CDR    =  7.0           # Crude death rate per 1,000/year

# beta = R0 / D_I = 6 / 28 ≈ 0.2143
# LASER FOI: lambda_i = beta * seasonality[tick] * I_i / N_i
# R0 = beta * D_I = 0.2143 * 28 = 6.0
BETA = 6.0 / 28.0
R0   = 6.0

params = PropertySet({
    "prng_seed":              42,
    "nticks":                 NTICKS,
    "beta":                   BETA,        # ~0.2143 — R0=6, D_I=28 days
    "capacity_safety_factor": 2.5,         # Pre-allocate 2.5x initial pop
})

print("=== Pakistan Polio SEIR — 10 Patch Spatial Model ===")
print(f"  beta           = {BETA:.6f}  (R0 = beta x D_I = {BETA*28:.2f})")
print(f"  Latent period  = 3 days  [gamma(shape=3, scale=1.0), mean=3d, std≈1.7d]")
print(f"  Infect period  = 28 days [normal(loc=28, scale=3), mean=28d, std=3d]")
print(f"  CBR/CDR        = {CBR}/{CDR} per 1,000/year")
print(f"  Simulation     = {NTICKS} ticks ({NTICKS // 365} years)")

# ============================================================================
# Build Scenario (GeoDataFrame)
#
# Initial conditions per patch:
#   population = 100,000
#   I = 5            (5 infectious)
#   R = 95,000       (95% recovered/immune)
#   E = 0
#   S = 100,000 - 5 - 0 - 95,000 = 4,995
#
# R_eff(t=0) = R0 * S/N = 6 * 4995/100000 ≈ 0.30 < 1
# Disease cannot spread — births replenish susceptibles over time.
# ============================================================================

geometry = [Point(p["lon"], p["lat"]) for p in patches]
scenario = gpd.GeoDataFrame(patches, geometry=geometry, crs="EPSG:4326")
nnodes   = len(scenario)

scenario["I"] = np.full(nnodes, 5, dtype=np.uint32)
scenario["E"] = np.zeros(nnodes, dtype=np.uint32)
scenario["R"] = np.full(nnodes, int(0.95 * 100_000), dtype=np.uint32)
scenario["S"] = (scenario["population"]
                 - scenario["E"] - scenario["I"] - scenario["R"]).astype(np.uint32)

# Validate initial conditions
assert (scenario.S + scenario.E + scenario.I + scenario.R
        == scenario.population).all(), \
    "S + E + I + R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch must have initial infections"

S0 = int(scenario["S"].iloc[0])
print(f"\nInitial per patch: S={S0:,}  E=0  I=5  R=95,000")
print(f"  R_eff(t=0) = R0 x S/N = {R0} x {S0/100_000:.5f} = {R0*S0/100_000:.3f} < 1")

# ============================================================================
# Vital Dynamics — Birth and Death Rates
# Units: per 1,000 per year (LASER divides by 1,000 internally)
# ============================================================================

birthrate_array = np.full((NTICKS, nnodes), CBR, dtype=np.float32)
mortality_array = np.full((NTICKS, nnodes), CDR, dtype=np.float32)

# Guard: must be per-1000/year — passing daily per-capita is the #1 silent error
assert np.all(birthrate_array >= 1) and np.all(birthrate_array <= 60), \
    f"Birthrates must be per-1000/year, got {birthrate_array.min():.6f}"
assert np.all(mortality_array >= 1) and np.all(mortality_array <= 60), \
    f"Mortality rates must be per-1000/year, got {mortality_array.min():.6f}"

# Stable age distribution for newborn initialisation (~65-year life expectancy)
ages           = np.arange(100)
stable_age_dist = np.array(1000 * np.exp(-ages / 65.0), dtype=np.int64)
stable_age_dist = np.maximum(stable_age_dist, 1)
pyramid        = AliasedDistribution(stable_age_dist)

# ============================================================================
# Seasonal Forcing
#
# Pakistan polio peaks during/after monsoon season (Jul-Oct) driven by
# WASH disruption and population displacement. Model a monsoon peak centred
# on day 245 (early September) with cosine amplitude ±30%.
#
#   season(t) = 1.0 + 0.30 * cos(2π*(t - 245)/365)
#   Peak (day 245, early Sep): 1.30   Trough (winter): 0.70
#   Annual mean: 1.0 (normalized)
# ============================================================================

days       = np.arange(365)
peak_day   = 245      # Early September (monsoon peak)
amplitude  = 0.30
season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
season_365 /= season_365.mean()    # Normalize to mean = 1.0

assert abs(season_365.mean() - 1.0) < 0.01, \
    f"Seasonal profile mean={season_365.mean():.4f}, must be ~1.0"

print(f"\nSeasonal forcing (monsoon peak, day={peak_day}):")
print(f"  Peak:   {season_365.max():.3f}x  (day {int(np.argmax(season_365))})")
print(f"  Trough: {season_365.min():.3f}x  (day {int(np.argmin(season_365))})")
print(f"  Mean:   {season_365.mean():.4f}")

# Tile across 10-year simulation; ValuesMap broadcasts to all nodes
season_tiled = np.tile(season_365, NTICKS // 365 + 1)[:NTICKS]
seasonality  = ValuesMap.from_timeseries(season_tiled, nnodes)

# ============================================================================
# Build Model
#
# Model.__init__ calls calc_capacity(birthrates, initial_pop, safety_factor)
# to pre-allocate agent slots. capacity_safety_factor=2.5 ensures enough
# headroom for population growth over 10 years (CBR=29, CDR=7 → net ~2.2%/yr).
# ============================================================================

model = Model(scenario, params, birthrates=birthrate_array)

# ============================================================================
# Pairwise Haversine Distance Matrix
# ============================================================================

lats = np.array([p["lat"] for p in patches])
lons = np.array([p["lon"] for p in patches])

dist_matrix = np.zeros((nnodes, nnodes))
for i in range(nnodes):
    for j in range(nnodes):
        if i != j:
            dist_matrix[i, j] = distance(lats[i], lons[i], lats[j], lons[j])

print(f"\nSelected inter-patch distances (km):")
print(f"  D01-D02 (Peshawar–Islamabad):  {dist_matrix[0, 1]:>5.0f} km")
print(f"  D02-D03 (Islamabad–Lahore):    {dist_matrix[1, 2]:>5.0f} km")
print(f"  D09-D10 (Hyderabad–Karachi):   {dist_matrix[8, 9]:>5.0f} km")

# ============================================================================
# Gravity Migration Network
#
# Convention in this project: a=0 (source pop does not drive outward flow)
# M_{i,j} = k * p_j^b / d_{ij}^c  with k=1, then normalise so mean
# outward fraction equals gravity_k. Cap at 20% per row via row_normalizer.
#
# Parameters: gravity_k=0.005, b=0.5, c=1.5
# ============================================================================

GRAVITY_K = 0.005
GRAVITY_B = 0.5
GRAVITY_C = 1.5

pops    = np.array(scenario.population, dtype=np.float64)
network = gravity(pops, dist_matrix, k=1, a=0, b=GRAVITY_B, c=GRAVITY_C)

# Normalise so gravity_k represents average outward flow fraction
avg_export = np.mean(network.sum(axis=1))
if avg_export > 0:
    network = network / avg_export * GRAVITY_K

network      = row_normalizer(network, max_fraction=0.2)
model.network = network

# Validate
assert model.network.sum() > 0, \
    "Gravity network is all zeros — check gravity parameters"
row_sums = model.network.sum(axis=1)
assert row_sums.max() < 0.3, \
    f"Max row sum {row_sums.max():.4f} exceeds safety limit"

print(f"\nGravity network (k={GRAVITY_K}, b={GRAVITY_B}, c={GRAVITY_C}):")
print(f"  Network total: {model.network.sum():.6f}")
print(f"  Row sum range: [{row_sums.min():.6f}, {row_sums.max():.6f}]")

# ============================================================================
# Duration Distributions
#
# expdurdist — Exposed → Infectious transition timer
#   gamma(shape=3, scale=1.0): mean = 3 days, std ≈ 1.73 days
#
# infdurdist — Infectious → Recovered transition timer
#   normal(loc=28, scale=3): mean = 28 days, std = 3 days
#   infdurmin=1 clamps rare negative samples to 1 tick.
# ============================================================================

expdurdist = dists.gamma(shape=3, scale=1.0)    # Latent:     mean=3d
infdurdist = dists.normal(loc=28, scale=3.0)    # Infectious: mean=28d

# ============================================================================
# Endemic Importation Patches
#
# R_eff(t=0) ≈ 0.30 < 1 — disease fades without importation.
# Model persistent foci in 3 geographically distinct patches representing
# endemic corridors (e.g., KP border, Balochistan corridor, Sindh).
# 2 infections seeded every 30 ticks in each endemic patch.
# ============================================================================

ENDEMIC_PATCHES = [0, 4, 9]   # Patches 1, 5, 10 — endemic corridors

# ============================================================================
# Assemble Components
#
# Execution order per tick:
#   1. Susceptible    — propagate S[t] → S[t+1]  (count bookkeeping)
#   2. Exposed        — decrement etimer; E→I transitions
#   3. Infectious     — decrement itimer; I→R transitions
#   4. Recovered      — propagate R[t] → R[t+1]  (count bookkeeping)
#   5. BirthsByCBR    — add newborns (CBR=29/1000/yr)
#   6. MortalityByCDR — remove deaths (CDR=7/1000/yr)
#   7. PatchImportation — seed 2 susceptibles in endemic patches every 30 ticks
#   8. Transmission   — compute FOI with gravity network + seasonality → S→E
#
# Births/deaths before importation so population counts are current.
# Transmission last so it operates on updated compartment sizes.
# ============================================================================

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortality_array),
    PatchImportation(model, infdurdist, patch_ids=ENDEMIC_PATCHES,
                     period=30, count=2),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
]

# ============================================================================
# Run Simulation
# ============================================================================

total_pop = scenario.population.sum()
print(f"\nRunning simulation...")
print(f"  Patches:         {nnodes} districts, {total_pop:,} total population")
print(f"  R0:              {R0}  (beta={BETA:.4f}, D_I=28d)")
print(f"  Initial S/N:     {S0/100_000:.5f}  →  R_eff={R0*S0/100_000:.3f}")
print(f"  Importation:     2 agents/30 ticks in patches {ENDEMIC_PATCHES}")
print(f"  Agent capacity:  {model.people.capacity:,}")

model.run("Pakistan Polio SEIR 10-Patch")
print("Simulation complete.")

# ============================================================================
# Post-Run Verification
# ============================================================================

print("\n--- Post-Run Verification ---")

# 1. Population trajectory — must grow (CBR > CDR)
pop_start = (model.nodes.S[1].sum()  + model.nodes.E[1].sum()
             + model.nodes.I[1].sum() + model.nodes.R[1].sum())
pop_end   = (model.nodes.S[NTICKS - 1].sum()  + model.nodes.E[NTICKS - 1].sum()
             + model.nodes.I[NTICKS - 1].sum() + model.nodes.R[NTICKS - 1].sum())
ratio = pop_end / max(pop_start, 1)
print(f"  Population: {int(pop_start):,} → {int(pop_end):,}  (ratio={ratio:.3f})")
if abs(ratio - 1.0) < 0.001:
    print("  WARNING: Population static — check birthrate units (must be per-1000/yr)")
else:
    print("  Population growing: OK (CBR > CDR)")

# 2. Compartment non-negativity
neg_found = False
for attr in ["S", "E", "I", "R"]:
    arr = getattr(model.nodes, attr)[:NTICKS]
    if np.any(arr < 0):
        print(f"  WARNING: Negative {attr} values detected!")
        neg_found = True
if not neg_found:
    print("  All compartments non-negative: OK")

# 3. Epidemic dynamics
total_infections = int(model.nodes.newly_infected[:NTICKS].sum())
print(f"  Total infections (10 yr): {total_infections:,}")
if total_infections == 0:
    print("  WARNING: Zero incidence — check beta, initial I, and importation")

# 4. Spatial coupling
patch_totals   = model.nodes.newly_infected[:NTICKS].sum(axis=0)
active_patches = int(np.sum(patch_totals > 0))
print(f"  Patches with infections: {active_patches}/{nnodes}")

# 5. Network check
net_sum = model.network.sum()
max_row = model.network.sum(axis=1).max()
print(f"  Network sum={net_sum:.6f}, max_row={max_row:.6f}: "
      + ("OK" if net_sum > 0 and max_row < 0.3 else "WARNING"))

# 6. Compartment table at key timepoints
print(f"\n  Compartment counts (national totals):")
print(f"  {'Tick':<8} {'Year':<6} {'S':>10} {'E':>6} {'I':>8} {'R':>10} {'N':>10}")
print(f"  {'-' * 60}")
for tick in [0, NTICKS // 4, NTICKS // 2, 3 * NTICKS // 4, NTICKS - 1]:
    S  = int(model.nodes.S[tick].sum())
    E  = int(model.nodes.E[tick].sum())
    I  = int(model.nodes.I[tick].sum())
    R  = int(model.nodes.R[tick].sum())
    yr = tick / 365
    print(f"  {tick:<8} {yr:<6.1f} {S:>10,} {E:>6,} {I:>8,} {R:>10,} "
          f"{S+E+I+R:>10,}")

# 7. Annual incidence summary
print(f"\n  Annual infections by patch:")
print(f"  {'District':<12} {'Total (10yr)':>12} {'Avg/yr':>10} {'per 100k':>10}")
print(f"  {'-' * 46}")
names = [p["name"] for p in patches]
for i in range(nnodes):
    yearly = [int(model.nodes.newly_infected[y*365:(y+1)*365, i].sum())
              for y in range(NTICKS // 365)]
    total  = sum(yearly)
    avg    = np.mean(yearly)
    print(f"  {names[i]:<12} {total:>12,} {avg:>10,.1f} {avg:>10,.1f}")

# 8. R_eff in years 5-10
sample_late = np.arange(5 * 365, NTICKS, 30)
S_nat = model.nodes.S[sample_late].sum(axis=1).astype(np.float64)
N_nat = (model.nodes.S[sample_late].sum(axis=1)
         + model.nodes.E[sample_late].sum(axis=1)
         + model.nodes.I[sample_late].sum(axis=1)
         + model.nodes.R[sample_late].sum(axis=1)).astype(np.float64)
Reff_late = R0 * S_nat / np.maximum(N_nat, 1)
print(f"\n  R_eff (years 5-10): mean={Reff_late.mean():.3f}, "
      f"range=[{Reff_late.min():.3f}, {Reff_late.max():.3f}]")

# ============================================================================
# Diagnostic Plots — 4-panel figure
# ============================================================================

cmap   = plt.cm.tab10
colors = [cmap(i / nnodes) for i in range(nnodes)]

fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# (A) Monsoon seasonal forcing profile
ax = axes[0, 0]
ax.plot(days, season_365, "b-", lw=2)
ax.axhline(1.0, color="gray", ls="--", alpha=0.5, label="Mean = 1.0")
ax.axvspan(182, 304, alpha=0.12, color="blue", label="Monsoon (Jul–Oct)")
ax.axhline(season_365.max(), color="red",  ls=":", alpha=0.6,
           label=f"Peak = {season_365.max():.2f}x")
ax.axhline(season_365.min(), color="navy", ls=":", alpha=0.6,
           label=f"Trough = {season_365.min():.2f}x")
ax.set_xlabel("Day of Year")
ax.set_ylabel("Seasonal Multiplier")
ax.set_title("(A) Monsoon Seasonal Forcing")
ax.legend(fontsize=8)

# (B) Weekly incidence by patch
ax = axes[0, 1]
nweeks = NTICKS // 7
for i in range(nnodes):
    weekly  = model.nodes.newly_infected[:nweeks*7, i].reshape(nweeks, 7).sum(axis=1)
    weeks_x = np.arange(nweeks) / 52.0
    ax.plot(weeks_x, weekly, color=colors[i], lw=0.8, alpha=0.85, label=names[i])
ax.set_xlabel("Year")
ax.set_ylabel("Weekly Infections")
ax.set_title("(B) Weekly Incidence by Patch")
ax.legend(fontsize=7, ncol=2, loc="upper right")

# (C) Susceptible fraction over time
ax = axes[1, 0]
sample30 = np.arange(0, NTICKS, 30)
time30   = sample30 / 365.0
for i in range(nnodes):
    S = model.nodes.S[sample30, i].astype(np.float64)
    N = (model.nodes.S[sample30, i] + model.nodes.E[sample30, i]
         + model.nodes.I[sample30, i] + model.nodes.R[sample30, i]).astype(np.float64)
    N = np.maximum(N, 1.0)
    ax.plot(time30, S / N, color=colors[i], lw=0.8, alpha=0.85, label=names[i])
ax.axhline(1.0 / R0, color="black", ls="--", lw=1.5,
           label=f"S* = 1/R0 = {1/R0:.3f}")
ax.set_xlabel("Year")
ax.set_ylabel("S / N")
ax.set_title("(C) Susceptible Fraction by Patch")
ax.legend(fontsize=7, ncol=2, loc="upper right")

# (D) Effective reproduction number R_eff = R0 × S/N
ax = axes[1, 1]
for i in range(nnodes):
    S = model.nodes.S[sample30, i].astype(np.float64)
    N = (model.nodes.S[sample30, i] + model.nodes.E[sample30, i]
         + model.nodes.I[sample30, i] + model.nodes.R[sample30, i]).astype(np.float64)
    N = np.maximum(N, 1.0)
    ax.plot(time30, R0 * S / N, color=colors[i], lw=0.8, alpha=0.85)
ax.axhline(1.0, color="black", ls="-", lw=2, alpha=0.5, label="R_eff = 1")
ax.set_xlabel("Year")
ax.set_ylabel("R_eff = R0 × S/N")
ax.set_title("(D) Effective Reproduction Number")
ax.set_ylim(0, 4)
ax.legend(fontsize=9)

plt.suptitle(
    f"Pakistan Polio SEIR — 10 Districts, R0=6, D_I=28d, 10-Year Simulation",
    fontsize=13, fontweight="bold",
)
plt.tight_layout()

outdir  = Path(__file__).parent / "eval" / "outputs"
outdir.mkdir(parents=True, exist_ok=True)
outpath = outdir / "polio_pakistan_seir.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nDiagnostic plot saved to {outpath}")
plt.close(fig)
```

---

**What the code does section-by-section:**

**Parameters**
- `beta = 6/28 ≈ 0.2143` from `R0 = beta × D_I`
- `expdurdist = gamma(3, 1.0)` → mean 3-day latent
- `infdurdist = normal(28, 3)` → mean 28-day infectious
- `infdurmin=1` (built into `Exposed`/`Infectious`) clamps rare near-zero samples

**Initial conditions** — `S=4,995 / E=0 / I=5 / R=95,000` per patch. `R_eff(t=0) ≈ 0.30`, so births drive the long-run dynamics: at CBR=29, susceptibles grow by ~2,900/year/patch; R_eff crosses 1 around year 3–4.

**Vital dynamics** — `BirthsByCBR(CBR=29)` and `MortalityByCDR(CDR=7)`, both in per-1,000/year with the assertion guard to prevent the #1 silent failure.

**Gravity network** — `k=0.005, a=0, b=0.5, c=1.5` normalised so `k` represents average outward fraction, then capped at 20% with `row_normalizer`.

**Seasonal forcing** — Monsoon peak day 245 (early September), ±30% cosine profile normalised to mean=1.0, wrapped in `ValuesMap.from_timeseries` and passed to `SEIR.Transmission`.

**Importation** — `PatchImportation` seeds 2 susceptible agents every 30 ticks in patches 0, 4, 9 (representing endemic KP/Balochistan/Sindh corridors), sustaining viral presence until susceptibles accumulate enough for self-sustaining transmission.

**Plots** — Four diagnostic panels: seasonal profile, weekly incidence by patch, susceptible fraction (with `S* = 1/R0 = 0.167` threshold line), and R_eff over time.
