Here is the complete, self-contained model. It differs from the existing file in three ways: no seasonal forcing, no pulse importation, and `calc_capacity` is called explicitly.

```python
#!/usr/bin/env python3
"""
Republic of Eastland — 4-Patch Spatial SEIR (Respiratory Illness)

Patches     : 100k, 200k, 150k, 80k at (30N,50E),(30.5N,50.5E),(31N,51E),(31.5N,51.5E)
Disease     : R0~5, latent gamma(4,1), infectious gamma(5,2)
Init        : 90% S, 0% E, 1% I, 9% R
Components  : Susceptible, Exposed, Infectious, Recovered, Transmission,
              BirthsByCBR(CBR=30), MortalityByCDR(CDR=10)
Capacity    : calc_capacity for 10 years, safety_factor=3
Distances   : Haversine via laser.core.migration.distance
Network     : gravity(k=0.01, a=1, b=1, c=1.5) + row_normalizer(0.15)
Duration    : 10 years = 3650 ticks
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.utils import calc_capacity
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# 1. Patch Configuration
# =============================================================================

PATCHES = [
    {"nodeid": 0, "name": "Eastport",   "population": 100_000, "lat": 30.0, "lon": 50.0},
    {"nodeid": 1, "name": "Northton",   "population": 200_000, "lat": 30.5, "lon": 50.5},
    {"nodeid": 2, "name": "Centerburg", "population": 150_000, "lat": 31.0, "lon": 51.0},
    {"nodeid": 3, "name": "Westham",    "population":  80_000, "lat": 31.5, "lon": 51.5},
]

# =============================================================================
# 2. Parameters
# =============================================================================

NTICKS = 10 * 365   # 10 years
CBR    = 30.0        # per 1000/year
CDR    = 10.0        # per 1000/year

# R0 = beta * infectious_period
# infectious_period = gamma(shape=5, scale=2).mean() = 10 days
# => beta = 5 / 10 = 0.5 /day
params = PropertySet({
    "prng_seed":              42,
    "nticks":                 NTICKS,
    "beta":                   0.5,
    "capacity_safety_factor": 3.0,
})

# =============================================================================
# 3. Scenario GeoDataFrame
# =============================================================================

geometry = [Point(p["lon"], p["lat"]) for p in PATCHES]
scenario = gpd.GeoDataFrame(PATCHES, geometry=geometry, crs="EPSG:4326")

# Initial conditions: 90% S, 0% E, 1% I, 9% R
scenario["E"] = np.zeros(len(scenario), dtype=np.uint32)
scenario["I"] = np.round(0.01 * scenario["population"]).astype(np.uint32)
scenario["R"] = np.round(0.09 * scenario["population"]).astype(np.uint32)
scenario["S"] = (
    scenario["population"] - scenario["E"] - scenario["I"] - scenario["R"]
).astype(np.uint32)

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch needs initial infections"

nnodes = len(scenario)
names  = list(scenario["name"])
pops   = scenario["population"].values.astype(np.float64)

print("Initial conditions:")
for _, row in scenario.iterrows():
    print(f"  {row['name']:<12}: pop={row['population']:>7,}  "
          f"S={row['S']:>6,}  E={int(row['E']):>4}  "
          f"I={row['I']:>5,}  R={row['R']:>5,}")

# =============================================================================
# 4. Vital-Dynamics Rate Arrays  (per-1000/year — LASER divides internally)
# =============================================================================

birthrate_array = np.full((NTICKS, nnodes), CBR, dtype=np.float32)
mortality_array = np.full((NTICKS, nnodes), CDR, dtype=np.float32)

# Guard against silent unit-error (daily per-capita would be ~0.00008)
assert np.all(birthrate_array >= 1) and np.all(birthrate_array <= 60), \
    f"Birthrates must be per-1000/year"
assert np.all(mortality_array >= 1) and np.all(mortality_array <= 60), \
    f"Death rates must be per-1000/year"

# Age pyramid for newborns (exponential, ~70yr life expectancy)
ages = np.arange(100)
stable_age_dist = np.maximum((1000 * np.exp(-ages / 70.0)).astype(np.int64), 1)
pyramid = AliasedDistribution(stable_age_dist)

# =============================================================================
# 5. calc_capacity — pre-allocate agent slots for 10 years of population growth
# =============================================================================

capacity_per_node = calc_capacity(
    birthrate_array,
    scenario["population"].values,
    safety_factor=params.capacity_safety_factor,
)

print(f"\ncalc_capacity (10 yr, safety_factor={params.capacity_safety_factor}):")
for name, cap in zip(names, capacity_per_node):
    print(f"  {name:<12}: {int(cap):,}")
print(f"  {'Total':<12}: {int(capacity_per_node.sum()):,}")

# =============================================================================
# 6. Build LASER Model
#    Model.__init__ also calls calc_capacity internally using birthrate_array
#    and params.capacity_safety_factor.
# =============================================================================

model = Model(scenario, params, birthrates=birthrate_array)
print(f"\nAgent capacity allocated by Model: {model.people.capacity:,}")

# =============================================================================
# 7. Pairwise Distances — Haversine via laser.core.migration.distance  (km)
# =============================================================================

lats = np.array([p["lat"] for p in PATCHES])
lons = np.array([p["lon"] for p in PATCHES])

dist_matrix = np.zeros((nnodes, nnodes))
for i in range(nnodes):
    for j in range(nnodes):
        if i != j:
            dist_matrix[i, j] = distance(lats[i], lons[i], lats[j], lons[j])

print("\nPairwise distances (km):")
print(f"  {'':>12}  " + "  ".join(f"{n:>11}" for n in names))
for i in range(nnodes):
    row = "  ".join(f"{dist_matrix[i, j]:>11.1f}" for j in range(nnodes))
    print(f"  {names[i]:>12}  {row}")

# =============================================================================
# 8. Gravity Network   M_{i,j} = k * pop_i^a * pop_j^b / dist_{i,j}^c
#    k=0.01, a=1, b=1, c=1.5
#
#    row_normalizer(0.15) caps each row's sum at 0.15, converting the raw
#    gravity values into fractional FOI-transfer rates usable by TransmissionSE.
# =============================================================================

network = gravity(pops, dist_matrix, 0.01, 1, 1, 1.5)
network = row_normalizer(network, 0.15)
model.network = network

assert model.network.sum() > 0, \
    "Network is all zeros — check k, populations, or distances"
assert model.network.sum(axis=1).max() < 0.3, \
    f"Network row sum {model.network.sum(axis=1).max():.4f} exceeds 0.3 safety limit"

row_sums = model.network.sum(axis=1)
print(f"\nGravity network (k=0.01, a=1, b=1, c=1.5; cap=0.15):")
for i, name in enumerate(names):
    print(f"  {name:<12}: outward fraction = {row_sums[i]:.6f}")

# =============================================================================
# 9. Duration Distributions
#    Latent     gamma(shape=4, scale=1) → mean = 4 days
#    Infectious gamma(shape=5, scale=2) → mean = 10 days
# =============================================================================

expdurdist = dists.gamma(shape=4, scale=1.0)
infdurdist = dists.gamma(shape=5, scale=2.0)

# =============================================================================
# 10. Assemble Components
#
#  Execution order per tick:
#    Susceptible  — propagates S[t+1] = S[t]
#    Exposed      — decrements etimer; fires E→I when etimer = 0
#    Infectious   — decrements itimer; fires I→R when itimer = 0
#    Recovered    — propagates R[t+1] = R[t]
#    Transmission — computes FOI with gravity coupling; drives S→E
#    BirthsByCBR  — Poisson births at CBR=30/1000/yr; newborns enter S
#    MortalityByCDR — removes agents at CDR=10/1000/yr from all compartments
# =============================================================================

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist),
    BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortality_array),
]

# =============================================================================
# 11. Run — 10 years (3650 ticks)
# =============================================================================

print(f"\nRunning 10-year Eastland respiratory SEIR ...")
print(f"  Patches: {nnodes}  ({', '.join(names)})")
print(f"  Total population: {int(pops.sum()):,}")
print(f"  R0=5  beta=0.5  latent=gamma(4,1)  infectious=gamma(5,2)")
print(f"  CBR={CBR}/1000/yr  CDR={CDR}/1000/yr")

model.run("Eastland Respiratory SEIR")
print("Simulation complete.")

# =============================================================================
# 12. Post-Run Verification
# =============================================================================

print("\n--- Verification ---")

# Population trajectory (CBR > CDR → expect ~20% growth)
pop_start = int(model.nodes.S[0].sum() + model.nodes.E[0].sum()
                + model.nodes.I[0].sum() + model.nodes.R[0].sum())
pop_end   = int(model.nodes.S[NTICKS-1].sum() + model.nodes.E[NTICKS-1].sum()
                + model.nodes.I[NTICKS-1].sum() + model.nodes.R[NTICKS-1].sum())
growth_pct = (pop_end / pop_start - 1) * 100
print(f"Population: {pop_start:,} → {pop_end:,}  ({growth_pct:+.1f}%)")
if pop_end > pop_start:
    print("  Population growing: OK (CBR > CDR)")
else:
    print("  WARNING: Population not growing — check birth/death rate units!")

# Compartment non-negativity
neg_flag = False
for arr, label in [(model.nodes.S, "S"), (model.nodes.E, "E"),
                   (model.nodes.I, "I"), (model.nodes.R, "R")]:
    if np.any(arr[:NTICKS] < 0):
        print(f"  WARNING: Negative {label} values detected!")
        neg_flag = True
if not neg_flag:
    print("  All compartments non-negative: OK")

# Epidemic dynamics and spatial coupling
total_inf = int(model.nodes.newly_infected[:NTICKS].sum())
patch_inf = (model.nodes.newly_infected[:NTICKS].sum(axis=0) > 0).sum()
print(f"  Total infections (10 yr): {total_inf:,}")
print(f"  Patches with infections:  {int(patch_inf)}/{nnodes}")

# Compartment table
print(f"\nNational compartment counts:")
print(f"  {'Tick':<8} {'S':>10} {'E':>8} {'I':>8} {'R':>10} {'N':>10}")
print(f"  {'-' * 56}")
for tick in [0, NTICKS // 2, NTICKS - 1]:
    S = int(model.nodes.S[tick].sum())
    E = int(model.nodes.E[tick].sum())
    I = int(model.nodes.I[tick].sum())
    R = int(model.nodes.R[tick].sum())
    print(f"  {tick:<8} {S:>10,} {E:>8,} {I:>8,} {R:>10,} {S+E+I+R:>10,}")

# Per-patch annual incidence
nyears = NTICKS // 365
print(f"\nMean annual incidence by patch:")
print(f"  {'Patch':<12} {'Population':>11} {'Cases/yr':>10} {'Per 100k':>10}")
print(f"  {'-' * 47}")
for i, name in enumerate(names):
    annual = np.array([
        model.nodes.newly_infected[y*365:(y+1)*365, i].sum()
        for y in range(nyears)
    ])
    avg = annual.mean()
    print(f"  {name:<12} {int(pops[i]):>11,} {avg:>10,.0f} {avg/pops[i]*100_000:>10,.0f}")

# =============================================================================
# 13. Diagnostic Plots
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 9))
colors = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12"]

# (A) Incidence by patch — 7-day rolling mean
ax = axes[0, 0]
kernel = np.ones(7) / 7
for i in range(nnodes):
    daily    = model.nodes.newly_infected[:NTICKS, i].astype(np.float64)
    smoothed = np.convolve(daily, kernel, mode="valid")
    ax.plot(np.arange(len(smoothed)) / 365, smoothed,
            color=colors[i], lw=0.9, label=names[i], alpha=0.85)
ax.set_xlabel("Year"); ax.set_ylabel("New infections (7-day avg)")
ax.set_title("(A) Incidence by Patch"); ax.legend(fontsize=8)

# (B) Susceptible fraction S/N
ax = axes[0, 1]
sample = np.arange(0, NTICKS, 7)
for i in range(nnodes):
    S = model.nodes.S[sample, i].astype(np.float64)
    N = np.maximum(S + model.nodes.E[sample, i]
                     + model.nodes.I[sample, i]
                     + model.nodes.R[sample, i], 1.0)
    ax.plot(sample / 365, S / N, color=colors[i], lw=1.0, label=names[i])
ax.axhline(1/5, color="black", ls=":", lw=1.5, alpha=0.6, label="1/R0 = 0.2")
ax.set_xlabel("Year"); ax.set_ylabel("S / N")
ax.set_title("(B) Susceptible Fraction"); ax.legend(fontsize=7)

# (C) Population trajectory
ax = axes[1, 0]
for i in range(nnodes):
    N_i = (model.nodes.S[:NTICKS, i].astype(np.float64)
           + model.nodes.E[:NTICKS, i]
           + model.nodes.I[:NTICKS, i]
           + model.nodes.R[:NTICKS, i])
    ax.plot(np.arange(NTICKS) / 365, N_i, color=colors[i], lw=1.0, label=names[i])
ax.set_xlabel("Year"); ax.set_ylabel("Population")
ax.set_title("(C) Population by Patch"); ax.legend(fontsize=8)

# (D) Weekly incidence heatmap
ax = axes[1, 1]
nweeks = NTICKS // 7
weekly = model.nodes.newly_infected[:nweeks*7].reshape(nweeks, 7, nnodes).sum(axis=1)
im = ax.imshow(weekly.T, aspect="auto", cmap="hot",
               extent=[0, nweeks/52, -0.5, nnodes-0.5], origin="lower")
ax.set_yticks(range(nnodes)); ax.set_yticklabels(names)
ax.set_xlabel("Year"); ax.set_title("(D) Weekly Incidence Heatmap")
fig.colorbar(im, ax=ax, shrink=0.8, label="Weekly infections")

plt.suptitle("Republic of Eastland — Respiratory SEIR  (R0=5, 10 yr)",
             fontsize=13, fontweight="bold")
plt.tight_layout()

outdir = Path(__file__).parent / "eval" / "outputs"
outdir.mkdir(parents=True, exist_ok=True)
outpath = outdir / "eastland_respiratory_seir.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nDiagnostic plot saved to: {outpath}")
plt.close(fig)
```

---

Key design decisions, annotated:

**Beta derivation** — `beta = R0 / infectious_period = 5 / (5×2) = 0.5`. The gamma(5,2) mean is 10 days. Beta is set directly on `params`; LASER uses it inside `TransmissionSE`.

**`calc_capacity`** (step 5) — called explicitly with `birthrate_array` (shape `NTICKS×nnodes`) and `safety_factor=3`. `Model.__init__` then calls it again internally; both are idempotent and the result is the same.

**Distance matrix** (step 7) — built by calling `distance(lat_i, lon_i, lat_j, lon_j)` for each off-diagonal pair. Diagonal stays 0; the `gravity()` function treats zero-distance (self-coupling) as zero flow.

**Gravity normalization** (step 8) — the raw values from `gravity(pops, dists, 0.01, 1, 1, 1.5)` are in absolute units (population × population / km^1.5) and would produce FOI-transfer fractions in the millions. `row_normalizer(network, 0.15)` rescales each row so no patch exports more than 15% of its local FOI — the standard LASER practice for turning gravity into fractional coupling.

**Component order** (step 10) — `Susceptible` and `Recovered` must bracket the transition steps so node-level counts satisfy S+E+I+R=N every tick. `Transmission` goes after `Recovered` so it operates on the fully-propagated state. Vital dynamics go last.
