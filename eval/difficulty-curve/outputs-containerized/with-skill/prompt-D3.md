The model runs cleanly. Here is the complete code — since `/workspace` is mounted read-only, paste this into `eastland_respiratory_model.py`:

```python
"""
Spatial SEIR model for a respiratory illness in the Republic of Eastland.

Specification
─────────────
  4 patches:           Eastport (100k), Northton (200k), Centerburg (150k), Westham (80k)
  R0:                  5
  Latent period:       gamma(shape=4, scale=1),  mean = 4 days
  Infectious period:   gamma(shape=5, scale=2),  mean = 10 days
  Initial conditions:  90% S, 0% E, 1% I, 9% R
  Component order:     Susceptible → Transmission → Exposed → Infectious → Recovered
  Transmission:        expdurdist = infdurdist (gamma infectious-period dist, per spec)
  Duration:            1 year (365 ticks, daily time-step)
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# 1. Patch definitions
# ─────────────────────────────────────────────────────────────────────────────

patches = [
    {"nodeid": 0, "name": "Eastport",   "population": 100_000, "lat": 30.0, "lon": 50.0},
    {"nodeid": 1, "name": "Northton",   "population": 200_000, "lat": 30.5, "lon": 50.5},
    {"nodeid": 2, "name": "Centerburg", "population": 150_000, "lat": 31.0, "lon": 51.0},
    {"nodeid": 3, "name": "Westham",    "population":  80_000, "lat": 31.5, "lon": 51.5},
]

NTICKS = 365
NNODES = len(patches)
populations = np.array([p["population"] for p in patches], dtype=np.int32)


# ─────────────────────────────────────────────────────────────────────────────
# 2. Initial conditions: 90% S, 0% E, 1% I, 9% R
# ─────────────────────────────────────────────────────────────────────────────

I_init = np.round(0.01 * populations).astype(np.int32)
R_init = np.round(0.09 * populations).astype(np.int32)
E_init = np.zeros(NNODES, dtype=np.int32)
S_init = (populations - I_init - R_init - E_init).astype(np.int32)

scenario = gpd.GeoDataFrame(
    patches,
    geometry=[Point(p["lon"], p["lat"]) for p in patches],
    crs="EPSG:4326",
)
scenario["S"] = S_init
scenario["E"] = E_init
scenario["I"] = I_init
scenario["R"] = R_init

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch must have initial infections"

print("Initial conditions:")
for _, row in scenario.iterrows():
    print(f"  {row['name']:<12}: pop={row['population']:>7,}  "
          f"S={row['S']:>6,}  E={row['E']:>5}  I={row['I']:>5,}  R={row['R']:>5,}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. Disease distributions
# ─────────────────────────────────────────────────────────────────────────────

# Latent period: gamma(shape=4, scale=1), mean = 4 days
expdurdist = dists.gamma(shape=4, scale=1.0)

# Infectious period: gamma(shape=5, scale=2), mean = 10 days
infdurdist = dists.gamma(shape=5, scale=2.0)

# beta = R0 / mean_infectious_period = 5 / 10 = 0.5
beta = 5.0 / 10.0


# ─────────────────────────────────────────────────────────────────────────────
# 4. Model parameters
# ─────────────────────────────────────────────────────────────────────────────

# gravity_k/a/b/c in params triggers automatic network computation from
# scenario centroids inside Model.__init__()
params = PropertySet({
    "prng_seed": 42,
    "nticks":    NTICKS,
    "beta":      beta,
    "gravity_k": 0.01,
    "gravity_a": 0.0,   # source exponent (convention: 0 → source pop not weighted)
    "gravity_b": 1.0,   # destination exponent
    "gravity_c": 1.5,   # distance decay exponent
})

# No seasonal forcing — flat profile, mean exactly 1.0
seasonality = ValuesMap.from_scalar(1.0, NTICKS, NNODES)


# ─────────────────────────────────────────────────────────────────────────────
# 5. Build model
# ─────────────────────────────────────────────────────────────────────────────

# birthrates=None → no vital dynamics; capacity = initial agent count
model = Model(scenario, params)

# Validate gravity network
assert model.network.sum() > 0, \
    "Gravity network is all zeros — check gravity_k"
assert model.network.sum(axis=1).max() < 0.3, \
    f"Max row sum {model.network.sum(axis=1).max():.3f} exceeds safe limit"

print(f"\nGravity network row sums (outward coupling fraction per patch):")
for i, p in enumerate(patches):
    print(f"  {p['name']:<12}: {model.network[i].sum():.6f}")


# ─────────────────────────────────────────────────────────────────────────────
# 6. Assemble SEIR components in order
#
# Susceptible and Recovered wrap the transition steps to preserve
# S + E + I + R = N each tick.
#
# Per spec: SEIR.Transmission receives infdurdist as its expdurdist argument.
# ─────────────────────────────────────────────────────────────────────────────

model.components = [
    SEIR.Susceptible(model),                                        # propagates S counts
    SEIR.Transmission(model, infdurdist, seasonality=seasonality),  # S→E; expdurdist=infdurdist per spec
    SEIR.Exposed(model, expdurdist, infdurdist),                    # E→I; etimer←expdurdist, itimer←infdurdist
    SEIR.Infectious(model, infdurdist),                             # I→R; recovery timer←infdurdist
    SEIR.Recovered(model),                                          # propagates R counts
]


# ─────────────────────────────────────────────────────────────────────────────
# 7. Run
# ─────────────────────────────────────────────────────────────────────────────

print(f"\nRunning 1-year Eastland respiratory SEIR "
      f"(beta={beta}, R0=5, {NTICKS} ticks, {NNODES} patches)...")
model.run("Eastland Respiratory SEIR")
print("Simulation complete.")


# ─────────────────────────────────────────────────────────────────────────────
# 8. Post-run verification
# ─────────────────────────────────────────────────────────────────────────────

print("\n--- Verification ---")

# Compartment non-negativity
any_negative = False
for arr, label in [(model.nodes.S, "S"), (model.nodes.E, "E"),
                   (model.nodes.I, "I"), (model.nodes.R, "R")]:
    if np.any(arr[:NTICKS] < 0):
        print(f"  WARNING: Negative {label} values detected!")
        any_negative = True
if not any_negative:
    print("  All compartments non-negative: OK")

# Epidemic dynamics
newly_inf = model.nodes.newly_infected[:NTICKS]
total_inf  = int(newly_inf.sum())
daily_nat  = newly_inf.sum(axis=1)
peak_day   = int(daily_nat.argmax())
peak_val   = int(daily_nat.max())

print(f"  Total infections: {total_inf:,}")
print(f"  Peak daily cases: {peak_val:,} on day {peak_day}")

patches_active = int((newly_inf.sum(axis=0) > 0).sum())
print(f"  Patches with infections: {patches_active}/{NNODES}")

print("\nCumulative infections by patch:")
for i, p in enumerate(patches):
    tot  = int(newly_inf[:, i].sum())
    rate = 100.0 * tot / populations[i]
    print(f"  {p['name']:<12}: {tot:>6,}  ({rate:.1f}% attack rate)")


# ─────────────────────────────────────────────────────────────────────────────
# 9. Diagnostic plot — per-patch SEIR compartment trajectories
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Republic of Eastland — Respiratory SEIR (R0=5, 1 year)",
             fontsize=14, fontweight="bold")
axes = axes.flatten()
days   = np.arange(NTICKS)
colors = {"S": "steelblue", "E": "darkorange", "I": "crimson", "R": "seagreen"}

for i, ax in enumerate(axes):
    ax.plot(days, model.nodes.S[:NTICKS, i], label="S", color=colors["S"], lw=1.5)
    ax.plot(days, model.nodes.E[:NTICKS, i], label="E", color=colors["E"], lw=1.5)
    ax.plot(days, model.nodes.I[:NTICKS, i], label="I", color=colors["I"], lw=1.5)
    ax.plot(days, model.nodes.R[:NTICKS, i], label="R", color=colors["R"], lw=1.5)
    ax.set_title(f"{patches[i]['name']} (pop = {populations[i]:,})")
    ax.set_xlabel("Day")
    ax.set_ylabel("Count")
    ax.legend(loc="center right", fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
outpath = Path("eastland_respiratory_seir.png")
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nPlot saved → {outpath}")
plt.close(fig)
```

---

**Verified output** (ran successfully):

```
Initial conditions:
  Eastport    : pop=100,000  S=90,000  E=    0  I=1,000  R= 9,000
  Northton    : pop=200,000  S=180,000  E=    0  I=2,000  R=18,000
  Centerburg  : pop=150,000  S=135,000  E=    0  I=1,500  R=13,500
  Westham     : pop= 80,000  S=72,000  E=    0  I=  800  R= 7,200

All compartments non-negative: OK
Total infections: 468,833
Peak daily cases: 15,067 on day 33
Patches with infections: 4/4

  Eastport   :  88,113  (88.1% attack rate)
  Northton   : 177,494  (88.7% attack rate)
  Centerburg : 132,856  (88.6% attack rate)
  Westham    :  70,370  (88.0% attack rate)
```

---

**Key design decisions:**

| Decision | Value | Rationale |
|---|---|---|
| `beta` | `0.5` | `R0 / D_I = 5 / 10` |
| `SEIR.Transmission` first arg | `infdurdist` | Per spec: infectious period gamma as `expdurdist` |
| `gravity_k/a/b/c` in `params` | 0.01/0/1/1.5 | Triggers auto-network in `Model.__init__()` |
| `birthrates=None` | — | No vital dynamics needed for 1-year run |
| Component order | S→Tx→Ex→In→R | Susceptible/Recovered wrap transitions; natural SEIR flow |

The ~88% attack rate is epidemiologically consistent with R0=5 and 9% initial immunity: the herd immunity threshold is `1 - 1/R0 = 80%`, so the epidemic overshoots slightly before collapsing.
