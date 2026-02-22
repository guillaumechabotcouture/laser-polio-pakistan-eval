The filesystem is read-only, so I'll provide the complete code directly:

```python
#!/usr/bin/env python3
"""
Republic of Eastland — Basic 4-Patch Spatial SEIR Model for Respiratory Illness

Specification:
  - 4 patches: Northburg (100k), Eastholm (200k), Southfield (150k), Westgate (80k)
  - All susceptible at t=0 (S = population, E = I = R = 0)
  - PropertySet: nticks=365, prng_seed=42
  - Respiratory illness: R0 ≈ 2.5, latent ~3 days, infectious ~5 days
  - Gravity network auto-computed from patch centroids
  - 10 seeded infections/day in Eastholm (days 0–7) to initiate the epidemic
  - 1-year simulation via model.run()
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.importation import Infect_Agents_In_Patch


# ── 1. Geographic Scenario ──────────────────────────────────────────────────
patches = [
    {"nodeid": 0, "name": "Northburg",  "population": 100_000, "lat": 51.0, "lon": 17.0},
    {"nodeid": 1, "name": "Eastholm",   "population": 200_000, "lat": 51.5, "lon": 18.0},
    {"nodeid": 2, "name": "Southfield", "population": 150_000, "lat": 50.5, "lon": 18.5},
    {"nodeid": 3, "name": "Westgate",   "population":  80_000, "lat": 50.0, "lon": 17.5},
]

df = pd.DataFrame(patches)
df["geometry"] = [Point(row["lon"], row["lat"]) for _, row in df.iterrows()]

# All susceptible: S = population, E = I = R = 0
n = len(df)
df["S"] = df["population"].astype(np.uint32)
df["E"] = np.zeros(n, dtype=np.uint32)
df["I"] = np.zeros(n, dtype=np.uint32)
df["R"] = np.zeros(n, dtype=np.uint32)

scenario = gpd.GeoDataFrame(df, crs="EPSG:4326")

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"


# ── 2. Parameters ───────────────────────────────────────────────────────────
# Respiratory illness (influenza-like):
#   R0 ≈ 2.5, beta = R0 / inf_mean = 2.5 / 5 = 0.5
#   Latent period:     gamma(shape=3, scale=1) → mean = 3 days
#   Infectious period: normal(mean=5, std=1)   → mean = 5 days
#
# Gravity network (auto-computed by Model.__init__ when gravity_k/a/b/c
# are present in params):
#   M_{i,j} = k * p_j^b / d_{ij}^c  (a=0: source-population-independent)
#   k=0.001 keeps row sums well below the 0.3 safety cap for patches
#   ~60–180 km apart.
#
# Importation: Infect_Agents_In_Patch reads these keys from model.params.
# Seeding in the scenario itself is left at zero to satisfy "all susceptible".

params = PropertySet({
    "prng_seed":  42,
    "nticks":     365,
    # Transmission
    "beta":       0.5,
    # Exposed duration
    "exp_shape":  3.0,
    "exp_scale":  1.0,
    # Infectious duration
    "inf_mean":   5.0,
    "inf_sigma":  1.0,
    # Gravity network
    "gravity_k":  0.001,
    "gravity_a":  0.0,   # LASER convention: outward flow is source-pop-independent
    "gravity_b":  1.0,
    "gravity_c":  2.0,
    # Importation: 10 infections/day in Eastholm (patch 1) for the first week
    "importation_period":    1,
    "importation_count":     10,
    "importation_start":     0,
    "importation_end":       7,
    "importation_patchlist": [1],
})

nticks = params.nticks
nnodes = len(scenario)


# ── 3. Duration Distributions ───────────────────────────────────────────────
expdurdist = dists.gamma(shape=params.exp_shape, scale=params.exp_scale)
infdurdist = dists.normal(loc=params.inf_mean,   scale=params.inf_sigma)


# ── 4. Seasonality — no forcing for this basic model ───────────────────────
# Constant 1.0 multiplier: effective beta = params.beta at every tick.
seasonality = ValuesMap.from_scalar(1.0, nticks, nnodes)


# ── 5. Build Model ──────────────────────────────────────────────────────────
# Model.__init__ automatically:
#   - Seeds the PRNG (params.prng_seed = 42)
#   - Allocates a LaserFrame for all agents (capacity = total population;
#     no births so no growth buffer is needed)
#   - Computes the gravity network from scenario centroids because
#     gravity_k/a/b/c are present in params

model = Model(scenario, params)


# ── 6. Validate Gravity Network ─────────────────────────────────────────────
assert model.network.sum() > 0, \
    "Gravity network is all zeros — check gravity_k and patch coordinates"
assert model.network.sum(axis=1).max() < 0.3, (
    f"Network row sum {model.network.sum(axis=1).max():.4f} exceeds safety "
    "cap — reduce gravity_k or increase distance decay (gravity_c)"
)


# ── 7. Assemble Components ──────────────────────────────────────────────────
# Execution order each tick:
#   Susceptible  — propagates S[t] → S[t+1]
#   Exposed      — decrements etimer; E→I on expiry
#   Infectious   — decrements itimer; I→R on expiry
#   Recovered    — propagates R[t] → R[t+1]
#   Importation  — seeds 10 infections/day into patch 1 during days 0–7
#   Transmission — computes FOI with spatial coupling, transitions S→E

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    Infect_Agents_In_Patch(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
]


# ── 8. Run ──────────────────────────────────────────────────────────────────
total_pop = int(scenario.population.sum())

print("Republic of Eastland — 1-Year Respiratory SEIR Simulation")
print(f"Patches:          {list(scenario['name'])}")
print(f"Total population: {total_pop:,} (all susceptible at t=0)")
print(f"R0 ≈ {params.beta * params.inf_mean:.1f}  |  "
      f"beta={params.beta}  |  latent~{params.exp_shape * params.exp_scale:.0f}d  |  "
      f"infectious~{params.inf_mean:.0f}d")
print(f"Importation:      {params.importation_count} infections/day into Eastholm "
      f"(days {params.importation_start}–{params.importation_end})")
print(f"Network max row sum: {model.network.sum(axis=1).max():.4f}")
print()

model.run("Eastland Respiratory Simulation")


# ── 9. Summary Output ───────────────────────────────────────────────────────
print("\n=== Simulation Complete ===\n")
print(f"{'Patch':<12} {'S_final':>9} {'E_final':>8} {'I_final':>8} "
      f"{'R_final':>9} {'Attack%':>8}")
print("─" * 58)

total_ever_infected = 0
for i, row in scenario.iterrows():
    s_f  = int(model.nodes.S[nticks - 1, i])
    e_f  = int(model.nodes.E[nticks - 1, i])
    i_f  = int(model.nodes.I[nticks - 1, i])
    r_f  = int(model.nodes.R[nticks - 1, i])
    pop  = int(row["population"])
    ever = int(model.nodes.newly_infected[:nticks, i].sum())
    total_ever_infected += ever
    print(f"{row['name']:<12} {s_f:>9,} {e_f:>8,} {i_f:>8,} "
          f"{r_f:>9,} {ever / pop * 100:>7.1f}%")

print("─" * 58)
print(f"{'National':<12} attack rate: {total_ever_infected / total_pop * 100:.1f}%  |  "
      f"infected: {total_ever_infected:,} / {total_pop:,}")

# Compartment non-negativity check
for arr, label in [(model.nodes.S, "S"), (model.nodes.E, "E"),
                   (model.nodes.I, "I"), (model.nodes.R, "R")]:
    if np.any(arr[:nticks] < 0):
        print(f"WARNING: Negative {label} values detected — check component logic")
```

---

**Key design decisions:**

**All-susceptible with importation.** Setting `S = population, E = I = R = 0` in the scenario honors the spec. Since the epidemic can't ignite from I=0, `Infect_Agents_In_Patch` seeds 10 infections/day into Eastholm (the largest patch) for days 0–7 via `model.params` — without pre-seeding any I counts in the scenario GeoDataFrame.

**Disease parameters (respiratory, R0 ≈ 2.5).**
- `beta = 0.5`: LASER computes daily force of infection as `beta × I/N`, so `R0 ≈ beta × inf_mean = 0.5 × 5 = 2.5`
- Latent period: `gamma(shape=3, scale=1)` → mean 3 days
- Infectious period: `normal(mean=5, std=1)` → mean 5 days

**Gravity network auto-setup.** Including `gravity_k/a/b/c` in `params` triggers `Model.__init__` to compute Haversine distances from the scenario centroids and build `model.network` automatically. `gravity_k=0.001` keeps row sums ~0.05–0.08 for patches 60–180 km apart, well inside the < 0.3 safety cap (validated by assertion after `Model()`).

**Component ordering** follows the reference pattern: `Susceptible → Exposed → Infectious → Recovered → Importation → Transmission`. No vital dynamics are needed for a 1-year, closed-population simulation.
