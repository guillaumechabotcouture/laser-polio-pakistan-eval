Now I have everything needed from the actual source code. Here is the complete example:

---

**Units:** `MortalityByCDR` expects `mortalityrates[tick]` to be an array of **crude death rates in deaths per 1,000 population per year**, one value per spatial node. It converts this internally to a per-tick (daily) death probability:

```
annual_survival = 1 - CDR / 1000
daily_survival  = annual_survival ** (1 / 365)
daily_p_death   = -expm1(-(1 - daily_survival))
```

So pass `13.0` to represent 13 deaths per 1,000 per year — **not** a rate per day, per tick, or per capita directly.

---

```python
"""
LASER SEIR model with MortalityByCDR.

MortalityByCDR unit contract
-----------------------------
mortalityrates[tick] → array of shape (num_nodes,)
  Values are **crude death rates in deaths per 1,000 population per year**.
  The component converts to a daily probability each tick via:

      annual_survival = 1 - CDR / 1000
      daily_survival  = annual_survival ** (1 / 365)
      daily_p_death   = -expm1(-(1 - daily_survival))

  For CDR = 13.0:
      annual_survival ≈ 0.987
      daily_survival  ≈ 0.999964
      daily_p_death   ≈ 3.57 × 10⁻⁵
"""

import numpy as np
import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.utils import grid
from laser.generic import SEIR, Model, State
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import MortalityByCDR

# ── Simulation settings ───────────────────────────────────────────────────────

NTICKS    = 3 * 365   # 3 years of daily ticks
CDR       = 13.0      # crude death rate: 13 deaths per 1,000 per year
INIT_POP  = 100_000
SEED      = 42

# ── Scenario: single spatial node ────────────────────────────────────────────
# grid() returns a GeoDataFrame with columns: nodeid, population, geometry.

scenario = grid(
    M=1, N=1,
    node_size_degs=0.08983,               # ~10 km cells
    population_fn=lambda r, c: INIT_POP,
    origin_x=-119.204167,                 # Black Rock Desert, NV (lon)
    origin_y=40.786944,                   # Black Rock Desert, NV (lat)
)

# Distribute initial population across SEIR compartments.
scenario["E"] = (scenario.population * 0.05).astype(np.int32)
scenario["I"] = (scenario.population * 0.05).astype(np.int32)
scenario["R"] = (scenario.population * 0.20).astype(np.int32)
scenario["S"] = (scenario.population - scenario.E - scenario.I - scenario.R).astype(np.int32)

# ── Parameters ────────────────────────────────────────────────────────────────

params = PropertySet({
    "nticks": NTICKS,
    "seed":   SEED,
})

# ── Mortality rates ───────────────────────────────────────────────────────────
# ValuesMap.from_scalar(value, nticks, num_nodes) creates a constant-rate array
# where mortalityrates[tick] returns a shape-(num_nodes,) array every tick.
# The value MUST be in **deaths per 1,000 per year** — MortalityByCDR converts
# it to a daily probability internally (see module docstring above).

NUM_NODES      = len(scenario)
mortalityrates = ValuesMap.from_scalar(CDR, NTICKS, NUM_NODES)

# ── Disease duration distributions ───────────────────────────────────────────

expdurdist = dists.normal(loc=5.0,  scale=1.0)   # latent (exposure) period, days
infdurdist = dists.normal(loc=10.0, scale=2.0)   # infectious period, days

# ── Build model and attach components ────────────────────────────────────────

model = Model(scenario, params, birthrates=None)

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    # Pass mortalityrates in **deaths per 1,000 per year**.
    MortalityByCDR(model, mortalityrates),
]

# ── Run ───────────────────────────────────────────────────────────────────────

print(f"Starting population : {INIT_POP:,}")
print(f"CDR                 : {CDR} per 1,000 per year")
print(f"Simulation duration : {NTICKS} days ({NTICKS / 365:.1f} years)\n")

model.run()

# ── Results ───────────────────────────────────────────────────────────────────

alive_end    = int((model.people.state != State.DECEASED.value).sum())
total_deaths = int(model.nodes.deaths.sum())

print(f"\nFinal living population : {alive_end:,}")
print(f"Total deaths recorded   : {total_deaths:,}")

# Compute observed CDR for each calendar year of the run.
# N[t] = total active (non-deceased) agents at tick t, summed across compartments.
N = model.nodes.S + model.nodes.E + model.nodes.I + model.nodes.R  # shape (nticks+1, num_nodes)

print("\nObserved CDR by year (per 1,000 per year):")
for year in range(NTICKS // 365):
    t_start = year * 365
    t_end   = t_start + 365
    n_start = int(N[t_start, 0])
    n_end   = int(N[t_end,   0])
    obs_cdr = 1000 * (n_start - n_end) / n_start
    print(f"  Year {year + 1}: {obs_cdr:.2f}  (target: {CDR:.1f})")
```

---

**Key assembly notes:**

| Item | Detail |
|------|--------|
| `mortalityrates` type | `ValuesMap` (supports `[tick]` indexing), one entry per tick returning shape `(num_nodes,)` |
| Input units | **Deaths per 1,000 per year** — the only unit `MortalityByCDR` accepts |
| CDR = 13 → daily prob | ≈ 3.57 × 10⁻⁵ per agent per day (converted internally) |
| `model.nodes.deaths` | Populated by `MortalityByCDR`; shape `(nticks+1, num_nodes)`, raw counts per tick |
| Component order | Mortality comes after the SEIR transition components so deaths reduce the correct compartment counts for tick `t+1` |
