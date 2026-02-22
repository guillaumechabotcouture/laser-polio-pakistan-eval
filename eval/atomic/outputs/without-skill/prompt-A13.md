Here's the complete code and explanation:

```python
import numpy as np
from laser.core.utils import calc_capacity

# --- Simulation parameters ---
initial_pop   = 500_000       # agents at t=0
cbr            = 44.0         # crude birth rate: 44 per 1,000 per year
years          = 10
ticks_per_year = 365
nticks         = years * ticks_per_year   # 3,650 daily time steps

# --- Build the inputs calc_capacity expects ---

# birthrates: 2D array of shape (nticks, nnodes), units = CBR per 1,000/year
# We have a single node, so nnodes=1
birthrates = np.full((nticks, 1), cbr, dtype=np.float32)

# initial_pop: 1D array of shape (nnodes,), units = number of agents
population = np.array([initial_pop], dtype=np.int32)

# --- Calculate capacity ---
capacity = calc_capacity(birthrates, population, safety_factor=1.0)

print(f"Pre-allocated agent capacity: {capacity[0]:,}")
```

## What units does `calc_capacity` expect?

| Parameter | Shape | Units |
|---|---|---|
| `birthrates` | `(nticks, nnodes)` — 2D | **Crude birth rate per 1,000 individuals per year**. The value 44 means 44 births per 1,000 people per year. The array is 2D so you can vary the rate across both time steps and spatial nodes. |
| `initial_pop` | `(nnodes,)` — 1D | **Absolute count of agents** (people) per node at t=0. |
| `safety_factor` | scalar | **Dimensionless** (0–6). Default `1.0`. Adds a stochastic buffer above expected growth. |

**Return value:** `np.ndarray` of shape `(nnodes,)`, dtype `int32` — the number of agent slots to pre-allocate per node.

## How the calculation works internally

1. **Daily growth rate** — Converts the annual CBR to a daily compound rate:
   ```
   λ = (1 + CBR/1000)^(1/365) - 1
   ```

2. **Expected growth multiplier** — Sums the daily rates over all ticks and exponentiates (a Geometric Brownian Motion approximation):
   ```
   exp_μt = exp(Σ λ)
   ```

3. **Safety buffer** — Accounts for stochastic variability in births:
   ```
   safety_multiplier = 1 + safety_factor × (√exp_μt - 1)
   capacity = round(initial_pop × safety_multiplier × exp_μt)
   ```

With CBR=44 over 10 years, the population roughly grows by a factor of ~1.55, so you'd expect a capacity around 775,000–800,000 agents (exact value depends on the safety factor). This ensures the `LaserFrame` agent array never runs out of slots during the simulation.
