Now I have the full picture. Here's the complete answer:

---

## `MortalityByCDR` — Units and Usage

### Units: **per 1000 population per year**

`MortalityByCDR` expects its `mortalityrates` parameter as **crude death rates in units of deaths per 1,000 population per year**. Internally (see `vitaldynamics.py:186-190`), it converts this annual rate to a daily mortality probability each tick:

```
annual_survival = 1 - CDR/1000
daily_survival  = annual_survival^(1/365)
daily_rate      = 1 - daily_survival
daily_prob      = 1 - exp(-daily_rate)     # = -expm1(-daily_rate)
```

So you pass `13` (not `0.013`, not a daily value) for 13 deaths per 1,000 per year.

### Complete Example

```python
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.generic.model import Model
from laser.generic import SIR
from laser.generic.vitaldynamics import MortalityByCDR

# --- 1. Define simulation parameters ---
nticks = 365  # 1 year, each tick = 1 day
params = PropertySet(
    nticks=nticks,
    beta=0.3,       # transmission rate
)

# --- 2. Build a scenario (one node, 100k people) ---
population = 100_000
scenario = gpd.GeoDataFrame({
    "nodeid":     [0],
    "population": [population],
    "S":          [population - 10],   # susceptible at t=0
    "I":          [10],                # infectious at t=0
    "R":          [0],
    "E":          [0],
    "geometry":   [Point(0, 0)],
})

# --- 3. Set up mortality rates ---
# MortalityByCDR expects: deaths per 1000 population per YEAR.
# A CDR of 13 means 13 deaths per 1,000 people per year.
# Pass one value per tick (per node). Here: constant across all ticks.
cdr = 13.0  # per 1000 per year
mortalityrates = np.full((nticks, 1), cdr, dtype=np.float64)

# --- 4. Create the model ---
model = Model(scenario, params)

# --- 5. Assemble components (order matters) ---
model.components = [
    SIR.Susceptible(model),
    SIR.Infectious(model, infdurdist=lambda tick, nid: np.random.exponential(7)),
    SIR.Recovered(model),
    SIR.Transmission(model, infdurdist=lambda tick, nid: np.random.exponential(7)),
    MortalityByCDR(model, mortalityrates),  # <-- mortality process
]

# --- 6. Run ---
model.run()

# --- 7. Inspect results ---
total_deaths = model.nodes.deaths.sum()
print(f"Total deaths over {nticks} days: {total_deaths}")
# Expected: ~100,000 * (13/1000) * (365/365) ≈ 1,300 deaths/year
```

### Key points about `MortalityByCDR`

| Parameter | Type | Description |
|---|---|---|
| `model` | `Model` | The simulation model instance |
| `mortalityrates` | `np.ndarray` shape `(nticks, nnodes)` | **CDR in deaths per 1,000 per year** for each tick and node |
| `mappings` | `list[(int, str)]` | State-to-name tuples for decrementing counts. Defaults to S, E, I, R |
| `validating` | `bool` | Enable validation checks (default `False`) |

- Each tick, agents are probabilistically marked `DECEASED` based on the converted daily probability.
- Deaths are subtracted from the appropriate compartment counts (S, E, I, R) via the `mappings`.
- The component records per-node deaths each tick in `model.nodes.deaths`.
