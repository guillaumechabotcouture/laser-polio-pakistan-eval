Here is the complete code with explanation:

```python
import numpy as np
from laser.core.utils import calc_capacity

# ---------------------------------------------------------------------------
# Scenario parameters
# ---------------------------------------------------------------------------
initial_pop = 500_000      # persons
cbr         = 44.0         # per-1,000 per year  <-- units calc_capacity expects
sim_years   = 10
nticks      = sim_years * 365   # 3,650 daily ticks
nnodes      = 1                 # single-patch

# ---------------------------------------------------------------------------
# Build the birthrates array
# Shape: (nticks, nnodes) — calc_capacity infers nticks from shape[0]
# Values: per-1,000/year.  A value of 44 means 44 births per 1,000 per year.
#
# WRONG: 44 / 1000 / 365 = 0.000121  (daily per-capita) → frame severely
#         under-allocated → LaserFrame.add() silently returns no slots →
#         no births occur for the entire simulation.
# CORRECT: 44.0  (the framework divides by 1000 and converts to daily internally)
# ---------------------------------------------------------------------------
birthrates = np.full((nticks, nnodes), cbr, dtype=np.float32)

# Unit guard — catches the most common silent failure
assert np.all(birthrates >= 1) and np.all(birthrates <= 100), (
    f"birthrates must be per-1,000/year. "
    f"Got min={birthrates.min():.6f}, max={birthrates.max():.6f}"
)

# ---------------------------------------------------------------------------
# Per-node initial population — shape (nnodes,)
# ---------------------------------------------------------------------------
initial_pop_arr = np.array([initial_pop], dtype=np.int32)

# ---------------------------------------------------------------------------
# calc_capacity
#
# Signature:
#   capacity = calc_capacity(birthrates, initial_pop, safety_factor=1.0)
#
# What it does:
#   Projects cumulative births over nticks using the daily compound formula
#     births_per_day = N × ((1 + CBR/1000)^(1/365) - 1)
#   Sums those births across the simulation, adds the initial population,
#   then multiplies by safety_factor to give headroom for stochastic variation.
#   Returns np.int32 array of shape (nnodes,).
#
# safety_factor guidance:
#   1.0  — bare minimum; any above-average birth run overflows
#   2–3  — recommended for CBR > 30 or multi-year runs
#   4    — conservative choice for > 20-year simulations
# ---------------------------------------------------------------------------
safety_factor = 2.5
capacity = calc_capacity(birthrates, initial_pop_arr, safety_factor=safety_factor)

print(f"Initial population   : {initial_pop:>10,}")
print(f"CBR                  : {cbr:.1f} per 1,000/year")
print(f"Simulation length    : {sim_years} years ({nticks:,} daily ticks)")
print(f"Safety factor        : {safety_factor}x")
print(f"Allocated capacity   : {capacity[0]:>10,}")
print(f"Headroom             : {capacity[0] / initial_pop:.2f}× initial population")

# ---------------------------------------------------------------------------
# How this connects to Model construction
#
# Model.__init__() calls calc_capacity() internally when you pass birthrates.
# You do NOT need to call it yourself in normal usage — it is shown here to
# make the allocation logic explicit.
#
# Equivalent Model usage (single-patch, minimal):
# ---------------------------------------------------------------------------
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.generic import Model

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     [0],
        "name":       ["node_0"],
        "population": [initial_pop],
        "S":          [initial_pop - 10],
        "E":          [0],
        "I":          [10],
        "R":          [0],
        "geometry":   [Point(74.0, 30.0)],   # lon, lat — Pakistan centroid
    },
    crs="EPSG:4326",
)

params = PropertySet({
    "prng_seed": 42,
    "nticks":    nticks,
    "beta":      14.0,    # R0 ~ 14 for polio-like transmission
    "cbr":       cbr,
    "capacity_safety_factor": safety_factor,
})

# Model.__init__() uses birthrates.shape[0] as nticks and calls calc_capacity
# with the same safety_factor from params.  The LaserFrame for model.people
# is pre-allocated to capacity[0] slots before any agents are created.
model = Model(scenario, params, birthrates=birthrates)

print(f"\nmodel.people capacity : {model.people.capacity:>10,}")
print(f"model.people count    : {model.people.count:>10,}")
print(f"Matches calc_capacity : {model.people.capacity == capacity[0]}")
```

---

## What units `calc_capacity` expects

| Parameter | Expected units | Example value | Common wrong value |
|---|---|---|---|
| `birthrates` | **per-1,000 per year** | `44.0` | `0.000121` (daily per-capita) |
| `initial_pop` | raw agent count | `500_000` | — |
| `safety_factor` | dimensionless multiplier | `2.5` | `1.0` (too tight) |

The framework converts the per-1000/year CBR to a daily compound rate internally using:

```
births_per_day = N × ((1 + CBR/1000)^(1/365) - 1)
```

So `44.0` is the correct value to pass — not `44/1000` and not `44/1000/365`.

**Why the safety factor matters here:** With CBR = 44 per 1,000/year over 10 years and no deaths, an initial population of 500,000 grows to roughly 770,000 (≈54% increase). A safety factor of 2.5 allocates ~1.9M slots, giving ample headroom for stochastic birth clustering. If `LaserFrame.add()` runs out of capacity during the run it raises a `ValueError` — so over-allocating is far preferable to under-allocating.
