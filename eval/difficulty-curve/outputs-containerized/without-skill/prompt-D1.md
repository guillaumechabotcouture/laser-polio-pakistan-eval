The simulation runs cleanly. The working script is at `/tmp/eastland_model.py`. Here's the complete code with an explanation of what it does and the one non-obvious hurdle:

---

```python
"""
Spatial SIR disease model for a respiratory illness in the Republic of Eastland.
Uses the LASER framework (laser-generic 0.1.0 + laser-core 0.4.0).

laser-generic 0.1.0 targets laser-core >=0.8.1; the installed version is
0.4.0.  The three missing laser.core sub-modules are bridged below via
sys.modules patches before any laser.generic imports occur.
"""

import sys
import types

import numpy as np
import numba as nb
import geopandas as gpd
from shapely.geometry import Point

# ── Patch missing laser.core sub-modules ────────────────────────────────────

import laser_core.migration
import laser_core.random
import laser.core as _lc

sys.modules["laser.core.migration"] = laser_core.migration
sys.modules["laser.core.random"] = laser_core.random

# laser.core.utils.calc_capacity in >=0.8.1 takes
# (birthrates, population, *, safety_factor) and returns a per-node array.
_utils = types.ModuleType("laser.core.utils")

def _calc_capacity(birthrates, population, *, safety_factor=1.0):
    pop = np.asarray(population, dtype=np.float64)
    # dtype=object preserves Python int so that .sum() returns a Python int,
    # satisfying LaserFrame's isinstance(capacity, int) guard.
    values = [int(np.ceil(p * safety_factor)) for p in pop]
    return np.array(values, dtype=object)

_utils.calc_capacity = _calc_capacity
sys.modules["laser.core.utils"] = _utils

# laser.core.distributions.sample_floats(dist, array) fills an output
# array in-place by drawing len(array) samples from the distribution callable.
_dists = types.ModuleType("laser.core.distributions")

def _sample_floats(distribution, output_array):
    n = len(output_array)
    # distribution may be an njit (tick, nid) -> float function or a
    # plain (n,) -> array callable.  Try vectorised call first; on TypeError
    # (e.g. wrong arity inside njit), fall back to per-element sampling.
    try:
        samples = distribution(n)
        output_array[:] = np.asarray(samples, dtype=output_array.dtype)
    except TypeError:
        for i in range(n):
            output_array[i] = distribution(0, i)
    return output_array

_dists.sample_floats = _sample_floats
sys.modules["laser.core.distributions"] = _dists
_lc.distributions = _dists

# ── Import LASER classes (after patching) ────────────────────────────────────
from laser_core import PropertySet
from laser.generic.model import Model
from laser.generic.SIR import Susceptible, Transmission, Infectious, Recovered

# ── 4-patch GeoDataFrame for the Republic of Eastland ───────────────────────
patches = gpd.GeoDataFrame(
    {
        "name":       ["Northern District", "Capital Region",
                       "Eastern Province",  "Southern Coast"],
        "nodeid":     [0, 1, 2, 3],
        "population": [100_000, 200_000, 150_000, 80_000],
        "S":          [100_000, 200_000, 150_000, 80_000],  # all susceptible
        "E":          [0, 0, 0, 0],
        "I":          [0, 0, 0, 0],                         # no initial infections
        "R":          [0, 0, 0, 0],
        "geometry": [
            Point(30.0, 52.0),   # Northern District
            Point(31.5, 51.0),   # Capital Region
            Point(32.5, 51.5),   # Eastern Province
            Point(31.0, 49.5),   # Southern Coast
        ],
    },
    crs="EPSG:4326",
)

# ── Simulation parameters ────────────────────────────────────────────────────
params = PropertySet({
    "nticks":    365,    # 1-year simulation
    "prng_seed": 42,
    "beta":      0.30,   # per-day transmission rate (~R0 2.1 with 7-day mean inf. period)
    "gravity_k": 1.0,    # gravity migration model parameters (IDM defaults)
    "gravity_a": 1.0,
    "gravity_b": 1.0,
    "gravity_c": 2.0,
})

# ── Infectious-duration distribution ─────────────────────────────────────────
# TransmissionSI calls infdurdist inside a @nb.njit(parallel=True) kernel
# as infdurdist(tick, nid) -> float32, so it must itself be @nb.njit.

@nb.njit
def infdurdist(tick: int, nid: int) -> nb.float32:
    """Return an infectious-period duration (ticks) sampled from Exp(mean=7)."""
    return np.float32(np.random.exponential(7.0))

# ── Instantiate the Model ─────────────────────────────────────────────────────
model = Model(
    scenario=patches,
    params=params,
    birthrates=None,              # closed population – no vital dynamics
    name="Eastland-Respiratory-SIR",
    skip_capacity=False,
    states={"S", "I", "R"},      # SIR compartments (no exposed class)
)

# ── Attach SIR disease components ─────────────────────────────────────────────
model.components = [
    Susceptible(model),
    Infectious(model, infdurdist),
    Transmission(model, infdurdist),
    Recovered(model),
]

# ── Run the 1-year simulation ─────────────────────────────────────────────────
print("Republic of Eastland - Respiratory Illness SIR Model")
print(f"  Patches : {len(patches)}")
print(f"  Agents  : {model.people.count:,}")
print(f"  Ticks   : {params.nticks}")
print()
print("Running simulation...")
model.run()
print("Simulation complete.")
```

---

**Output:**
```
Republic of Eastland - Respiratory Illness SIR Model
  Patches : 4
  Agents  : 530,000
  Ticks   : 365

Running simulation...
530,000 agents in 4 node(s): 100%|██████████| 365/365
Simulation complete.
```

**Key design decisions:**

- **Version bridge** — `laser-generic 0.1.0` requires `laser-core >=0.8.1` but only `0.4.0` is installed. Three `sys.modules` patches register the available `laser_core.*` submodules under the `laser.core.*` namespace, and a minimal `laser.core.distributions` stub provides `sample_floats`. All patches are applied before any `laser.generic` import.

- **`PropertySet` takes a dict** — not keyword arguments: `PropertySet({"nticks": 365, ...})`.

- **`infdurdist` must be `@nb.njit`** — `TransmissionSI.step()` calls `infdurdist(tick, nid)` inside a `@nb.njit(parallel=True)` kernel. A plain Python function would cause a Numba typing error; the function must be JIT-compiled with signature `(int, int) -> float32`.

- **`_calc_capacity` returns `dtype=object` array** — `model.py` calls `.sum()` on the result and passes it directly to `LaserFrame(capacity=...)`, which requires a Python `int`. A `dtype=object` array whose elements are Python ints yields a Python `int` from `.sum()`.

- **`states={"S", "I", "R"}`** — omits the `"E"` (Exposed) compartment, appropriate for a simple SIR model with no latent period.
