The complete code is at `/home/evaluser/eastland_seir.py`. Here it is with a walkthrough of every design decision:

---

```python
"""
Spatial SEIR model for the Republic of Eastland using the LASER framework.
...
"""

import numba as nb
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core import PropertySet
from laser.generic import Model
from laser.generic.SEIR import Exposed, Infectious, Recovered, Susceptible, Transmission

# ── Reproducibility ─────────────────────────────────────────────────────────
PRNG_SEED = 42

# ── Disease parameters ───────────────────────────────────────────────────────
R0 = 5.0
INF_SHAPE, INF_SCALE = 5.0, 2.0   # Gamma → mean = 10 days
EXP_SHAPE, EXP_SCALE = 4.0, 1.0   # Gamma → mean =  4 days
MEAN_INFECTIOUS = INF_SHAPE * INF_SCALE   # 10 days
beta = R0 / MEAN_INFECTIOUS               # 0.5 / day

# ── Duration distribution callables ─────────────────────────────────────────
# Must be @nb.njit: called inside numba-parallel kernels (nb_timer_update_timer_set,
# nb_transmission_step) AND from Python (laser.core.distributions.sample_floats).
# Signature: f(tick: int, node_id: int) -> float

@nb.njit
def expdurdist(tick, node_id):
    """Latent period: Gamma(4, 1), mean = 4 days."""
    return np.random.gamma(EXP_SHAPE, EXP_SCALE)

@nb.njit
def infdurdist(tick, node_id):
    """Infectious period: Gamma(5, 2), mean = 10 days."""
    return np.random.gamma(INF_SHAPE, INF_SCALE)

# ── Scenario: 4 patches ──────────────────────────────────────────────────────
populations = [100_000, 200_000, 150_000, 80_000]
lons = [28.0, 29.0, 28.0, 29.0]
lats = [12.0, 12.0, 13.0, 13.0]

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     [0, 1, 2, 3],
        "population": populations,
        "S": [int(p * 0.90) for p in populations],
        "E": [int(p * 0.00) for p in populations],
        "I": [int(p * 0.01) for p in populations],
        "R": [int(p * 0.09) for p in populations],
    },
    geometry=[Point(lon, lat) for lon, lat in zip(lons, lats)],
    crs="EPSG:4326",
)

# ── Parameters ───────────────────────────────────────────────────────────────
params = PropertySet({"nticks": 365, "beta": beta, "seed": PRNG_SEED})

# ── Model & components ───────────────────────────────────────────────────────
model = Model(scenario, params, name="eastland_seir")

model.components = [
    Susceptible(model),                                    # 1. marks all → S; S[0]
    Exposed(model, expdurdist=expdurdist,
                   infdurdist=infdurdist),                 # 2. seeds E; E→I via etimer
    Infectious(model, infdurdist=infdurdist),              # 3. seeds I; I→R via itimer
    Recovered(model),                                      # 4. seeds R; passive tracker
    Transmission(model, expdurdist=infdurdist),            # 5. S→E per tick (per spec)
]

model.run()
```

---

**Key design decisions:**

| Choice | Rationale |
|---|---|
| `beta = R0 / mean_infectious = 0.5` | Standard SEIR relationship β = R₀ / T_I |
| `@nb.njit` on both distributions | LASER's step kernels (`nb_timer_update_timer_set`, `nb_transmission_step`) are `@nb.njit(parallel=True)` and invoke the distribution as a first-class function — it must be JIT-compiled |
| `np.random.gamma` (legacy API) | The only numpy random API supported in numba nopython mode |
| `EXP_SHAPE`/`EXP_SCALE` as module-level constants | numba can lift these as compile-time constants |
| Component order: S, E, I, R, Tx | Each compartment's `__init__` seeds agents from the remaining susceptible pool; Transmission only runs at step-time, so it goes last |
| `Transmission(model, expdurdist=infdurdist)` | Per specification: uses the Gamma(5, 2) infectious-period distribution as `expdurdist`; this sets `etimer` for newly exposed agents from that distribution |
| `nodeid` as a DataFrame column | `Model.__init__` reads `scenario.nodeid` as a column to compute `num_nodes = max(unique(nodeid)) + 1` |
| `crs="EPSG:4326"` with lon/lat | `Model` reprojects to EPSG:3857 internally to compute gravity-model distances between patch centroids |

**Note:** `laser.generic 1.0.0` declares `Requires: laser.core >= 1.0.0`. That version of `laser.core` (which provides `laser.core.distributions.sample_floats` used in component `__init__` methods) is not yet on PyPI; install it from the IDM development channel or a private index alongside `laser.generic`.
