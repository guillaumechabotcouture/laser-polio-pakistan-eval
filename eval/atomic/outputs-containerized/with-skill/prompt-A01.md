Here is the complete code:

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.generic import Model

# ── 1. Four-patch GeoDataFrame ───────────────────────────────────────────────
# Arrange patches on a simple 2×2 degree grid so distance calculations
# (used internally for gravity network setup) are well-defined.
patch_data = {
    "nodeid":     [0, 1, 2, 3],
    "name":       ["Patch_0", "Patch_1", "Patch_2", "Patch_3"],
    "population": [50_000, 50_000, 50_000, 50_000],
    "geometry": [
        Point(0.0, 0.0),   # lon, lat
        Point(1.0, 0.0),
        Point(0.0, 1.0),
        Point(1.0, 1.0),
    ],
    # All agents start susceptible
    "S": [50_000, 50_000, 50_000, 50_000],
    "E": [     0,      0,      0,      0],
    "I": [     0,      0,      0,      0],
    "R": [     0,      0,      0,      0],
}

scenario = gpd.GeoDataFrame(patch_data, crs="EPSG:4326")

# Sanity check: compartments must sum to population in every patch
assert (
    scenario.S + scenario.E + scenario.I + scenario.R == scenario.population
).all(), "S+E+I+R must equal population in every patch"

# ── 2. Parameter set ─────────────────────────────────────────────────────────
# Model.__init__ searches params for prng_seed → prngseed → seed (in that order).
# nticks and beta are required at minimum by the Model constructor.
params = PropertySet({
    "nticks":     365,
    "prng_seed":  20250101,
    "beta":       0.0,   # placeholder — set >0 and add a Transmission component to run
})

# ── 3. Instantiate the Model ──────────────────────────────────────────────────
# birthrates=None: no vital dynamics; Model skips capacity growth projection.
# No gravity_k/a/b/c in params, so model.network is initialised to zeros.
model = Model(scenario, params)

# ── 4. Verify the instance ────────────────────────────────────────────────────
nnodes = len(scenario)

print(f"Patches:          {nnodes}")
print(f"nticks:           {model.params.nticks}")
print(f"prng_seed:        {model.params.prng_seed}")
print(f"People capacity:  {model.people.capacity:,}")
print(f"Node S[tick=0]:   {model.nodes.S[0].tolist()}")
print(f"Network shape:    {model.network.shape}")
print(f"Network nonzero:  {np.count_nonzero(model.network)}")  # 0 — no gravity params

assert model.nodes.S[0].sum() == 200_000, "Total initial S must equal 4 × 50,000"
assert model.network.shape == (nnodes, nnodes), "Network must be (4, 4)"
```

**Expected output:**
```
Patches:          4
nticks:           365
prng_seed:        20250101
People capacity:  200000
Node S[tick=0]:   [50000, 50000, 50000, 50000]
Network shape:    (4, 4)
Network nonzero:  0
```

A few things to note:

- `beta=0.0` is a required stub — add a `SEIR.Transmission` component and raise beta to a physiologically meaningful value (e.g. 0.3–0.5 for daily rate) before calling `model.run()`.
- `model.network` is all zeros because no `gravity_k/a/b/c` are in params. Add those keys (or set `model.network` manually) to enable spatial coupling.
- Because all agents are susceptible with no initial infections (`I=0` everywhere), any non-zero beta still produces no epidemic without seeding — use `seed_infections_in_patch(model, patch=0, ninfections=5)` from `laser.generic.utils` before `model.run()`.
