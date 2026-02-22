It works. Here's the complete code:

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core import PropertySet
from laser.generic import Model

# Build a 4-patch scenario GeoDataFrame
n_patches = 4
pop_per_patch = 50_000

scenario = gpd.GeoDataFrame(
    {
        "nodeid": np.arange(n_patches),
        "population": np.full(n_patches, pop_per_patch),
    },
    geometry=[
        Point(0.0, 0.0),
        Point(1.0, 0.0),
        Point(0.0, 1.0),
        Point(1.0, 1.0),
    ],
    crs="EPSG:4326",
)

# All susceptible — no initial infected columns needed (defaults to 0 infected)

# Create parameters
params = PropertySet({"nticks": 365, "prng_seed": 20250101})

# Create the model
model = Model(scenario, params)

print(f"Patches:          {model.nodes.count}")
print(f"Total population: {model.people.count}")
print(f"Ticks:            {model.params.nticks}")
print(f"PRNG seed:        {model.params.prng_seed}")
print(f"Model name:       {model.name}")
```

**Key details:**

- **Scenario GeoDataFrame** — requires `nodeid` (0-indexed patch IDs), `population`, and a geometry column with a CRS set (`EPSG:4326`). The model internally reprojects to EPSG:3857 for distance calculations and builds a gravity-based connectivity network between patches.
- **PropertySet** — wraps a dict into an attribute-accessible object. `nticks` and `prng_seed` are consumed directly by the Model constructor.
- **All susceptible** — with no initial infected counts in the scenario, the entire population of 200,000 (4 × 50,000) starts susceptible with the default SEIR state set.
