Here's the complete code:

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.generic import Model

# 4-patch scenario: each 50,000 population, all susceptible
scenario = gpd.GeoDataFrame(
    {
        "nodeid": [0, 1, 2, 3],
        "name": ["Patch_A", "Patch_B", "Patch_C", "Patch_D"],
        "population": np.uint32([50_000, 50_000, 50_000, 50_000]),
        "S": np.uint32([50_000, 50_000, 50_000, 50_000]),
        "E": np.uint32([0, 0, 0, 0]),
        "I": np.uint32([0, 0, 0, 0]),
        "R": np.uint32([0, 0, 0, 0]),
    },
    geometry=[Point(0, 0), Point(1, 0), Point(0, 1), Point(1, 1)],
    crs="EPSG:4326",
)

# Validate initial conditions
assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all()

# Parameters
params = PropertySet({"nticks": 365, "prng_seed": 20250101})

# Create model
model = Model(scenario, params)
```

A few notes:

- **All susceptible** means `S = population` and `E = I = R = 0`. No initial infections, so nothing will happen dynamically unless you later add importation or seed some infections.
- The `geometry` points are arbitrary placeholders — replace with real coordinates if you need a gravity-based migration network.
- `Model.__init__` calls `calc_capacity()` internally to pre-allocate agent slots based on the scenario populations.
