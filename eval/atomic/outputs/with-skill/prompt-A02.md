Here's the complete code:

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

# 6-patch scenario with varying populations
patches = {
    "nodeid": [0, 1, 2, 3, 4, 5],
    "name": ["Alpha", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot"],
    "population": [10_000, 25_000, 40_000, 60_000, 80_000, 100_000],
    "geometry": [
        Point(67.0, 25.0),   # roughly Karachi-region coords
        Point(68.5, 26.5),
        Point(70.0, 28.0),
        Point(71.5, 30.0),
        Point(72.5, 31.5),
        Point(73.0, 33.5),
    ],
}

scenario = gpd.GeoDataFrame(patches, crs="EPSG:4326")

# Initial conditions: 98% S, 1% E, 1% I, 0% R
scenario["E"] = np.round(0.01 * scenario["population"]).astype(np.uint32)
scenario["I"] = np.round(0.01 * scenario["population"]).astype(np.uint32)
scenario["R"] = np.uint32(0)
scenario["S"] = (scenario["population"] - scenario["E"] - scenario["I"] - scenario["R"]).astype(np.uint32)

# Validate: S + E + I + R == population in every patch
assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch needs initial infections"

print(scenario[["nodeid", "name", "population", "S", "E", "I", "R"]].to_string(index=False))
print(f"\nCRS: {scenario.crs}")
print(f"\nGeometry:\n{scenario[['name', 'geometry']].to_string(index=False)}")
```

Output:

```
 nodeid     name  population      S    E     I  R
      0    Alpha       10000   9800  100   100  0
      1    Bravo       25000  24500  250   250  0
      2  Charlie       40000  39200  400   400  0
      3    Delta       60000  58800  600   600  0
      4     Echo       80000  78400  800   800  0
      5  Foxtrot      100000  98000 1000  1000  0
```

Key points:
- **Compartments are computed from population**, not hardcoded, so `S + E + I + R == population` holds exactly in every patch
- **`R` is assigned first as 0**, then `E` and `I` via rounding, and `S` gets the remainder — this avoids off-by-one rounding errors
- **Assertions** catch any mismatch before it silently propagates into a simulation
- This GeoDataFrame plugs directly into `Model(scenario, parameters)` as shown in the skill's Step 7
