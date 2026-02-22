The filesystem is read-only, so I'll present the complete code directly:

```python
"""
Six-patch GeoDataFrame scenario for a LASER spatial SEIR model.

Each patch represents a distinct geographic community node. The scenario
DataFrame follows LASER conventions (laser-core / laser-measles):
  - population      : total patch population
  - S, E, I, R      : initial compartment counts (must sum to population)
  - latitude/longitude : patch centroid coordinates
  - geometry        : shapely Point (makes this a GeoDataFrame)

Initial conditions: 98% S, 1% E, 1% I, 0% R.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# ---------------------------------------------------------------------------
# 1. Define the six patches
# ---------------------------------------------------------------------------

PATCH_DATA = {
    "patch_id":   [0,       1,       2,       3,       4,        5],
    "name":       ["Alpha", "Beta",  "Gamma", "Delta", "Epsilon","Zeta"],
    # Populations spaced between 10,000 and 100,000
    "population": [10_000,  25_000,  40_000,  58_000,  80_000,  100_000],
    # Approximate lat/lon coordinates (representative US cities)
    "latitude":   [41.85,   38.63,   44.98,   40.80,   39.74,    33.45],
    "longitude":  [-87.65,  -90.20,  -93.27,  -96.68, -104.98, -112.07],
}

df = pd.DataFrame(PATCH_DATA)

# ---------------------------------------------------------------------------
# 2. Compute initial SEIR compartment counts
# ---------------------------------------------------------------------------
# Floor E and I to avoid fractional counts; assign the remainder to S
# so that S + E + I + R == population exactly for every patch.

pop = df["population"].values

df["E"] = np.floor(0.01 * pop).astype(int)
df["I"] = np.floor(0.01 * pop).astype(int)
df["R"] = 0
df["S"] = pop - df["E"] - df["I"] - df["R"]

# Sanity check
assert (df["S"] + df["E"] + df["I"] + df["R"] == df["population"]).all(), \
    "Compartment counts do not sum to population for all patches."

# ---------------------------------------------------------------------------
# 3. Attach Point geometry and promote to GeoDataFrame
# ---------------------------------------------------------------------------

geometry = [
    Point(lon, lat)
    for lon, lat in zip(df["longitude"], df["latitude"])
]

scenario: gpd.GeoDataFrame = gpd.GeoDataFrame(
    df,
    geometry=geometry,
    crs="EPSG:4326",   # WGS-84 geographic CRS
)

# Reorder columns for readability
scenario = scenario[
    ["patch_id", "name", "population", "S", "E", "I", "R",
     "latitude", "longitude", "geometry"]
]

# ---------------------------------------------------------------------------
# 4. Display the scenario
# ---------------------------------------------------------------------------

print("=" * 70)
print("LASER 6-Patch Spatial SEIR Scenario")
print("=" * 70)
print(scenario.to_string(index=False))
print()
print("CRS:", scenario.crs)
print("Geometry type:", scenario.geometry.geom_type.unique().tolist())
print()

print("-" * 70)
print(f"{'Patch':<10} {'Pop':>8}  {'S%':>7}  {'E%':>7}  {'I%':>7}  {'R%':>7}")
print("-" * 70)
for _, row in scenario.iterrows():
    n = row["population"]
    print(
        f"{row['name']:<10} {n:>8,}  "
        f"{100*row['S']/n:>6.2f}%  "
        f"{100*row['E']/n:>6.2f}%  "
        f"{100*row['I']/n:>6.2f}%  "
        f"{100*row['R']/n:>6.2f}%"
    )
print("-" * 70)
print(f"{'TOTAL':<10} {scenario['population'].sum():>8,}")

# ---------------------------------------------------------------------------
# 5. Optional: load into laser_core LaserFrame (if installed)
# ---------------------------------------------------------------------------

try:
    from laser_core.laserframe import LaserFrame

    npatches = len(scenario)
    patches = LaserFrame(npatches)

    patches.add_scalar_property("population", dtype=np.int64)
    patches.add_scalar_property("S",           dtype=np.int64)
    patches.add_scalar_property("E",           dtype=np.int64)
    patches.add_scalar_property("I",           dtype=np.int64)
    patches.add_scalar_property("R",           dtype=np.int64)
    patches.add_scalar_property("latitude",    dtype=np.float64)
    patches.add_scalar_property("longitude",   dtype=np.float64)

    patches.population[:] = scenario["population"].values
    patches.S[:]          = scenario["S"].values
    patches.E[:]          = scenario["E"].values
    patches.I[:]          = scenario["I"].values
    patches.R[:]          = scenario["R"].values
    patches.latitude[:]   = scenario["latitude"].values
    patches.longitude[:]  = scenario["longitude"].values

    print(f"\nlaser_core LaserFrame initialised: {patches.count} patches")
    print(f"  S counts: {patches.S[:]}")

except ImportError:
    print(
        "\nlaser_core not installed. The GeoDataFrame `scenario` is complete\n"
        "and ready to pass to a LASER Model constructor."
    )
```

---

**Expected output for the summary table:**

```
----------------------------------------------------------------------
Patch          Pop      S%      E%      I%      R%
----------------------------------------------------------------------
Alpha       10,000   98.00%   1.00%   1.00%   0.00%
Beta        25,000   98.00%   1.00%   1.00%   0.00%
Gamma       40,000   98.00%   1.00%   1.00%   0.00%
Delta       58,000   98.00%   1.00%   1.00%   0.00%
Epsilon     80,000   98.00%   1.00%   1.00%   0.00%
Zeta       100,000   98.00%   1.00%   1.00%   0.00%
----------------------------------------------------------------------
TOTAL      313,000
```

---

**Key design decisions:**

| Point | Rationale |
|---|---|
| `np.floor` for E and I | Avoids fractional agents; remainder goes to S so counts sum exactly to `population` |
| `crs="EPSG:4326"` | WGS-84 is the standard for lat/lon; required by LASER's `scenario.plot()` and distance utilities (`calc_distances`) |
| `geometry = Point(lon, lat)` | Shapely convention is `(x, y)` = `(longitude, latitude)` |
| `LaserFrame` block is try/except | Makes the script self-contained — the GeoDataFrame itself is the scenario artifact; LASER integration is additive |
| Columns match LASER conventions | `population`, `latitude`, `longitude`, `geometry` are the names referenced directly in `laser_measles.model` and `laser_cholera.metapop.model` source code |

Sources:
- [InstituteforDiseaseModeling/laser](https://github.com/InstituteforDiseaseModeling/laser)
- [laser_measles.model — laser_measles 0.5.0 docs](https://docs.idmod.org/projects/laser-measles/en/latest/_modules/laser_measles/model.html)
- [laser_cholera.metapop.model — LASIK 0.8.0 docs](https://docs.idmod.org/projects/laser-cholera/en/v0.8.0/_modules/laser_cholera/metapop/model.html)
