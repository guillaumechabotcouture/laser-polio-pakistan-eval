```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

# ── 1. Patch populations (10k → 100k) ──────────────────────────────────────
populations = np.array(
    [10_000, 25_000, 40_000, 55_000, 70_000, 100_000],
    dtype=np.uint32,
)
n_patches = len(populations)

# ── 2. Point geometry (lon, lat) ────────────────────────────────────────────
# Six illustrative locations arranged along a corridor
coordinates = [
    (67.01, 24.86),   # Patch 0 — smallest
    (67.35, 25.20),   # Patch 1
    (67.70, 25.55),   # Patch 2
    (68.05, 25.90),   # Patch 3
    (68.40, 26.25),   # Patch 4
    (68.75, 26.60),   # Patch 5 — largest
]

# ── 3. Initial compartment counts (98% S, 1% E, 1% I, 0% R) ────────────────
# Round E and I; derive S as remainder so S+E+I+R == population exactly.
E = np.round(0.01 * populations).astype(np.uint32)
I = np.round(0.01 * populations).astype(np.uint32)
R = np.zeros(n_patches, dtype=np.uint32)
S = populations - E - I - R          # guarantees exact sum

# ── 4. Assemble GeoDataFrame ────────────────────────────────────────────────
scenario = gpd.GeoDataFrame(
    {
        "nodeid":      np.arange(n_patches, dtype=np.uint32),
        "name":        [f"Patch_{i}" for i in range(n_patches)],
        "population":  populations,
        "S":           S,
        "E":           E,
        "I":           I,
        "R":           R,
        "geometry":    [Point(lon, lat) for lon, lat in coordinates],
    },
    crs="EPSG:4326",
)

# ── 5. Validation ────────────────────────────────────────────────────────────
assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "S + E + I + R must equal population in every patch"
assert (scenario.I > 0).all(), \
    "Every patch must have at least one infectious individual"
assert len(scenario) == 6, "Expected exactly 6 patches"

# ── 6. Summary ───────────────────────────────────────────────────────────────
print(scenario[["nodeid", "name", "population", "S", "E", "I", "R"]].to_string(index=False))

fracs = scenario[["S", "E", "I", "R"]].div(scenario["population"], axis=0)
print("\nCompartment fractions (per patch):")
print(fracs.round(4).to_string(index=False))

print(f"\nCRS: {scenario.crs}")
print(f"Geometry type: {scenario.geometry.geom_type.unique().tolist()}")
```

**Key design decisions:**

- **Population** — six steps from 10 k to 100 k: `[10k, 25k, 40k, 55k, 70k, 100k]`.
- **Compartments** — `E` and `I` are computed with `np.round`, then `S = population − E − I − R`. This ensures `S + E + I + R == population` exactly in every row, avoiding off-by-one errors from independent rounding.
- **Geometry** — `Point(lon, lat)` in `EPSG:4326`; six points spaced ~0.5° apart along a northeast corridor.
- **Column names** — `nodeid`, `name`, `population`, `S`, `E`, `I`, `R`, `geometry` match what `laser_generic.Model` expects when ingesting a scenario GeoDataFrame (see skill Step 2).

For the smallest patch (10 k), the fractions work out exactly: `E = I = 100`, `S = 9 800` (98.0 % / 1.0 % / 1.0 % / 0 %). Larger patches (multiples of 100) also round cleanly.
