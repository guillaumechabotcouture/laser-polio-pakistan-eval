Here is a focused example showing how to compute a gravity network for 4 patches and assign it to `model.network` after construction.

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.core.migration import gravity, row_normalizer, distance
from laser.generic import Model

# ── 1. Build a 4-patch scenario ──────────────────────────────────────────────
# Coordinates (lat, lon) for 4 cities in Pakistan
lats = np.array([33.72, 31.55, 24.86, 34.01])
lons = np.array([73.04, 74.35, 67.01, 71.57])
pops = np.array([2_200_000, 13_000_000, 16_500_000, 680_000], dtype=np.int64)

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(4),
        "name":       ["Islamabad", "Lahore", "Karachi", "Peshawar"],
        "population": pops,
        "S": pops - 10,
        "E": np.zeros(4, dtype=int),
        "I": np.full(4, 10, dtype=int),
        "R": np.zeros(4, dtype=int),
    },
    geometry=[Point(lon, lat) for lat, lon in zip(lats, lons)],
    crs="EPSG:4326",
)

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all()

# ── 2. Params — deliberately omit gravity_* so Model does NOT auto-build ──────
params = PropertySet({
    "prng_seed": 42,
    "nticks":    365,
    "beta":      0.3,
})

# ── 3. Instantiate the model ──────────────────────────────────────────────────
model = Model(scenario, params)

# ── 4. Compute pairwise Haversine distances (km) ─────────────────────────────
n = len(scenario)
dist_matrix = distance(lats, lons, lats, lons)   # (4, 4) matrix, diagonal = 0

# Avoid self-coupling: set diagonal to inf so gravity gives 0 self-flow
np.fill_diagonal(dist_matrix, np.inf)

# ── 5. Compute raw gravity flows: M_{i,j} = 1 * pop_j^b / d_{ij}^c ──────────
#   k=1  → we normalise manually afterward
#   a=0  → source population does not scale outward flow (LASER convention)
#   b=0.5, c=1.5 → typical measles / polio parameters
raw = gravity(pops, dist_matrix, k=1, a=0, b=0.5, c=1.5)

# ── 6. Normalise so gravity_k represents average export fraction ──────────────
gravity_k = 0.02                            # ~2 % of residents travel per day
avg_export = raw.sum(axis=1).mean()
if avg_export > 0:
    raw = raw / avg_export * gravity_k

# Cap maximum row-sum at 20 % (prevents unrealistic mass movement)
network = row_normalizer(raw, max_fraction=0.2)

# ── 7. Assign to the model ───────────────────────────────────────────────────
model.network = network

# ── 8. Verify ────────────────────────────────────────────────────────────────
assert model.network.sum() > 0, \
    "Network is all zeros — check gravity_k or distance matrix"
assert model.network.sum(axis=1).max() < 0.3, (
    f"Max row-sum {model.network.sum(axis=1).max():.3f} too high — "
    "agents may exceed patch capacity"
)
assert np.allclose(np.diag(model.network), 0), \
    "Diagonal should be zero (no self-migration)"

print("network shape:", model.network.shape)
print("row sums (export fractions):\n", model.network.sum(axis=1).round(4))
print("\nnetwork matrix:\n", np.round(model.network, 5))
```

**Key decisions explained:**

| Choice | Why |
|--------|-----|
| Omit `gravity_k/a/b/c` from params | Prevents `Model.__init__()` from auto-building the network, so the manual assignment is the authoritative one |
| `k=1`, then manual rescale | Lets you set `gravity_k` to mean "average daily export fraction" directly, which is easier to reason about than raw gravity units |
| `row_normalizer(..., 0.2)` | Caps each row at 20 % export so no patch loses more than a fifth of its residents per day |
| `np.fill_diagonal(dist_matrix, np.inf)` | The `gravity()` function divides by distance; zero distance produces `inf` flows on the diagonal, which `row_normalizer` would cap badly — explicit fill is cleaner |

After assigning, `SEIR.Transmission` (i.e., `TransmissionSE`) reads `model.network` automatically during its force-of-infection computation — no further wiring is needed.
