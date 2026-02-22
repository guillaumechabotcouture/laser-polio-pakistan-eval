Here is how to compute a gravity migration network for 4 patches and assign it to `model.network`, based on the actual LASER framework source:

```python
import numpy as np

# --- LASER imports (available after: pip install laser-core laser-generic) ---
from laser.core.migration import distance, gravity, row_normalizer
# from laser.generic.model import Model  # uncomment when using a full Model

# =============================================================================
# 1. Define 4 patches: populations and geographic centroids (lat, lon)
# =============================================================================

populations = np.array([500_000, 1_200_000, 300_000, 800_000], dtype=np.float64)

# Example centroids (latitude, longitude in decimal degrees)
lats = np.array([40.71, 34.05, 41.88, 29.76])   # New York, LA, Chicago, Houston
lons = np.array([-74.01, -118.24, -87.63, -95.37])

# =============================================================================
# 2. Compute the NxN pairwise distance matrix (km) via Haversine formula
# =============================================================================
# laser.core.migration.distance() returns a symmetric NxN matrix.
# Passing the same arrays for both origin and destination gives the full matrix.

dist_matrix = distance(lats, lons, lats, lons)   # shape: (4, 4), units: km

print("Distance matrix (km):")
print(np.round(dist_matrix, 1))

# =============================================================================
# 3. Compute the raw gravity network
#
#   network[i, j] = k * (pop[i]^a) * (pop[j]^b) * (dist[i,j]^(-c))
#
#   k : overall scaling constant
#   a : origin  population exponent
#   b : destination population exponent
#   c : distance-decay exponent
#   Diagonal is set to 0 (no self-migration).
# =============================================================================

k, a, b, c = 500.0, 1.0, 1.0, 2.0

raw_network = gravity(populations, dist_matrix, k=k, a=a, b=b, c=c)

print("\nRaw gravity network (flows):")
print(np.round(raw_network, 4))

# =============================================================================
# 4. Row-normalize so that no patch exports more than max_fraction
#    of its population per tick.
#    row_normalizer() rescales any row whose sum exceeds max_rowsum.
# =============================================================================

max_rowsum = (1/16) + (1/32)   # ≈ 0.0938 (matches laser-generic default)

network = row_normalizer(raw_network, max_rowsum)

print("\nNormalized gravity network (migration rates per tick):")
print(np.round(network, 6))
print("\nRow sums (should all be <= {:.4f}):".format(max_rowsum))
print(np.round(network.sum(axis=1), 6))

# =============================================================================
# 5. Assign to model.network
#
#    The Model constructor computes its own network automatically from the
#    scenario GeoDataFrame.  To override it with a custom gravity network,
#    simply assign after construction.
# =============================================================================

# model = Model(scenario, params)   # normal construction
# model.network = network           # replace with manually computed network

# Standalone demo (no Model instance required):
class _StubModel:
    """Minimal stub showing the attribute contract."""
    network: np.ndarray = None

model = _StubModel()
model.network = network

print("\nmodel.network shape:", model.network.shape)
print("model.network dtype:", model.network.dtype)
print("\nFull model.network:")
print(model.network)
```

---

Key points:

| Concern | Detail |
|---|---|
| `distance()` | Haversine great-circle distance; diagonal is set to 1 internally to avoid division-by-zero, then zeroed |
| `gravity()` | `network[i,j] = k·pop_i^a·pop_j^b·dist_ij^(-c)`; diagonal forced to 0 |
| `row_normalizer()` | Scales any row whose sum exceeds `max_rowsum` so the value represents a per-tick emigration fraction in [0, 1] |
| `model.network` | Plain `np.ndarray` of shape `(n_patches, n_patches)` — directly assignable after `Model.__init__` |
| Default params | `gravity_k=500`, `gravity_a=1`, `gravity_b=1`, `gravity_c=2` (overridable via `params` object) |

Special-case parameterizations from the LASER docs:
- `a=0` → Xia's model (destination-only gravity, no origin-size dependence)
- `a=0, b=0` → pure spatial diffusion (distance decay only)
- `c=0, a=1, b=1` → mean-field / fully-mixed model (distance ignored)
