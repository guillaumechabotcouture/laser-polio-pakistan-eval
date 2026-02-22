Here is the complete code:

```python
import numpy as np
from laser.core.migration import distance, gravity, row_normalizer

# ── Patch setup ─────────────────────────────────────────────────────────────
populations = np.array([50_000, 100_000, 75_000, 200_000], dtype=float)
n_patches = len(populations)

# Four patches in a line, 100 km apart along a meridian.
# 1 degree latitude ≈ 111.32 km → 100 km ≈ 0.8983 degrees.
km_per_deg = 111.32
spacing_deg = 100.0 / km_per_deg

lats = np.array([i * spacing_deg for i in range(n_patches)])
lons = np.zeros(n_patches)

# ── Distance matrix (Haversine) ──────────────────────────────────────────────
# distance(lats, lons, lats, lons) returns an (n x n) pairwise matrix in km.
dist_matrix = distance(lats, lons, lats, lons)

print("Distance matrix (km):")
print(np.round(dist_matrix, 1))
# Expected:
#   [[  0.  100.  200.  300.]
#    [100.    0.  100.  200.]
#    [200.  100.    0.  100.]
#    [300.  200.  100.    0.]]

# ── Gravity network ──────────────────────────────────────────────────────────
# Formula: M_{i,j} = k * pop_i^a * pop_j^b / d_{ij}^c
# With k=0.01, a=1, b=1, c=2 → M_{i,j} = 0.01 * pop_i * pop_j / d_{ij}^2
# Diagonal (i==i, d=0) is set to 0 by the LASER implementation.

k, a, b, c = 0.01, 1.0, 1.0, 2.0

network = gravity(populations, dist_matrix, k, a, b, c)

print("\nRaw gravity network (flow fractions):")
print(network)

print("\nRow sums (total outflow fraction per patch):")
print(network.sum(axis=1))

# ── Validation ───────────────────────────────────────────────────────────────
assert network.sum() > 0, "Network is all zeros — check gravity_k or distances"
assert network.sum(axis=1).max() < 0.3, (
    f"Max row sum {network.sum(axis=1).max():.3f} exceeds 30% — "
    "agents may leave patches faster than they can be replenished"
)
print("\nValidation passed.")

# ── Optional: cap outflow per patch at 20% ───────────────────────────────────
network_capped = row_normalizer(network, max_fraction=0.2)

print("\nCapped network (max 20% outflow per patch):")
print(network_capped)
print("Row sums after capping:")
print(network_capped.sum(axis=1))
```

**Key points:**

| Item | Detail |
|---|---|
| `distance(lats, lons, lats, lons)` | Returns full (n×n) Haversine matrix in km; diagonal is 0 |
| `gravity(pops, dists, k, a, b, c)` | Computes M_{i,j} = k·pop_i^a·pop_j^b / d_{i,j}^c; diagonal is 0 |
| `row_normalizer(net, 0.2)` | Caps each row sum at 20% so no patch loses >20% of its population per tick |
| `a=1` here | Source population scales outflow — larger patches send *more* migrants |
| `a=0` convention | LASER's default convention; source pop has no effect on outward flow |

With `c=2` the network is distance-squared decay, so patches 300 km apart (patches 0↔3) are coupled 9× more weakly than adjacent patches at 100 km. The large destination population (patch 3, 200k) pulls disproportionately high inflow from all sources.
