Here's the complete code:

```python
import numpy as np
from laser.core.migration import gravity, row_normalizer, distance

# --- 4 patches: name, population, lat, lon ---
patches = {
    "name":       ["CityA", "CityB", "TownC", "TownD"],
    "population": [5_000_000, 2_000_000, 500_000, 100_000],
    "lat":        [33.69, 31.55, 34.01, 30.20],
    "lon":        [73.04, 74.35, 71.58, 71.47],
}

populations = np.array(patches["population"], dtype=np.float64)
lats = np.array(patches["lat"])
lons = np.array(patches["lon"])

# Pairwise Haversine distances (km) — LASER built-in
distances = distance(lats, lons)
print("Distance matrix (km):\n", np.round(distances, 1))

# Raw gravity network: M_{i,j} = k * p_i^a * p_j^b / d_{ij}^c
# a=0 (source pop doesn't affect outflow), b=0.5, c=1.5
k, a, b, c = 1.0, 0, 0.5, 1.5
network = gravity(populations, distances, k, a, b, c)

print("\nRaw gravity matrix (before normalization):")
print(np.round(network, 6))
print("Row sums:", np.round(network.sum(axis=1), 4))

# Row-normalize so no patch exports more than 20% of its FOI
network = row_normalizer(network, max_fraction=0.2)

print("\nNormalized gravity matrix:")
print(np.round(network, 6))
print("Row sums:", np.round(network.sum(axis=1), 4))

# Verify: no row exceeds the 20% cap
assert network.sum(axis=1).max() <= 0.2 + 1e-9, \
    f"Row sum exceeds cap: {network.sum(axis=1).max():.4f}"
print("\nAll row sums <= 0.20: PASS")
```

**What this does:**

1. **`distance(lats, lons)`** — computes the 4x4 Haversine distance matrix in km (diagonal = 0).
2. **`gravity(populations, distances, k, a, b, c)`** — produces the raw gravity flow matrix. With `a=0`, source population doesn't affect outflow; only destination size (`b=0.5`) and distance (`c=1.5`) matter.
3. **`row_normalizer(network, max_fraction=0.2)`** — scales down each row so that `row.sum() <= 0.2`. Rows already below the cap are left unchanged. This ensures no patch exports more than 20% of its force of infection to other patches.

The row sums represent the fraction of each patch's FOI that "leaks" to neighbors. Capping at 20% keeps most transmission local while still allowing spatial coupling (traveling waves, reimportation into small patches).
