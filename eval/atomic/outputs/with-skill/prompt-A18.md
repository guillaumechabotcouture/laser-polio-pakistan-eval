Here's how to compute a gravity network for 4 patches and assign it to the model after construction:

```python
import numpy as np
from laser.core.migration import gravity, row_normalizer, distance

# --- 4 patches: lat/lon coordinates and populations ---
lats = np.array([30.2, 31.5, 25.4, 33.7])
lons = np.array([71.0, 74.3, 68.4, 73.1])
populations = np.array([5_000_000, 3_000_000, 2_000_000, 1_000_000])

# Haversine distance matrix (km) between all patch pairs
distances = distance(lats, lons)  # shape (4, 4), diagonal = 0

# Gravity parameters
gravity_k = 0.01   # overall coupling strength
gravity_b = 0.5    # destination population exponent
gravity_c = 1.5    # distance decay exponent

# Compute raw gravity matrix: M_ij = pop_j^b / d_ij^c
network = gravity(populations, distances, 1, 0, gravity_b, gravity_c)

# Normalize so gravity_k controls the average export fraction directly
avg_export = np.mean(network.sum(axis=1))
network = network / avg_export * gravity_k

# Cap per-patch export at 20% to prevent unrealistic mass migration
network = row_normalizer(network, 0.2)

# Assign to the model (after Model() construction)
model.network = network

# --- Verify ---
assert model.network.shape == (4, 4)
assert model.network.sum() > 0, "Network is all zeros — check gravity_k"
assert model.network.sum(axis=1).max() < 0.3, \
    f"Max row sum too high: {model.network.sum(axis=1).max():.3f}"

print("Row sums (export fractions per patch):", model.network.sum(axis=1))
print("Network:\n", np.round(model.network, 6))
```

Key points:

- **`gravity()`** computes raw flows as `pop_j^b / d_ij^c`. The `1, 0` arguments set `k=1` and `a=0` (source population exponent — convention is to leave it at 0).
- **`distance()`** returns a Haversine distance matrix in km from lat/lon arrays.
- **`row_normalizer(network, max_frac)`** clips each row so no patch exports more than `max_frac` of its population per tick.
- **`gravity_k`** controls overall coupling strength. Too small (< 1e-4) and patches evolve independently with no traveling waves. Too large (> 0.1) and you get unrealistic mass movement.
- If you include `gravity_k`, `gravity_b`, and `gravity_c` in the model's `PropertySet` params, `Model.__init__()` will auto-compute the network — manual assignment is only needed when you want custom normalization like above.
