Here's the complete code:

```python
import numpy as np
from laser.core.migration import gravity, row_normalizer

# --- Patch setup ---
populations = np.array([50_000, 100_000, 75_000, 200_000])
n = len(populations)
spacing_km = 100.0  # patches in a line, 100 km apart

# --- Distance matrix (km) ---
# Patches at positions 0, 100, 200, 300 km along a line
positions = np.arange(n) * spacing_km
distances = np.abs(positions[:, None] - positions[None, :])

print("Distance matrix (km):")
print(distances)

# --- Gravity network ---
# M_{i,j} = k * p_i^a * p_j^b / d_{ij}^c
k, a, b, c = 0.01, 1, 1, 2
network = gravity(populations, distances, k, a, b, c)

print("\nRaw gravity matrix:")
print(np.array2string(network, precision=4, suppress_small=True))
print(f"\nRow sums (export fractions): {network.sum(axis=1)}")

# --- Normalize rows (cap at 20% export per patch) ---
network = row_normalizer(network, max_fraction=0.2)

print("\nNormalized gravity matrix:")
print(np.array2string(network, precision=6, suppress_small=True))
print(f"\nRow sums after normalization: {network.sum(axis=1)}")

# --- Validation (from Layer 1) ---
assert network.sum() > 0, "Network is all zeros — check gravity_k"
assert network.sum(axis=1).max() < 0.3, \
    f"Max row sum = {network.sum(axis=1).max():.3f} — may be too high"
print("\nValidation passed.")
```

**What this does:**

1. **Distance matrix** — Since patches are on a line with 100 km spacing, `d(i,j) = |i - j| * 100`. This gives a 4x4 symmetric matrix with zeros on the diagonal and distances of 100, 200, 300 km.

2. **`gravity()`** — Computes raw flows as `M_{i,j} = k * p_i^a * p_j^b / d_{i,j}^c`. With `c=2`, flow decays as the square of distance, so nearby large patches dominate.

3. **`row_normalizer()`** — Caps each row sum at 20% so no patch exports more than 20% of its population per tick. The diagonal stays zero (no self-migration).

**Expected output structure:**
- Patch 3 (pop 200k) attracts the most inflow due to its large population
- Nearest-neighbor flows dominate because `c=2` means rapid distance decay
- The 50k→200k flow at 300 km distance is ~44x weaker than at 100 km (distance ratio squared: `(300/100)^2 = 9`, further reduced by the `1/d^2` term)
