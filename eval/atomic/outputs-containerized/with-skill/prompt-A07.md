Here is a complete, self-contained script demonstrating `row_normalizer` on a 4-patch gravity network.

```python
import numpy as np
from laser.core.migration import gravity, row_normalizer, distance

# ------------------------------------------------------------------
# 1. Define 4 patches: populations and geographic coordinates
# ------------------------------------------------------------------
patch_names = ["Karachi", "Lahore", "Peshawar", "Quetta"]
populations = np.array([15_000_000, 13_000_000, 2_200_000, 1_100_000], dtype=float)

# (latitude, longitude) in decimal degrees
lats = np.array([24.86, 31.55, 34.01, 30.19])
lons = np.array([67.01, 74.34, 71.57, 66.99])

# ------------------------------------------------------------------
# 2. Compute pairwise Haversine distance matrix (km)
# ------------------------------------------------------------------
# distance(lats, lons, lats, lons) returns a (4, 4) matrix
distances = distance(lats, lons, lats, lons)

# Set diagonal to 1 to avoid divide-by-zero in gravity()
np.fill_diagonal(distances, 1.0)

print("Distance matrix (km):")
for i, name in enumerate(patch_names):
    row = "  ".join(f"{distances[i, j]:8.1f}" for j in range(len(patch_names)))
    print(f"  {name:<12} {row}")

# ------------------------------------------------------------------
# 3. Build gravity network  M_{i,j} = k * p_i^a * p_j^b / d_{ij}^c
#    Convention: a=0 (source pop does not affect outward flow)
# ------------------------------------------------------------------
k, a, b, c = 1.0, 0.0, 0.5, 1.5
raw_network = gravity(populations, distances, k, a, b, c)

# Scale so k=1 means ~average export fraction of 1.0 before normalization,
# then apply a coupling constant (gravity_k) to set the actual export level.
gravity_k = 0.05  # 5% average export fraction target
avg_export = np.mean(raw_network.sum(axis=1))
network = raw_network / avg_export * gravity_k

print("\n--- BEFORE row_normalizer ---")
row_sums = network.sum(axis=1)
for i, name in enumerate(patch_names):
    print(f"  {name:<12}  row_sum = {row_sums[i]:.4f}  "
          f"({row_sums[i]*100:.1f}% of FOI exported)")

# ------------------------------------------------------------------
# 4. Apply row_normalizer: cap each row at 20% (max_fraction=0.2)
#
#    row_normalizer(network, max_fraction) scales down any row whose
#    sum exceeds max_fraction so that sum == max_fraction exactly.
#    Rows already below the cap are left unchanged.
# ------------------------------------------------------------------
normalized_network = row_normalizer(network, max_fraction=0.2)

print("\n--- AFTER row_normalizer (cap = 20%) ---")
row_sums_norm = normalized_network.sum(axis=1)
for i, name in enumerate(patch_names):
    capped = row_sums[i] > 0.2
    flag = " <- capped" if capped else ""
    print(f"  {name:<12}  row_sum = {row_sums_norm[i]:.4f}  "
          f"({row_sums_norm[i]*100:.1f}% of FOI exported){flag}")

# ------------------------------------------------------------------
# 5. Full normalized matrix
# ------------------------------------------------------------------
print("\nNormalized network matrix (M[i,j] = export fraction from i to j):")
header = "             " + "  ".join(f"{n:>12}" for n in patch_names)
print(header)
for i, src in enumerate(patch_names):
    row = "  ".join(f"{normalized_network[i, j]:12.6f}"
                    for j in range(len(patch_names)))
    print(f"  {src:<12} {row}")

# ------------------------------------------------------------------
# 6. Validate: no row exceeds 20%
# ------------------------------------------------------------------
max_export = normalized_network.sum(axis=1).max()
assert max_export <= 0.2 + 1e-9, \
    f"Row normalization failed: max export = {max_export:.4f}"
assert normalized_network.sum() > 0, \
    "Network is all zeros — check gravity_k or distances"
print(f"\nValidation passed: max export fraction = {max_export:.4f} (<= 0.20)")
```

**What each step does:**

| Step | Purpose |
|------|---------|
| `distance(lats, lons, lats, lons)` | Builds a 4x4 Haversine distance matrix (km) |
| `gravity(populations, distances, k=1, a=0, b=0.5, c=1.5)` | Computes raw flow matrix $M_{i,j} = p_j^{0.5} / d_{ij}^{1.5}$ |
| Manual scaling by `gravity_k / avg_export` | Maps `k=1` to a meaningful coupling constant (5% here) |
| `row_normalizer(network, max_fraction=0.2)` | Scales down any row whose sum exceeds 0.2; leaves rows below the cap unchanged |

**`row_normalizer` semantics:** it divides each row by `row_sum / max_fraction` only when `row_sum > max_fraction`. The diagonal (self-coupling) stays at zero — only cross-patch export fractions are constrained. After normalization, each row sum represents the fraction of a patch's local force of infection that is redistributed to neighboring patches, capped at 20%.
