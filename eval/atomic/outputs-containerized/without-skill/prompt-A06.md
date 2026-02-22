Here is the complete code and an explanation of each step.

---

## Complete code

```python
import numpy as np

# ---------------------------------------------------------------------------
# Patch configuration
# ---------------------------------------------------------------------------
populations = np.array([50_000, 100_000, 75_000, 200_000], dtype=np.float64)
n_patches   = len(populations)

# Patches in a line, 100 km apart → positions 0, 100, 200, 300 km
positions_km = np.arange(n_patches) * 100.0

# Gravity model parameters
k = 0.01   # scaling constant
a = 1.0    # source population exponent
b = 1.0    # destination population exponent
c = 2.0    # distance friction exponent

# ---------------------------------------------------------------------------
# Step 1: Distance matrix  D[i,j] = |pos_i - pos_j|
# ---------------------------------------------------------------------------
distances = np.abs(positions_km[:, np.newaxis] - positions_km[np.newaxis, :])

# ---------------------------------------------------------------------------
# Step 2: Gravity flow matrix
#   F[i,j] = k * P_i^a * P_j^b / D[i,j]^c    (i != j)
#   F[i,i] = 0
# ---------------------------------------------------------------------------
with np.errstate(divide="ignore", invalid="ignore"):
    gravity = (
        k
        * np.outer(populations**a, populations**b)
        / distances**c
    )
np.fill_diagonal(gravity, 0.0)

# ---------------------------------------------------------------------------
# Step 3: Row-normalised migration rate matrix
#   rate[i,j] = F[i,j] / sum_j(F[i,j])  → fraction of patch i leaving to j
# ---------------------------------------------------------------------------
row_sums       = gravity.sum(axis=1, keepdims=True)
migration_rates = np.where(row_sums > 0, gravity / row_sums, 0.0)
```

---

## Output

**Distance matrix D (km)**
```
[[  0, 100, 200, 300],
 [100,   0, 100, 200],
 [200, 100,   0, 100],
 [300, 200, 100,   0]]
```

**Gravity flow matrix F** — `F[i,j] = 0.01 * P_i * P_j / d_ij²`
```
[[    0,    5000,   937.5,  1111.1],
 [ 5000,       0,  7500,    5000  ],
 [  937.5,  7500,     0,  15000  ],
 [ 1111.1,  5000, 15000,      0  ]]
```

**Row-normalised migration rates** (fraction of patch i moving to j)
```
[[0,      0.709, 0.133, 0.158],
 [0.286,  0,     0.429, 0.286],
 [0.040,  0.320, 0,     0.640],
 [0.053,  0.237, 0.711, 0    ]]
```

Total gravity flux: **69,097.2**  |  Max flow: **15,000** (patch 2 → 3)

---

## Key design notes

| Point | Explanation |
|---|---|
| **Formula** | `F_ij = k * P_i^a * P_j^b / d_ij^c` — the standard LASER gravity kernel |
| **`np.outer`** | Efficiently computes all `P_i^a * P_j^b` pairs in one call |
| **`errstate`** | Suppresses the divide-by-zero on the diagonal; `fill_diagonal(0)` then zeros it cleanly |
| **`c=2` (inverse-square)** | Distance decay is steep — flows between patches 0 and 3 (300 km) are ~16× smaller than patches 2 and 3 (100 km) |
| **Symmetric net flow** | With `a=b=1`, `F_ij = F_ji`, so net flow per patch is identically zero. Breaking symmetry (`a ≠ b`) or making `k` direction-dependent produces net migration |
| **Migration rate matrix** | Row-normalising gives the per-step fraction leaving each patch, which is what LASER's mobility layer ingests as the network kernel |
