The filesystem is read-only, so here is the complete code:

```python
"""
Apply row normalization to a 4-patch gravity network matrix using LASER,
capping each patch's export fraction at 20% of its force of infection.
"""

import numpy as np
from laser.core.migration import row_normalizer

# ---------------------------------------------------------------------------
# 4-patch gravity network matrix
#
# Entry [i, j] is the fraction of patch i's force of infection that is
# exported to patch j (off-diagonal). Diagonal entries are 0 (no
# self-coupling term included here).
#
# Two patches intentionally exceed the 20% export cap so the effect
# of row_normalizer is clearly visible.
# ---------------------------------------------------------------------------
network = np.array(
    [
        #  p0     p1     p2     p3      row sum
        [0.000, 0.120, 0.150, 0.080],  # patch 0: 0.350  (> 0.20 → rescaled)
        [0.080, 0.000, 0.060, 0.030],  # patch 1: 0.170  (≤ 0.20 → unchanged)
        [0.100, 0.180, 0.000, 0.140],  # patch 2: 0.420  (> 0.20 → rescaled)
        [0.050, 0.070, 0.040, 0.000],  # patch 3: 0.160  (≤ 0.20 → unchanged)
    ],
    dtype=np.float64,
)

print("=== Raw gravity network matrix ===")
print(network)
print()

print("Row sums before normalization:")
for i, s in enumerate(network.sum(axis=1)):
    status = "  <-- exceeds 20%, will be rescaled" if s > 0.20 else ""
    print(f"  Patch {i}: {s:.4f}{status}")

# ---------------------------------------------------------------------------
# row_normalizer(network, max_rowsum)
#
# For every row whose sum exceeds max_rowsum, the entire row is multiplied
# by (max_rowsum / row_sum), preserving relative weights while capping
# the total export fraction. Rows already within the limit are untouched.
# max_rowsum must be in [0, 1].
# ---------------------------------------------------------------------------
max_export = 0.20  # no patch may export more than 20% of its FoI

normalized = row_normalizer(network, max_rowsum=max_export)

print()
print(f"=== Normalized network (max_rowsum={max_export}) ===")
print(normalized)
print()

print("Row sums after normalization:")
for i, s in enumerate(normalized.sum(axis=1)):
    print(f"  Patch {i}: {s:.4f}")

print()
print("Verification – no row sum exceeds 20%:",
      bool((normalized.sum(axis=1) <= max_export + 1e-6).all()))
```

**Expected output:**

```
=== Raw gravity network matrix ===
[[0.    0.12  0.15  0.08 ]
 [0.08  0.    0.06  0.03 ]
 [0.1   0.18  0.    0.14 ]
 [0.05  0.07  0.04  0.   ]]

Row sums before normalization:
  Patch 0: 0.3500  <-- exceeds 20%, will be rescaled
  Patch 1: 0.1700
  Patch 2: 0.4200  <-- exceeds 20%, will be rescaled
  Patch 3: 0.1600

=== Normalized network (max_rowsum=0.2) ===
[[0.         0.06857143 0.08571429 0.04571429]
 [0.08       0.         0.06       0.03      ]
 [0.04761905 0.08571429 0.         0.06666667]
 [0.05       0.07       0.04       0.        ]]

Row sums after normalization:
  Patch 0: 0.2000
  Patch 1: 0.1700
  Patch 2: 0.2000
  Patch 3: 0.1600

Verification – no row sum exceeds 20%: True
```

**How `row_normalizer` works internally:**

```python
rowsums = network.sum(axis=1)
rows_to_renorm = rowsums > max_rowsum
network[rows_to_renorm] *= max_rowsum / rowsums[rows_to_renorm, np.newaxis]
```

Key points:
- **Import path**: `laser.core.migration` (the package namespace is `laser.core`, not `laser_core`)
- **`max_rowsum=0.20`** enforces the 20% cap; the argument must be in `[0, 1]`
- Rows whose sums already fall at or below the cap are **left unchanged** (patches 1 and 3 here)
- Rows that exceed the cap are **proportionally rescaled**, so relative weights between destinations are preserved (patches 0 and 2 here)
- The function returns a `float32` copy; the original matrix is not mutated
