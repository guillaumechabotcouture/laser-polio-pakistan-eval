Now I have the full source code. Here's the complete working example:

```python
import numpy as np
from laser.core.migration import gravity, row_normalizer

# --- 1. Define 4 patches ---
populations = np.array([50_000, 120_000, 30_000, 80_000])

# Pairwise distances (km) — symmetric, with 0 on diagonal
distances = np.array([
    [  0, 150, 300, 200],
    [150,   0, 250, 100],
    [300, 250,   0, 350],
    [200, 100, 350,   0],
], dtype=np.float64)

# --- 2. Build a gravity network ---
network = gravity(populations, distances, k=1e-4, a=1.0, b=1.0, c=2.0)

print("Raw gravity network:")
print(np.array2string(network, precision=4, suppress_small=True))
print("\nRow sums (fraction of FOI exported per patch):")
print(network.sum(axis=1).round(4))

# --- 3. Row-normalize so no patch exports more than 20% of its FOI ---
network_norm = row_normalizer(network, max_rowsum=0.20)

print("\nNormalized network (max 20% export per row):")
print(np.array2string(network_norm, precision=4, suppress_small=True))
print("\nRow sums after normalization:")
print(network_norm.sum(axis=1).round(4))
```

**How it works:**

- `gravity()` computes raw flows as `k * pop_i^a * pop_j^b / dist_ij^c`, with the diagonal zeroed out (no self-interaction).
- `row_normalizer(network, max_rowsum)` then checks each row's sum. Any row whose sum exceeds `max_rowsum` (here 0.20) gets scaled down proportionally so its sum equals exactly 0.20. Rows already at or below the threshold are left untouched. The key line in the implementation (`migration.py:137`):

  ```python
  network[rows_to_renorm] = network[rows_to_renorm] * max_rowsum / rowsums[rows_to_renorm, np.newaxis]
  ```

  This preserves the relative proportions between destinations within each row while capping the total outflow.

The `max_rowsum` parameter must be in `[0, 1]`, representing the fraction of force of infection a patch is allowed to export. Setting it to `0.20` ensures no patch sends more than 20% of its FOI to other patches.

Sources:
- [LASER GitHub Repository](https://github.com/InstituteforDiseaseModeling/laser)
- [LASER Migration Module Docs (v0.4.0)](https://docs.idmod.org/projects/laser/en/v0.4.0/migration.html)
