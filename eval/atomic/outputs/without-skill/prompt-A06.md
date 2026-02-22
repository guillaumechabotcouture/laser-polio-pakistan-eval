The script runs successfully. Here's a summary of what it does:

**Distance matrix** — 4 patches at positions 0, 100, 200, 300 km, so distances are simply `|i - j| * 100` km.

**Gravity formula** — `network[i,j] = k * (pop_i^a * pop_j^b) / dist[i,j]^c`, matching the LASER `laser.core.migration.gravity()` API. The diagonal is set to 0 (no self-migration).

**Key results:**
- The strongest link is Patch 2 → Patch 3 (15,000) — large destination population (200k) at short distance (100 km)
- The weakest link is Patch 2 → Patch 0 / Patch 0 → Patch 2 (937.5) — moderate populations separated by 200 km with `c=2` distance decay
- The proportional outflow shows Patch 2 sends 20% of its population toward Patch 3, while Patch 3 sends only 7.5% toward Patch 2 (asymmetry from differing origin populations)

The LASER framework equivalent is a single call:
```python
from laser.core.migration import gravity
network = gravity(populations, distances, k=0.01, a=1, b=1, c=2)
```

The code is in `gravity_network.py`.
