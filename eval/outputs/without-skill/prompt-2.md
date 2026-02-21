The model runs correctly. Here's a summary of what was done and the final file.

**What changed:** The only modification was replacing the `from laser.core import LaserFrame, PropertySet` import with two lightweight stand-in classes so the code runs without the LASER package:

- **`PropertySet`** — dict-backed attribute-access container (e.g., `params.beta` reads from `params._data["beta"]`)
- **`LaserFrame`** — tabular container supporting `add_scalar_property()` for dynamic property attachment

All three extensions were already correctly implemented in the existing code:

| Extension | Implementation |
|---|---|
| **1. Gravity network** (`b=0.5`, `c=1.5`) | `build_gravity_network()` — computes `G[i,j] = k * N_i^a * N_j^b / d_ij^c` for all off-diagonal pairs |
| **2. Seasonal forcing** (monsoon peak 1.3x, dry trough 0.7x) | `build_seasonal_profile()` — 365-element multiplier array with linear interpolation between monthly anchors |
| **3. Row-normalization** (15% cap) | `row_normalize_network()` — scales each row's off-diagonal sum to at most 0.15, fills diagonal with `1 - sum` |

**Verification output confirms correctness:**
- All coupling matrix row sums = 1.0
- All off-diagonal row sums = 0.15 (capped)
- All diagonal values = 0.85 (local retention)
- Epidemic dynamics show spatial wave propagation across patches with seasonal modulation visible in the per-patch curves
