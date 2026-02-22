Here's the complete code and explanation:

## `calc_capacity` — Units and Usage

`calc_capacity` expects birthrates in **per-1000/year** (e.g., 44 means "44 births per 1,000 people per year"). It uses these rates to project population growth over the simulation and pre-allocate enough agent slots in the `LaserFrame`. This is the same unit convention used by `BirthsByCBR`.

**The #1 silent failure**: if you pass daily per-capita rates (e.g., `0.00012` instead of `44`), `calc_capacity` sees near-zero growth, returns capacity ≈ initial population, and `LaserFrame.add()` silently has no free slots — births never happen.

## Complete Code

```python
import numpy as np
from laser.core.utils import calc_capacity
from laser.generic.utils import ValuesMap

# --- Parameters ---
initial_pop = np.array([500_000], dtype=np.int64)   # 1 patch, 500K people
cbr = 44                                              # per-1000/year (NOT daily per-capita)
nticks = 10 * 365                                     # 10-year simulation
nnodes = len(initial_pop)
safety_factor = 2.0                                   # extra headroom for stochastic variance

# --- Build birthrate array (nticks x nnodes), values in per-1000/year ---
birthrate_map = ValuesMap.from_scalar(cbr, nticks, nnodes)

# --- Validate units BEFORE calling calc_capacity ---
assert np.all(birthrate_map.values >= 1) and np.all(birthrate_map.values <= 60), \
    f"Birthrates must be per-1000/year, got {birthrate_map.values.min():.4f}"

# --- Pre-allocate agent storage ---
capacity = calc_capacity(birthrate_map.values, initial_pop, safety_factor=safety_factor)

print(f"Initial population:  {initial_pop[0]:>10,}")
print(f"Computed capacity:   {capacity[0]:>10,}")
print(f"Growth headroom:     {capacity[0] / initial_pop[0]:.2f}x")
```

## What to expect

With CBR = 44/1000/yr over 10 years, the population grows by roughly `500,000 × (1 + 44/1000)^10 ≈ 770,000`. With `safety_factor=2.0`, `calc_capacity` returns approximately **1,040,000–1,540,000** slots (depending on how it applies the factor to the projected growth). This ensures `LaserFrame.add()` always has free slots for newborns.

## Key parameters

| Parameter | Type | Expected Value |
|---|---|---|
| `birthrates` | ndarray `(nticks, nnodes)` | **per-1000/year**, range 10–50 typical |
| `initial_pop` | ndarray `(nnodes,)` | per-node population counts |
| `safety_factor` | float | 1.0 default; use 2–4 for high-growth or long runs |

The returned `capacity` is a `np.int32` array (one value per node) that gets passed to the `Model` constructor, which uses it to size the `LaserFrame` agent arrays. If growth exceeds capacity at runtime, `LaserFrame.add()` raises a `ValueError` — in that case, increase `safety_factor`.
