Here is how to create a normal distribution with mean=28, std=5 using `laser.core.distributions`:

```python
import numpy as np
import laser.core.distributions as dists

# Create a Numba-wrapped normal distribution: mean=28 days, std=5 days
dur_dist = dists.normal(loc=28, scale=5)

# The result is a Numba JIT-compiled function with signature (tick, node) -> float32
# Sample a single value (tick and node are unused by stationary distributions)
single_sample = dur_dist(0, 0)
print(f"Single sample: {single_sample:.2f} days")

# Sample many values in parallel using sample_floats
n = 1000
dest = np.zeros(n, dtype=np.float32)
samples = dists.sample_floats(dur_dist, dest)
print(f"Mean: {samples.mean():.2f}, Std: {samples.std():.2f}")
```

**Key details from the source** (`laser/core/distributions.py`):

| Aspect | Detail |
|---|---|
| Import | `import laser.core.distributions as dists` |
| Factory function | `dists.normal(loc, scale)` — cached via `@lru_cache` |
| `loc` | Mean of the distribution (28) |
| `scale` | Standard deviation (5) |
| Return type | `@nb.njit` function with signature `(tick: int, node: int) -> float32` |
| Batch sampling | `dists.sample_floats(fn, dest)` — parallelized with `nb.prange` |

The `tick` and `node` arguments exist to support future spatially- or temporally-varying distributions. For a stationary normal distribution like this one, they are ignored — pass `0, 0` when calling directly.

In the existing codebase, this pattern appears at `scripts/polio_seir_10patch.py:191`:
```python
infdurdist = dists.normal(loc=PARAMS.inf_mean, scale=PARAMS.inf_sigma)
```
followed by `dists.sample_floats(infdurdist, np.zeros(n, dtype=np.float32))` to populate timer arrays for newly-infected agents.
