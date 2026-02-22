Here's how to create a normal distribution with mean 28 days and standard deviation 5 days using LASER's `dists` module:

```python
import laser.core.distributions as dists

# Create a normal distribution: mean=28 days, std=5 days
recovery_dist = dists.normal(loc=28.0, scale=5.0)

# Sample a single value
value = recovery_dist(0, 0)  # returns np.float32
```

**Key details:**

- `dists.normal(loc, scale)` returns a Numba-compiled (`@nb.njit`) callable with signature `(tick: int, node: int) -> np.float32`
- `loc` is the mean (28.0) and `scale` is the standard deviation (5.0)
- The `tick` and `node` arguments exist for interface consistency but are unused by the built-in distributions

**Bulk sampling** (for filling an array in parallel via Numba):

```python
import numpy as np
import laser.core.distributions as dists

recovery_dist = dists.normal(loc=28.0, scale=5.0)

# Pre-allocate and fill with samples (parallelized with nb.prange)
samples = np.zeros(10_000, dtype=np.float32)
dists.sample_floats(recovery_dist, samples)
```

The distribution functions are LRU-cached, so calling `dists.normal(28.0, 5.0)` multiple times returns the same compiled function object.
