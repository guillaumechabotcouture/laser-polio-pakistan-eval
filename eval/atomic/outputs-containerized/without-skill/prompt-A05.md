Now I have a clear picture of the actual API. The module is `laser.core.distributions` (LASER's distribution utilities). Here's how to create a normal distribution with mean=28 days and standard deviation=5 days:

```python
import laser.core.distributions as dists
import numpy as np

# Create a Numba-wrapped normal distribution (mean=28 days, std=5 days)
duration_dist = dists.normal(loc=28.0, scale=5.0)

# Sample a single value (tick and node params are for temporal/spatial variation)
sample = duration_dist(tick=0, node=0)
print(f"Single sample: {sample:.2f} days")

# Bulk sampling into a pre-allocated float32 array
n_agents = 1000
samples = np.zeros(n_agents, dtype=np.float32)
dists.sample_floats(duration_dist, samples, tick=0, node=0)

print(f"Sample mean:  {samples.mean():.2f} days  (expected ~28)")
print(f"Sample std:   {samples.std():.2f} days   (expected ~5)")
```

**Key points:**

- The module is `laser.core.distributions`, which the LASER docs refer to as the distributions/dists utilities.
- `dists.normal(loc, scale)` uses `@lru_cache`, so calling it again with the same parameters returns the cached compiled function — no redundant JIT compilation.
- The returned function is a Numba `@njit`-compiled callable with signature `(tick: int, node: int) -> float32`. The `tick`/`node` arguments allow for future temporal or spatial parameterization but are unused in the standard distributions.
- `dists.sample_floats(fn, dest, tick, node)` uses `parallel=True` Numba for efficient bulk sampling into a pre-allocated `float32` array.
- For integer-valued distributions (e.g., `poisson`, `binomial`), use `dists.sample_ints` instead.

In a model component you'd typically assign the distribution once at initialization and call it per-agent during simulation:

```python
class InfectiousPeriod:
    def __init__(self):
        # Cached at construction; free to call repeatedly
        self.duration_dist = dists.normal(loc=28.0, scale=5.0)

    def draw(self, tick, node):
        return self.duration_dist(tick, node)
```
