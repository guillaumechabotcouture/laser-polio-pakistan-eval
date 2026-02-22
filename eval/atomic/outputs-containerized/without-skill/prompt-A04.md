Here's how to create and use a gamma distribution with LASER's `laser.core.distributions` module:

## Creating the Gamma Distribution

```python
import laser.core.distributions as dists
import numpy as np

# Gamma(shape=5, scale=2) => mean = shape * scale = 5 * 2 = 10 days
#                             variance = shape * scale^2 = 5 * 4 = 20
infectious_period_dist = dists.gamma(shape=5, scale=2)
```

`dists.gamma()` is decorated with `@lru_cache`, so calling it again with the same arguments returns the same compiled object — no redundant compilation.

## Sampling from It

```python
# Single sample (returns np.float32)
duration = infectious_period_dist(tick=0, node=0)

# Bulk sampling across a population (Numba-parallelized)
n_agents = 100_000
buf = np.zeros(n_agents, dtype=np.float32)
durations = dists.sample_floats(infectious_period_dist, buf)
```

The returned closure always takes `(tick: int, node: int)` — even if unused — so the same interface works for future spatially or temporally varying distributions.

## What It Represents in an SEIR Model

| SEIR Parameter | Distribution Choice | Rationale |
|---|---|---|
| Incubation period (E→I) | Gamma or Erlang | Delay has a realistic lower bound |
| **Infectious period (I→R)** | **Gamma(5, 2) = mean 10d** | Matches observed clinical data; not exponential |
| Symptom onset delay | Lognormal | Right-skewed, never negative |

**Why gamma over exponential for infectious period:**

An exponential distribution (memoryless) implies an agent who has been infectious for 9 days is just as likely to recover tomorrow as one who just became infectious. Gamma(shape=5, scale=2) concentrates mass around 10 days with a coefficient of variation of `1/√5 ≈ 0.45`, giving a more realistic bell-shaped duration. This is sometimes called an **Erlang** distribution when the shape is an integer (which shape=5 is).

## Usage Pattern in a Recovery Component

```python
import numba as nb
import numpy as np
import laser.core.distributions as dists

# Parameterize once at model setup
inf_dist = dists.gamma(shape=5, scale=2)  # mean=10 days

# Inside a Numba-compiled recovery component, sample when agent becomes infectious
@nb.njit(nogil=True, parallel=True, cache=True)
def update_infectious_timers(count, timers, state, shape, scale):
    for i in nb.prange(count):
        if timers[i] > 0:
            timers[i] -= 1
            if timers[i] == 0:
                state[i] = 3  # I -> R
        elif state[i] == 2:  # newly infectious
            # Draw infectious duration from gamma
            timers[i] = max(1, round(np.random.gamma(shape, scale)))

# Call each tick
update_infectious_timers(pop_count, inf_timers, epi_state,
                         shape=5, scale=2)
```

Note: Within `@nb.njit` functions, you call `np.random.gamma(shape, scale)` directly rather than the `dists` closure — the closure pattern is used at the model configuration layer to allow swapping distributions without changing component code.

**Key parameters summary:**
- `shape=5` (k): controls how "peaky" the distribution is — higher shape → tighter, more symmetric
- `scale=2` (θ): the scale in days — mean = k·θ = 10, variance = k·θ² = 20, std ≈ 4.5 days
