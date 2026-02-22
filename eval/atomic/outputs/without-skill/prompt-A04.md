Now I have the details. Here's how to use LASER's distribution utilities to create a gamma distribution:

## Creating a Gamma Distribution in LASER

```python
from laser.core import distributions as dist

# Gamma distribution: mean = shape * scale
# mean=10 days, shape=5 → scale = 10/5 = 2
infectious_period = dist.gamma(shape=5, scale=2)
```

The `dist.gamma(shape, scale)` function returns a **Numba-compiled callable** that you can use to draw samples. It accepts `tick` and `node` parameters (for time-step and spatial-node context in simulations):

```python
import numpy as np
from laser.core.distributions import gamma, sample_floats

# Create the distribution
infectious_period = gamma(shape=5, scale=2)

# Draw a single sample
duration = infectious_period(tick=0, node=0)  # e.g., 8.7 days

# Fill an array of durations for many agents at once
durations = np.empty(1000, dtype=np.float32)
sample_floats(infectious_period, durations, tick=0, node=0)
```

## What This Represents in an SEIR Model

In an **SEIR** (Susceptible → Exposed → Infectious → Recovered) model, this gamma distribution would model the **infectious period duration** — how long each agent stays in the **I** compartment before recovering.

**Why gamma instead of exponential?**

- An **exponential** distribution (shape=1) implies a constant hazard of recovery — agents are equally likely to recover at any moment, which is biologically unrealistic.
- A **gamma** with shape=5 concentrates durations around the mean of 10 days with moderate variance (variance = shape × scale² = 20, so SD ≈ 4.5 days). This better reflects reality: most people are infectious for roughly similar durations, with some natural variation.
- As shape increases, the distribution becomes more peaked. Shape=5 is a common choice — enough regularity to be realistic, enough spread to capture individual variation.

**Other SEIR uses for gamma distributions:**
- **Exposed (latent) period** — time in the E compartment before becoming infectious
- **Incubation period** — time from infection to symptom onset

The key advantage of LASER's implementation is that these distributions are **Numba JIT-compiled**, so sampling millions of agent durations in a spatial simulation remains fast.
