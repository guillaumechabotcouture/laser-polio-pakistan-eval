Here's how to create that gamma distribution and what it means in an SEIR context:

```python
import laser.core.distributions as dists
import numpy as np

# Gamma(shape=5, scale=2): mean = shape * scale = 5 * 2 = 10 days
inf_dur_dist = dists.gamma(shape=5, scale=2)
```

**Verification — confirm the parameterization:**

```python
# Sample a large batch to verify empirical mean and spread
offsets = np.zeros(100_000, dtype=np.int32)  # dummy offsets for batch sampling
samples = dists.sample_floats(inf_dur_dist, offsets)

print(f"Empirical mean:    {samples.mean():.2f} days  (expected: 10.0)")
print(f"Empirical std dev: {samples.std():.2f} days  (expected: {np.sqrt(5) * 2:.2f})")
print(f"5th–95th pct:      {np.percentile(samples, 5):.1f} – {np.percentile(samples, 95):.1f} days")
```

**How it's used in the SEIR model:**

In LASER's `SEIR.Infectious` component, each agent draws from this distribution when entering the **I** compartment. The sample sets that agent's individual infectious period — how many ticks until they transition to **R**.

```python
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, exp_dur_dist, inf_dur_dist),
    SEIR.Infectious(model, inf_dur_dist),   # <-- uses the gamma dist here
    SEIR.Recovered(model),
    ...
]
```

**Why gamma instead of exponential or constant:**

| Distribution | Shape | Epidemiological meaning |
|---|---|---|
| `dists.exponential(scale=10)` | Monotone-decreasing | Memoryless; most agents recover quickly, tail is long |
| `dists.gamma(shape=1, scale=10)` | Same as exponential | Identical |
| `dists.gamma(shape=5, scale=2)` | Unimodal, mean=10 | Most agents infect for 7–13 days; few very long/short |
| `dists.constant_int(10)` | Delta at 10 | All agents infect exactly 10 days |

The gamma with `shape=5` (coefficient of variation = `1/sqrt(5) ≈ 0.45`) is a standard epidemiological choice for infectious period: it avoids the unrealistic long tail of the exponential, and avoids the unrealistic rigidity of a constant, producing a biologically plausible spread of durations.

**Relationship to the SEIR `exp_shape`/`exp_scale` convention** you may see in existing parameters:

```python
# The skill's reference model uses this pattern for the exposed period:
parameters = PropertySet({
    "exp_shape": 40,    # High shape → very narrow (almost deterministic ~10 days)
    "exp_scale": 0.25,  # mean = 40 * 0.25 = 10 days
    ...
})
expdurdist = dists.gamma(shape=parameters.exp_shape, scale=parameters.exp_scale)

# Your infectious period with shape=5 has more variance — realistic for I compartment
infdurdist = dists.gamma(shape=5, scale=2)  # mean=10, std≈4.5 days
```

The exposed period often uses a high shape (narrow variance) because latency is biologically consistent across hosts, while the infectious period uses a lower shape to reflect variation in symptom severity and behavior.
