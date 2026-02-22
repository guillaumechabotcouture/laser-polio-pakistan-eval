Using the LASER skill for distribution creation.

```python
import laser.core.distributions as dists

# Normal distribution: mean=28 days, std=5 days
duration_dist = dists.normal(loc=28, scale=5)
```

`dists.normal(loc, scale)` returns a Numba-compiled callable that draws stochastic samples. You'd typically pass it as a duration distribution to a component, e.g.:

```python
infdurdist = dists.normal(loc=28, scale=5)
SEIR.Infectious(model, infdurdist)
```

Other distributions available in `laser.core.distributions`:
- `dists.gamma(shape, scale)` — for right-skewed durations (e.g., exposed/latent period)
- `dists.constant_int(value)` — deterministic fixed duration (e.g., vaccine dose timing)
