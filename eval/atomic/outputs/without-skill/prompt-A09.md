Here's the complete code for setting up an SEIR model with a gamma-distributed infectious period (shape=5, scale=5.6, mean=28 days) and seasonal forcing using the LASER framework:

```python
import numpy as np
from laser.core import PropertySet
from laser.core.random import seed
from laser.core.distributions import gamma
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap, grid, get_default_parameters

# Reproducibility
seed(42)

# --- Distribution Setup ---
# Gamma-distributed infectious period: mean = shape * scale = 5 * 5.6 = 28 days
infectious_dist = gamma(shape=5.0, scale=5.6)

# Incubation (exposure) period distribution (e.g., mean = 10 days)
incubation_dist = gamma(shape=4.0, scale=2.5)

# --- Seasonal ValuesMap ---
# Sinusoidal seasonal forcing over 365 days
nticks = 365
nnodes = 9
days = np.arange(nticks)
amplitude = 0.2
phase = 0.0
seasonal_curve = (1.0 + amplitude * np.cos(2 * np.pi * (days - phase) / 365.0)).astype(np.float32)
seasonality = ValuesMap.from_timeseries(seasonal_curve, nnodes=nnodes, nticks=nticks)

# --- Scenario (3x3 grid of nodes) ---
scenario = grid(M=3, N=3, population_fn=lambda i, j: 10_000)
scenario["I"] = [10, 0, 0, 0, 0, 0, 0, 0, 0]
scenario["S"] = scenario["population"] - scenario["I"]
scenario["E"] = 0
scenario["R"] = 0

# --- Parameters ---
params = get_default_parameters() | PropertySet({
    "nticks": nticks,
    "beta": 0.8,
    "seed": 42,
})

# --- Build and Run Model ---
model = Model(scenario, params)

model.components = [
    SEIR.Susceptible(model),
    SEIR.Transmission(model, expdurdist=incubation_dist, expdurmin=1, seasonality=seasonality),
    SEIR.Exposed(model, expdurdist=incubation_dist, infdurdist=infectious_dist, expdurmin=1, infdurmin=1),
    SEIR.Infectious(model, infdurdist=infectious_dist, infdurmin=1),
    SEIR.Recovered(model),
]

model.run()
```

**Key details:**

- **`gamma(shape=5.0, scale=5.6)`** — Returns a Numba-compiled callable with signature `(tick, node) -> float32`. The mean infectious period is `5 * 5.6 = 28 days`. The shape parameter controls variance; higher shape means more concentrated around the mean.

- **`SEIR.Transmission`** — This is `TransmissionSE` (S→E transition). It takes `expdurdist` (the *incubation* distribution, not the infectious one) because Transmission assigns the exposure timer when agents move from S to E. The infectious distribution is passed to `SEIR.Exposed` and `SEIR.Infectious`.

- **`ValuesMap.from_timeseries(curve, nnodes, nticks)`** — Applies the same seasonal curve to all nodes. The Transmission component uses it as `beta * seasonality[tick] * I/N` each timestep. Other factory methods:
  - `ValuesMap.from_scalar(1.0, ...)` — No seasonality (flat)
  - `ValuesMap.from_nodes(array, nticks)` — Per-node static multipliers
  - `ValuesMap.from_array(2d_array)` — Full `(nticks, nnodes)` custom matrix
