Now I have enough information to write accurate code. Let me put it together.

Based on the source code I've gathered from the LASER repository, here is the complete code:

```python
import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap

# ── 1. Gamma-distributed infectious period ─────────────────────────────────────
# shape=5, scale=5.6  →  mean = shape × scale = 28 days
# This callable has signature (tick, node) -> float, as required by InfectiousIR.
infdurdist = dists.gamma(shape=5, scale=5.6)

# ── 2. Incubation (exposure) duration distribution ────────────────────────────
# TransmissionSE requires an expdurdist for setting the E-timer on newly
# infected agents. Define it separately (adjust shape/scale for your disease).
expdurdist = dists.gamma(shape=4, scale=2.5)   # example: mean 10-day incubation

# ── 3. Seasonal ValuesMap ─────────────────────────────────────────────────────
# Build a sinusoidal annual seasonality multiplier for each tick.
# Amplitude 0.3 → transmission varies ±30 % around beta; peak at tick 0.
nticks     = model.params.nticks
nnodes     = model.nodes.count
amplitude  = 0.3
phase      = 0.0   # shift peak (radians); e.g. np.pi shifts peak to day 182

t = np.arange(nticks, dtype=np.float32)
seasonal_values = (1.0 + amplitude * np.sin(2.0 * np.pi * t / 365.0 + phase)).astype(np.float32)

# from_timeseries broadcasts one time-series across all nodes.
seasonality = ValuesMap.from_timeseries(seasonal_values, nnodes=nnodes)

# ── 4. Assemble SEIR components ───────────────────────────────────────────────
# Component roles:
#   Transmission (TransmissionSE) – S→E transition; uses expdurdist + seasonality
#   Exposed      – E-timer countdown; uses expdurdist & infdurdist
#   Infectious   (InfectiousIR)   – I→R transition; uses infdurdist (gamma, mean=28 d)
#   Susceptible / Recovered       – bookkeeping only

susceptible  = SEIR.Susceptible(model)
exposed      = SEIR.Exposed(model, expdurdist, infdurdist)
infectious   = SEIR.Infectious(model, infdurdist)        # ← gamma(5, 5.6) here
recovered    = SEIR.Recovered(model)
transmission = SEIR.Transmission(                        # ← seasonality here
    model,
    expdurdist,
    expdurmin=1,
    seasonality=seasonality,
)

model.components = [susceptible, exposed, infectious, recovered, transmission]
```

**A few notes on the design:**

- `dists.gamma(shape, scale)` returns a `(tick, node) -> float` callable that both `InfectiousIR` (`infdurdist`) and `TransmissionSE` (`expdurdist`) expect.
- The **infectious period** distribution (`infdurdist`) is wired into `SEIR.Infectious` and also into `SEIR.Exposed` (which samples the I-timer when the agent moves E→I). `SEIR.Transmission` (`TransmissionSE`) takes only the **incubation** duration distribution; it is the component where seasonality is applied.
- `ValuesMap.from_timeseries` accepts a 1-D float32 array of length `nticks` and broadcasts it across all nodes. Use `ValuesMap.from_array(data_2d)` if you need per-node seasonal variation (shape `nticks × nnodes`).
- The seasonal multiplier scales β at each tick: `ft = beta * seasonality[tick] * I/N`. Values must be non-negative (the `ValuesMap` constructor asserts this), so an amplitude above 1.0 would require a different functional form (e.g., `(1 + a*sin(…)) / (1 + a)` normalised).
