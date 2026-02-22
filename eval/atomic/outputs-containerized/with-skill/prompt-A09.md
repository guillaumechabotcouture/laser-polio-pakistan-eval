Here's the complete code. Note that `SEIR.Transmission` (`TransmissionSE`) takes an `expdurdist` — the exposed-period distribution assigned at the moment of S→E transition. The gamma-distributed **infectious** period you defined is `infdurdist`, which goes into `SEIR.Exposed` and `SEIR.Infectious`. Both are shown below in context.

```python
import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR
from laser.generic.utils import ValuesMap

# ── 1. Infectious-period distribution: Gamma(shape=5, scale=5.6) ──────────────
# Mean = shape × scale = 5 × 5.6 = 28 days
# Shape=5 gives moderate dispersion (CV = 1/√5 ≈ 0.45)
infdurdist = dists.gamma(shape=5, scale=5.6)

# Sanity-check the mean via sampling
_samples = np.array([infdurdist(0, 0) for _ in range(10_000)])
assert abs(_samples.mean() - 28.0) < 1.0, \
    f"infdurdist mean={_samples.mean():.2f}, expected ~28"

# ── 2. Exposed-period distribution (required by SEIR.Transmission) ────────────
# Separate from infdurdist — define to match your disease; example: ~10-day incubation
expdurdist = dists.gamma(shape=4, scale=2.5)   # mean = 10 days

# ── 3. Seasonal ValuesMap ─────────────────────────────────────────────────────
# Build a cosine profile (peak at day 15, amplitude 0.3)
# Replace with your disease-specific profile as needed.
nticks = 20 * 365   # example: 20-year simulation
nnodes = 10         # example: 10 spatial patches

days = np.arange(365)
peak_day = 15
amplitude = 0.3
season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
season_365 /= season_365.mean()   # Normalize — MUST be ~1.0 to avoid biasing R_eff

assert abs(season_365.mean() - 1.0) < 0.01, \
    f"Seasonal profile mean={season_365.mean():.4f}, must be ~1.0"

# Tile to cover the full simulation length
season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]

# Wrap as ValuesMap (same profile broadcast to all nodes)
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

# ── 4. Initialize SEIR.Transmission ──────────────────────────────────────────
# TransmissionSE(model, expdurdist, expdurmin=1, seasonality=None, validating=False)
# expdurdist: assigns the exposed timer at S→E; infdurdist goes to Exposed/Infectious
transmission = SEIR.Transmission(model, expdurdist, seasonality=seasonality)

# ── 5. Wire infdurdist into the components that use it ───────────────────────
# Both Exposed and InfectiousIR need infdurdist to schedule the I→R transition.
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),   # E→I assigns itimer via infdurdist
    SEIR.Infectious(model, infdurdist),             # I→R uses infdurdist
    SEIR.Recovered(model),
    transmission,                                   # S→E with seasonal forcing
]
```

**Key points:**

| Distribution | Parameter | Used by | Mean |
|---|---|---|---|
| `infdurdist` | `gamma(shape=5, scale=5.6)` | `SEIR.Exposed`, `SEIR.Infectious` | 28 days |
| `expdurdist` | `gamma(shape=4, scale=2.5)` | `SEIR.Transmission`, `SEIR.Exposed` | 10 days |

- `SEIR.Transmission` is `TransmissionSE` under the hood — its distribution argument (`expdurdist`) sets the exposed timer at the moment of new infection, not the infectious timer.
- The seasonal profile **must normalize to mean ≈ 1.0** (`season_365 /= season_365.mean()`); an un-normalized profile systematically biases effective R.
- `ValuesMap.from_timeseries` broadcasts one 1-D time series to all `nnodes` — use `ValuesMap.from_array(array_2d)` if you need per-node seasonal profiles.
