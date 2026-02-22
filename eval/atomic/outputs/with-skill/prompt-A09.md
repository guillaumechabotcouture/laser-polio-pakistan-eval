Here's the complete code for creating the gamma-distributed infectious period and initializing the `SEIR.Transmission` component with seasonal forcing:

```python
import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR
from laser.generic.utils import ValuesMap

# --- 1. Exposed-period distribution (gamma, shape=5, scale=5.6, mean=28 days) ---
expdurdist = dists.gamma(shape=5, scale=5.6)  # mean = 5 * 5.6 = 28 days

# --- 2. Seasonal forcing profile (365-day cosine, peaked at day 0 / Jan 1) ---
days = np.arange(365)
peak_day = 0          # peak transmission day-of-year
amplitude = 0.3       # ±30% seasonal swing
season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
season_365 /= season_365.mean()  # normalize to mean == 1.0

assert abs(season_365.mean() - 1.0) < 0.01, \
    f"Seasonal profile mean={season_365.mean():.3f}, must be ~1.0"

# --- 3. Tile to simulation length and wrap in ValuesMap ---
nticks = 20 * 365     # example: 20-year simulation
nnodes = 13           # example: number of spatial patches

season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

# --- 4. Initialize the Transmission component ---
# SEIR.Transmission is an alias for TransmissionSE (S→E transition)
# It uses model.params.beta for the base transmission rate and
# model.network for spatial coupling (force of infection with migration).
transmission = SEIR.Transmission(model, expdurdist, seasonality=seasonality)
```

Key points:

- **`dists.gamma(shape=5, scale=5.6)`** — the `shape` and `scale` parameters follow numpy/scipy convention. Mean = shape × scale = 28 days. This is the *exposed* (latent) duration distribution passed to `TransmissionSE`, which draws from it to set the E→I timer when an agent becomes exposed.
- **`ValuesMap.from_timeseries()`** — takes a 1D array of length `nticks` and broadcasts it identically to all `nnodes` patches. Use `ValuesMap.from_array()` instead if you need per-patch seasonal profiles.
- **`SEIR.Transmission`** multiplies `model.params.beta` by the seasonality value at each tick to get the effective transmission rate. A mean of 1.0 in the profile ensures the annual-average R_eff matches your intended `beta`.
- The `expdurdist` argument tells the transmission component how long newly exposed agents stay in E before becoming infectious — it is **not** the infectious period. The infectious period distribution is passed separately to `SEIR.Exposed` and `SEIR.Infectious`.

If the gamma distribution is meant for the *infectious* period (I→R) rather than the exposed period, you'd use it like this instead:

```python
infdurdist = dists.gamma(shape=5, scale=5.6)  # infectious period, mean 28 days

# Pass to the Exposed and Infectious components (not Transmission):
SEIR.Exposed(model, expdurdist, infdurdist)
SEIR.Infectious(model, infdurdist)

# Transmission still needs the exposed-period distribution:
SEIR.Transmission(model, expdurdist, seasonality=seasonality)
```
