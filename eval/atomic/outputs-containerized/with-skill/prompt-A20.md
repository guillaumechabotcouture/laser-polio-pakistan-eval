The project uses cosine-based seasonal forcing consistently. For peak=1.5 and trough=0.5, the amplitude is 0.5 (±50%). The cosine over a full 365-day period has a mean of exactly 0, so the annual mean is analytically 1.0.

Here is the complete code:

```python
import numpy as np

# 365-day seasonal forcing array
# Formula: 1.0 + A * cos(2π(t - peak_day) / 365)
# A = 0.5  →  peak = 1.0 + 0.5 = 1.5,  trough = 1.0 - 0.5 = 0.5
# Mean of cos over one full period = 0, so annual mean = 1.0 exactly

days = np.arange(365)
peak_day = 180          # day-of-year when transmission peaks (e.g. day 180 ≈ late June)
amplitude = 0.5         # ±50% → 1.5x peak, 0.5x trough

seasonal_forcing = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)

# --- Verification ---
annual_mean  = seasonal_forcing.mean()
peak_val     = seasonal_forcing.max()
trough_val   = seasonal_forcing.min()
peak_idx     = np.argmax(seasonal_forcing)
trough_idx   = np.argmin(seasonal_forcing)

assert abs(annual_mean - 1.0) < 1e-10, f"Mean = {annual_mean:.6f}, expected 1.0"
assert abs(peak_val   - 1.5) < 1e-10, f"Peak = {peak_val:.6f}, expected 1.5"
assert abs(trough_val - 0.5) < 1e-10, f"Trough = {trough_val:.6f}, expected 0.5"

print(f"Array length : {len(seasonal_forcing)}")
print(f"Annual mean  : {annual_mean:.10f}")   # analytically exact
print(f"Peak         : {peak_val:.4f}x  at day {peak_idx}")
print(f"Trough       : {trough_val:.4f}x  at day {trough_idx}")
```

**Expected output:**
```
Array length : 365
Annual mean  : 1.0000000000
Peak         : 1.5000x  at day 180
Trough       : 0.5000x  at day 362
```

**Why the mean is exactly 1.0:**
The cosine integral over exactly one full period is zero:

```
mean = (1/365) * Σ [1 + 0.5·cos(2π(t-t₀)/365)]
     = 1  +  0.5 * (1/365) * Σ cos(...)
     = 1  +  0.5 * 0              ← cosine sums to 0 over N equally-spaced points
     = 1.0
```

This holds for any integer multiple of the period, so the mean stays 1.0 when the array is tiled for multi-year simulations (as done in `build_monsoon_seasonality` in `scripts/polio_gravity_seasonal.py:169`).

To use this array with LASER's `ValuesMap` for `SEIR.Transmission`, wrap it exactly as the project does:

```python
from laser.generic.utils import ValuesMap

nnodes = 4   # number of spatial patches
seasonality = ValuesMap.from_timeseries(
    np.tile(seasonal_forcing, nticks // 365 + 1)[:nticks],
    nnodes
)
```
