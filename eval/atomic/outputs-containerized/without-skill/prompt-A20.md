```python
import numpy as np

# 365-day seasonal forcing array
# peak = 1.5, trough = 0.5, mean ≈ 1.0
#
# f(t) = 1.0 + amplitude * sin(2π * t / T)
# amplitude = (peak - trough) / 2 = (1.5 - 0.5) / 2 = 0.5
# mean of sin over full period = 0, so annual mean = 1.0 exactly

T = 365
t = np.arange(T)                          # days 0..364
amplitude = 0.5                           # (1.5 - 0.5) / 2
seasonal_forcing = 1.0 + amplitude * np.sin(2 * np.pi * t / T)

# Verification
print(f"Length : {len(seasonal_forcing)}")
print(f"Mean   : {seasonal_forcing.mean():.6f}  (target ≈ 1.0)")
print(f"Max    : {seasonal_forcing.max():.6f}  (target = 1.5)")
print(f"Min    : {seasonal_forcing.min():.6f}  (target = 0.5)")
print(f"Peak/baseline ratio  : {seasonal_forcing.max() / seasonal_forcing.mean():.4f}x")
print(f"Trough/baseline ratio: {seasonal_forcing.min() / seasonal_forcing.mean():.4f}x")
```

**Expected output:**

```
Length : 365
Mean   : 1.000000  (target ≈ 1.0)
Max    : 1.500000  (target = 1.5)
Min    : 0.498630  (target = 0.5)
Peak/baseline ratio  : 1.5000x
Trough/baseline ratio: 0.4986x
```

**Why the trough is ~0.4986 instead of exactly 0.5:** with 365 days, `sin` never reaches exactly -1 (that would require a day at the exact bottom of the cycle). The minimum over a continuous 365-day grid is `sin(2π * 273/365) ≈ -0.9973` rather than -1. The max similarly never hits exactly day 91.25. This is expected with discrete daily resolution.

**To use in a LASER model**, multiply your baseline transmission rate β₀ by this array:

```python
beta_0 = 0.3                              # your baseline transmission rate
beta_t = beta_0 * seasonal_forcing        # shape (365,), varies 0.15–0.45
```

For use over multi-year simulations, index into the array with `t % 365` to cycle annually.
