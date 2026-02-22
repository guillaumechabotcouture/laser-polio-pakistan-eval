The complete code:

```python
import numpy as np

amplitude = 0.5
seasonal_forcing = 1.0 + amplitude * np.cos(2 * np.pi * np.arange(365) / 365)
```

Verified properties:
- **Peak**: 1.5000 (day 0)
- **Trough**: 0.5000 (day 182)
- **Mean**: 1.0000000000

The cosine form `1.0 + 0.5 * cos(2πt/365)` guarantees the mean is exactly 1.0 analytically since the integral of cosine over a full period is zero. The tiny numerical residual (~10⁻¹⁶) comes from using 365 discrete samples rather than a continuous integral.

To shift the peak to a different day (e.g., day `peak_day`), replace `np.arange(365)` with `np.arange(365) - peak_day`.
