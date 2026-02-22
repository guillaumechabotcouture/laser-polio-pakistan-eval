Everything works. Here is the complete code with explanation:

---

```python
import numpy as np
from laser.generic.utils import ValuesMap

NTICKS = 365
NNODES = 4
TROUGH = 0.7
PEAK   = 1.3

# Season boundary days (inclusive, 0-indexed)
WINTER1_END   = 90    # days   0– 90 → winter trough (0.7x)
SUMMER_START  = 180   # days  91–179 → spring ramp (0.7→1.3)
SUMMER_END    = 270   # days 180–270 → summer peak  (1.3x)
WINTER2_START = 330   # days 271–329 → autumn ramp (1.3→0.7), 330–364 → winter trough

days = np.arange(NTICKS, dtype=np.float32)

seasonal = np.piecewise(
    days,
    [
        days <= WINTER1_END,
        (days > WINTER1_END) & (days < SUMMER_START),
        (days >= SUMMER_START) & (days <= SUMMER_END),
        (days > SUMMER_END)   & (days < WINTER2_START),
        days >= WINTER2_START,
    ],
    [
        TROUGH,
        lambda d: TROUGH + (PEAK - TROUGH) * (d - WINTER1_END) / (SUMMER_START - WINTER1_END),
        PEAK,
        lambda d: PEAK   + (TROUGH - PEAK) * (d - SUMMER_END)  / (WINTER2_START - SUMMER_END),
        TROUGH,
    ],
).astype(np.float32)

# Broadcast the 1-D time series across all 4 patches
seasonality_map = ValuesMap.from_timeseries(seasonal, nnodes=NNODES)
```

**Key design decisions:**

| Segment | Days | Value |
|---|---|---|
| Winter trough (start) | 0–90 | 0.7x constant |
| Spring ramp | 91–179 | linear 0.7 → 1.3 |
| Summer peak | 180–270 | 1.3x constant |
| Autumn ramp | 271–329 | linear 1.3 → 0.7 |
| Winter trough (end) | 330–364 | 0.7x constant |

- **`ValuesMap.from_timeseries`** is the right factory here — it uses `np.broadcast_to` internally, so the `(365, 4)` shaped `.values` array costs no extra memory per patch; all 4 patches alias the same underlying row data.
- **`np.piecewise`** with lambda ramps gives exact boundary values at days 90, 180, 270, and 330 without any floating-point overshoot.
- At simulation tick `t`, `seasonality_map[t]` returns a `(4,)` array of the same multiplier for each patch, ready to be used by a `TransmissionSI` or similar component.
