The complete model is written to `/tmp/eastland_seir.py`. Here's a walkthrough of every design decision:

---

## Key decisions explained

### Beta / R0
`BETA = R0 / mean_infectious_period = 5.0 / 10 = 0.5 day⁻¹`
Mean infectious period = shape × scale = 5 × 2 = 10 days.

### Seasonal forcing (piecewise, exact mean = 1.0)
Rather than a cosine (which only peaks at a single day), the array uses **flat plateaus** that literally sustain the specified values across the required day ranges:

| Segment | Days | Value | Sum |
|---------|------|-------|-----|
| Winter plateau | 0–90 (91 d) | 1.3 | 118.3 |
| Spring ramp | 91–149 (59 d) | 1.3→0.7 | 59.0 |
| Summer plateau | 150–240 (91 d) | 0.7 | 63.7 |
| Fall ramp | 241–364 (124 d) | 0.7→1.3 | 124.0 |
| **Total** | **365 d** | | **365.0 → mean = 1.0** |

Mean of `linspace(a, b, n)` is always `(a+b)/2 = 1.0`, so the mean is exactly 1.0.

### ValuesMap.from_timeseries
```python
seasonality = ValuesMap.from_timeseries(seasonal_ts, nnodes=4, nticks=3650)
```
The 365-day array is broadcast identically to all 4 patches and cycled modularly across the 3650-tick run (10 repetitions).

### Gravity network
Built manually after `Model.__init__` and assigned to `model.network`, overriding the framework's default. This gives full control over `k=0.01, a=1, b=1, c=1.5` and `row_normalizer(0.15)`.

### Component order
`[Susceptible, Exposed, Infectious, Recovered, Transmission, BirthsByCBR, MortalityByCDR]` — SEIR state management first, then new exposures via Transmission, then demographic flows.

### calc_capacity
Called explicitly with `birthrates.values` (shape `(3650, 4)`) before model construction; also passed as `birthrates=birthrates` to `Model()` so the framework can size the agent `LaserFrame` accordingly.
