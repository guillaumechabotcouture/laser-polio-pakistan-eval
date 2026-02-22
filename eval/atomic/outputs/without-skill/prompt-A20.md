All constraints are satisfied:

- **Peak** = 1.5000 (day 0)
- **Trough** = 0.5000 (day 182)
- **Mean** = 1.000000

The core formula is:

```
seasonal_forcing = 1.0 + 0.5 * cos(2π * t / 365)
```

The cosine integrates to zero over a full period, so the mean is exactly the baseline (1.0). The amplitude of 0.5 swings the multiplier from 0.5 to 1.5. Day 0 represents the transmission peak (e.g., winter for respiratory diseases); you can shift the peak to any day by adding a phase offset like `cos(2π * (t - peak_day) / 365)`.
