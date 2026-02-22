The file runs correctly with all assertions passing. The implementation in `seasonal_valuesmap_demo.py` already satisfies all requirements:

**Approach:** A cosine centered at day 225 (midpoint of the summer window 180–270) with amplitude 0.3. This naturally places:
- **Peak 1.30x at day 225** — falls within summer window [180, 270] ✓
- **Trough 0.70x at day 42** — falls within winter window [0, 90] ✓
- **Mean = 1.0000** — normalized, no R_eff bias ✓
- **ValuesMap shape (365, 4)** — same profile broadcast to all 4 patches ✓

The cosine is the correct choice here rather than a piecewise function: it satisfies the peak/trough constraints exactly while giving a smooth, physically realistic interpolation between summer and winter windows — no discontinuous jumps at the boundary days.

Key call: `ValuesMap.from_timeseries(season_365, npatches)` takes the 1D 365-day array and broadcasts it uniformly across all `npatches` patches. To pass it to the built-in transmission component:

```python
SEIR.Transmission(model, expdurdist, seasonality=seasonality)
```

If you need **per-patch profiles** (e.g., different seasonal timing by latitude), construct a `(365, 4)` array and use `ValuesMap.from_array(season_array)` instead.
