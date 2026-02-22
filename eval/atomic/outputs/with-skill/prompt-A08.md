All assertions pass. Here's a summary:

- **Profile**: Cosine with amplitude 0.3, centered at day 225 (midpoint of summer window)
- **Peak**: 1.30x at day 225 (within summer days 180-270)
- **Trough**: 0.70x at day 42 (within winter days 0-90)
- **Mean**: exactly 1.0 (no normalization bias on R_eff)
- **ValuesMap**: `from_timeseries` broadcasts the 365-day profile identically across all 4 patches

The key call is `ValuesMap.from_timeseries(season_365, npatches)` — this takes a 1D time series and replicates it for each patch. To use it in a model, pass it as the `seasonality=` parameter to `SEIR.Transmission()`.
