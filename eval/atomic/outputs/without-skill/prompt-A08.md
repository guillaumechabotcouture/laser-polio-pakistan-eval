The script runs successfully. Here's what it does:

**Seasonal profile** (`seasonal_transmission.py`):
- **365-day array** with piecewise construction and smooth cosine transitions between seasons
- **Winter trough (0.7x)**: days 0-90 and 330-364
- **Summer peak (1.3x)**: days 180-270
- **Smooth transitions**: cosine interpolation for days 91-179 (rising) and 271-329 (falling), passing through 1.0x at the midpoints

**ValuesMap creation**:
- Uses `ValuesMap.from_timeseries(data=seasonal, nnodes=4)` which broadcasts the single 365-day seasonal profile identically across all 4 patches
- Result shape is `(365, 4)` — 365 time steps by 4 patches
- The `.values` property or `[tick, patch]` indexing gives the transmission multiplier for any day and patch
