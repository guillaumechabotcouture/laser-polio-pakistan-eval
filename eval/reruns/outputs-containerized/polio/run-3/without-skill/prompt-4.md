All files are at `/home/evaluser/`. Here's a summary of every design decision:

---

## Calibration framework — file layout

```
calibration/
    __init__.py          clean public API surface
    parameters.py        ParameterBounds, CalibrationParams, CalibrationResult
    scoring.py           AFP conversion + ZWP fitness metric
    sampler.py           Latin Hypercube and Monte Carlo samplers
    runner.py            CalibrationRunner — loop, scoring, ranking
calibrate.py             entry point with LASER model stub
```

---

## Key design decisions

### AFP → infection conversion (`scoring.py:expected_afp`)
`sim_afp = sim_infections × (1/200)`. Kept as real-valued expected counts rather than Poisson-sampled integers: sampling would introduce extra variance that inflates apparent zero-week counts and obscures parameter signal.

### Zero-week proportion metric (`scoring.py:score_zero_week_proportion`)
```
ZWP[d] = fraction of weeks where cases[d] < 0.5
score   = −mean_d |obs_ZWP[d] − sim_ZWP[d]|
```
Range is `(−1, 0]`; zero is a perfect fit. Negated so that *higher = better* — consistent convention throughout the runner. The 0.5 threshold treats both integer zeros (observed AFP) and floating-point near-zeros (simulated) uniformly.

### Parameter sampling (`sampler.py`)
**Latin Hypercube** (default) divides each axis into `n` equal-probability strata and draws one sample per stratum with randomised column permutations. For 3 parameters at 100 draws, LHS gives coverage equivalent to ~1000 pure Monte Carlo draws. `gravity_k` is sampled as `10^U(−3,−1)` (log-uniform) so the prior is flat in log-space, matching its multi-order-of-magnitude uncertainty.

### Calibration loop (`runner.py:CalibrationRunner.run`)
- Each model run is wrapped in a `try/except` — one bad draw cannot abort the ensemble.
- Failed runs receive `score = −inf` and are excluded from ranking by default.
- Progress is logged every `log_every` runs (default 10) with the current best score.

### LASER model interface
The only contract imposed on the LASER model is:

```python
def run_laser_model(params: CalibrationParams) -> np.ndarray:
    # shape must be (n_weeks, n_districts)
```

The `calibrate.py` stub shows exactly where to slot in a real `PolioModel` instantiation, with parameter wiring comments for `transmission.beta`, `transmission.seasonal_amplitude`, and `mobility.gravity_k`.
