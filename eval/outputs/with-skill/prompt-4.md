The existing `scripts/calibrate_polio.py` already implements everything you asked for. Let me walk through the architecture and key design decisions.

## Existing Calibration Framework Overview

### Architecture

```
calibrate_polio.py
├── Parameter Space (PARAM_RANGES, FIXED)
├── Observed Targets (compute_observed_targets)
│   └── MMWR AFP data → per-district zero-week proportions via Poisson
├── Scoring Functions
│   ├── score_afp_provincial()   — log-MSE of model vs MMWR paralytic cases
│   ├── score_ccs_zero_weeks()   — zero-AFP-week proportion match
│   └── score_combined()         — weighted sum (50/50 default)
├── Model Runner
│   ├── precompute_distance_matrix()
│   └── build_and_run()          — constructs full SEIRV model per sample
├── Sampling
│   └── generate_samples()       — Latin Hypercube Sampling (scipy.stats.qmc)
├── Calibration Loop
│   └── run_calibration()        — iterate, score, rank, save CSV
└── Diagnostics
    └── plot_calibration_summary() — 4-panel figure
```

### Scoring Functions

**1. AFP Provincial Loss** (`score_afp_provincial`, line 172):
- Aggregates post-burn-in infections by province using `DISTRICT_TO_PROVINCE` mapping
- Converts to expected paralytic cases via `× 1/200`
- Computes `MSE(log(model+1), log(obs+1))` across 5 provinces
- Log transform prevents large-population provinces from dominating

**2. CCS Zero-Week Loss** (`score_ccs_zero_weeks`, line 212):
- For each post-burn-in week and district, computes `P(0 AFP) = (199/200)^W` where W = weekly infections
- Averages across weeks to get model proportion of zero-AFP weeks per district
- Compares to observed targets derived from MMWR annual data via Poisson: `P(0) = exp(-annual_cases/52)`
- MSE across all 13 districts

**3. Combined** (`score_combined`, line 242): Weighted sum, default 50/50.

### Sampling

LHS with 3 dimensions (beta, gravity_k in log-space, seasonal_amplitude). Falls back to stratified random if `scipy.stats.qmc` is unavailable. Each sample gets a unique seed (`seed + idx`).

### Key Design Choices

- **Memory management**: Explicit `del model, scenario; gc.collect()` between runs — the LASER model with 20M agents is ~GBs
- **Distance matrix pre-computation**: Haversine distances computed once, reused across all samples
- **Seasonal amplitude > 1**: At amplitude=1.5, the seasonal multiplier `1 + 1.5*cos(...)` ranges from -0.5 to 2.5, floored at 0 — effectively creating a "silent season" where transmission drops to zero. This is a deliberate modeling choice.
- **Separate PRNG seeds per sample**: `seed + idx` ensures reproducibility while giving each run independent stochasticity

### Potential Improvements

1. **Calabaria integration**: The current script is self-contained. Wrapping as a `BaseModel` would enable Optuna TPE optimization (adaptive sampling, not just LHS) and cloud scaling via `mops`.

2. **Parallelization**: Currently sequential. Each run takes ~60-120s, so 32 samples = ~1 hour. Could use `multiprocessing` for local speedup.

3. **Additional metrics**: Wavelet phase similarity for traveling wave validation, per-district annual case count distributions.

4. **Robustness**: Multiple seeds per parameter set to average over stochastic variation, rather than 1 seed per sample.

The framework is complete and functional as-is. Would you like me to make any specific improvements, or would you prefer the calabaria `BaseModel` redesign for more structured optimization?
