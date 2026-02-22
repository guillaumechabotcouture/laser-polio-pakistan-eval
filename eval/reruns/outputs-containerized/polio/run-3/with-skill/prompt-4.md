The framework is fully implemented. Here's the complete map of your four requirements to the existing code:

---

**`scripts/polio_calibration.py`** — already satisfies all four requirements:

**1. AFP → infection conversion (1:200 ratio)**

`derive_observed_targets()` (lines 87–127):
```python
# AFP/year ÷ 0.005 ÷ 52 = infections/week per district
district_infections_wk = district_afp_yr / PARALYSIS_RATIO / 52.0
# P(zero) under Poisson
expected_prop_zero = np.exp(-district_infections_wk)
```
`score_afp_match()` inverts this: `model_infections × 0.005 → expected AFP`, then computes province-level RMSLE against `OBSERVED_WPV1` (2020–2025 MMWR data).

**2. Zero-incidence-weeks CCS metric**

`compute_zero_week_proportions()` (line 139–147) + `score_ccs_pattern()` (lines 150–186):
- Direct SSE between model and Poisson-derived target prop_zero per district
- Logistic CCS fit over log₁₀(population), with penalty if the fitted midpoint x₀ falls outside [4.5, 6.5] (outside plausible polio CCS range)
- Combined: `direct_SSE + 0.3 × logistic_score`

**3. Parameter sampling**

`sample_parameters()` (lines 416–448) with Latin Hypercube Sampling (better space-filling than pure random; stratified fallback if scipy < 1.7):
```python
beta          = 0.15 + u[0] * 0.15          # U(0.15, 0.30)
gravity_k     = 10 ** (-3.0 + u[1] * 2.0)  # 10^U(-3, -1)
seasonal_amp  = 1.0  + u[2] * 0.5           # U(1.0, 1.5)
```

**4. Ranking**

`run_calibration()` (lines 455–562): returns a DataFrame sorted ascending by `total_score = 0.6×CCS + 0.4×AFP_RMSLE`, with a `rank` column inserted at position 0.

---

**Usage:**

```bash
python scripts/polio_calibration.py --n-samples 32 --seed 42
# Results written to eval/outputs/calibration_results.csv
```

---

**One thing to be aware of:** `seasonal_amplitude` in the sampler maps to cosine amplitude as `amplitude = seasonal_amplitude - 1.0`, so the range [1.0, 1.5] gives cosine amplitudes [0.0, 0.5] — a flat profile at the lower bound and ±50% peak-to-mean modulation at the upper bound. The FIXED `peak_day = 245` (early September, monsoon peak) is baked in and not calibrated. If you want to calibrate peak timing, add it to `PARAM_RANGES` and pass it through `run_single_simulation`.
