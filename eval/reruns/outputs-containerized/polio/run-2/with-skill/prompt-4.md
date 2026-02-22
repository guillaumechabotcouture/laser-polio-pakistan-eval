`scripts/calibrate_polio.py` already implements all four requirements. Here's a precise mapping:

---

**Requirement 1 — AFP comparison with 1:200 ratio** → `score_afp_provincial()` (`calibrate_polio.py:172`)
```python
model_para = model_by_province[prov].mean() * PARALYSIS_RATIO  # × 1/200
loss += (np.log(model_para + 1) - np.log(obs_para + 1)) ** 2
```
Log-scale MSE, aggregated by province, compared against MMWR 2020–2025 averages.

---

**Requirement 2 — CCS zero-incidence metric** → `score_ccs_zero_weeks()` (`calibrate_polio.py:212`)

One design decision worth understanding: this doesn't count weeks where `infections == 0` (hard zeros). It computes the *expected* proportion of weeks with zero AFP detections given the 1:200 observation probability:
```python
p_zero_afp = (1.0 - PARALYSIS_RATIO) ** weekly   # (199/200)^W per week
model_prop_zero_afp = p_zero_afp.mean(axis=0)     # per district
```
This is the correct formulation when your observable is AFP, not infections — a week with 10 infections still has a 95% chance of zero AFP cases. Hard zeros (`weekly == 0`) would be a special case where `W=0`.

The observed targets (`obs_prop_zero_afp`) are derived from MMWR via a Poisson approximation in `compute_observed_targets()` (`calibrate_polio.py:106`).

---

**Requirement 3 — Parameter sampling** → `generate_samples()` (`calibrate_polio.py:447`)
```python
PARAM_RANGES = {
    "beta":                {"low": 0.15, "high": 0.30},
    "gravity_k_log10":     {"low": -3.0, "high": -1.0},   # 10^U(-3,-1)
    "seasonal_amplitude":  {"low": 1.0,  "high": 1.5},
}
```
Uses Latin Hypercube Sampling (`scipy.stats.qmc.LatinHypercube`) with a stratified-random fallback if `qmc` is unavailable. `gravity_k` is exponentiated after LHS so samples are log-uniform.

---

**Requirement 4 — Ranking** → `run_calibration()` (`calibrate_polio.py:488`)
```python
scores = score_combined(model, scenario, obs_prop_zero_afp,
                        afp_weight=0.5, ccs_weight=0.5)
...
results_df = results_df.sort_values("total_loss").reset_index(drop=True)
results_df["rank"] = range(1, len(results_df) + 1)
```
Equal-weighted combination of both loss components. Weights are CLI-adjustable via `--afp-weight` / `--ccs-weight`.

---

**What `polio_calibration.py` (655 lines) adds differently:** it uses hard Poisson approximation for the zero-week target (`exp(-weekly_rate)`) rather than the soft probabilistic form, and does not include diagnostic plots. The two files represent slightly different formulations of the same calibration problem.

**To run:**
```bash
python scripts/calibrate_polio.py --n-samples 32 --seed 42
# Outputs: eval/outputs/calibration_results.csv
#          eval/outputs/calibration_best_diagnostics.csv
#          eval/outputs/calibration_summary.png
```

The 4-panel summary plot (`plot_calibration_summary`, `calibrate_polio.py:701`) shows parameter space colored by loss, the CCS curve with logistic fit vs observed, provincial AFP bar comparison, and loss distribution — enough to diagnose whether the model is fitting the right regime.

If you want a standalone version decoupled from `polio_seir_10patch.py`, or a version wired to the calabaria `BaseModel` interface for Optuna-based optimization rather than a flat LHS sweep, let me know which direction is more useful.
