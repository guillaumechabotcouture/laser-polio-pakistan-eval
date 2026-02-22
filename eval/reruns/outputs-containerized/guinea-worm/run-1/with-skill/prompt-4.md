The calibration framework is already fully implemented at `scripts/guinea_worm_calibration.py` (committed in `a19f8cc`). Here's a detailed walkthrough of how it satisfies each requirement.

---

## Architecture

The file has five sections:

```
Section 1: Scoring Functions     — log_case_mse, ratio_penalty, trend_penalty,
                                   spatial_concentration_score, combined_fitness
Section 2: Model Infrastructure  — build_scenario, build_distance_matrix,
                                   build_gravity_network, build_age_pyramid,
                                   extract_data_period_annual
Section 3: CalibrationContext    — pre-computes shared data (distances, pyramids,
                                   duration distributions, base season)
           run_trial             — one Sobol sample → LASER dual-host run → scores
Section 4: calibrate()           — outer loop: sample → run → rank → print
Section 5: plot_calibration_results — 4-panel diagnostics + JSON export
```

---

## Objective 1 — Log-scale MSE for annual counts (weight 2.0 each)

`log_case_mse` at line 139:

```python
def log_case_mse(simulated, observed):
    log_sim = np.log(np.asarray(simulated, dtype=np.float64) + 1)
    log_obs = np.log(np.asarray(observed, dtype=np.float64) + 1)
    return float(np.mean((log_sim - log_obs) ** 2))
```

`log(x+1)` handles zeros and makes a factor-of-2 error at 10 cases cost the same as a factor-of-2 error at 1000 cases — appropriate for counts ranging from 6 to 1935.

---

## Objective 2 — Dog-to-human ratio penalty (weight 1.0)

`ratio_penalty` at line 158 operates on **aggregate totals** across all six data years rather than per-year ratios, which is more robust when annual human counts are as small as 6:

```python
ratio = total_dog / total_human     # target: 50–100
if target_low <= ratio <= target_high:
    return 0.0
# Quadratic log-scale penalty outside range:
return float((log_ratio - log_bound) ** 2)
```

---

## Objective 3 — Declining trend (weight 1.5 each species)

`trend_penalty` at line 190 fits a linear regression on log-counts over the six years, penalizing only *positive* slopes (increasing trend). A perfectly declining series scores 0.0:

```python
slope = np.sum((x - x_mean) * (y - y_mean)) / denom
return float(max(slope, 0.0) ** 2)
```

The trend is also driven mechanistically: `run_trial` at lines 553–565 ramps `h_containment.detection_rate` from 0.2 → `containment_efficacy` and `d_tethering.village_coverage` from 0.15 → `0.7×containment_efficacy` linearly over the 6 data years, mirroring the real Carter Center scale-up.

---

## Objective 4 — Spatial concentration (weight 1.0 each species)

`spatial_concentration_score` at line 216 uses the `(years × patches)` array directly:

```python
ENDEMIC_PATCHES = [0, 1, 2, 3]   # Chari-Baguirmi, Moyen-Chari, Salamat, Mandoul
frac = annual_by_patch[:, endemic_indices].sum() / annual_by_patch.sum()
if frac >= threshold:   # threshold=0.80
    return 0.0
return float((threshold - frac) ** 2)
```

Total absence (all zeros) triggers a hard penalty of 10.0 to prevent vacuous "spatial concentration".

---

## Composite fitness

`combined_fitness` at line 243:

```python
weights = {
    "human_log_mse": 2.0,   # primary magnitude match
    "dog_log_mse":   2.0,   # primary magnitude match
    "ratio":         1.0,   # species balance
    "human_trend":   1.5,   # intervention effect
    "dog_trend":     1.5,   # intervention effect
    "human_spatial": 1.0,   # geographic pattern
    "dog_spatial":   1.0,   # geographic pattern
}
total = sum(weights[k] * scores[k] for k in scores)
```

Log-MSE objectives carry 2× the weight of ratio/spatial because they're the primary quantitative targets. Trend is 1.5× — important but partly redundant with the mechanistic ramp.

---

## LASER conventions in `run_trial`

The trial runner follows the dual-host manual tick loop pattern from `guinea_worm_chad.py`:

```python
for tick in range(TOTAL_TICKS):
    human_model._initialize_flows(tick)   # propagate S/E/I/R from t to t+1
    dog_model._initialize_flows(tick)

    for c in human_model.components:
        c.step(tick)
    for c in dog_model.components:
        c.step(tick)

    cross_host.step(tick)   # IntervenedDualHostTransmission
    h_import.step(tick)
    d_import.step(tick)
```

Output extraction uses `model.nodes.newly_infected[t0:t1, :]` — the standard LASER incidence array written by `TransmissionSE`/`IntervenedWaterTransmission`.

`CalibrationContext` (line 382) pre-computes distance matrix, age pyramids, and duration distributions once and reuses them across all trials — the most expensive shared computation.

---

## Sampling

`sample_params` at line 600 uses `scipy.stats.qmc.Sobol` (power-of-2 samples for optimal space-filling) with automatic fallback to uniform random if scipy is unavailable — no additional install required.

---

## Usage

```bash
python scripts/guinea_worm_calibration.py              # 64 trials (Sobol)
python scripts/guinea_worm_calibration.py --n-trials 128
python scripts/guinea_worm_calibration.py --n-trials 32 --seed 7
```

Output: `eval/guinea-worm/outputs/calibration_results.png` (4-panel diagnostics) and `calibration_results.json` (top-20 trials with full scores and parameter values).
