All tests pass. Here's a summary of the calibration framework in `scripts/guinea_worm_calibration.py`:

## Architecture

The framework builds on the existing dual-host SEIS model (`guinea_worm_chad.py` + `guinea_worm_components.py`) and adds calibration-specific scoring and optimization.

### Scoring Functions (5 components, 7 weighted terms)

| Component | Weight | What it measures |
|-----------|--------|-----------------|
| `log_case_mse` (human) | 2.0 | log(x+1) MSE of annual human cases vs Carter Center data |
| `log_case_mse` (dog) | 2.0 | log(x+1) MSE of annual dog infections vs Carter Center data |
| `ratio_penalty` | 1.0 | Penalizes dog:human ratio outside 50:1–100:1 |
| `trend_penalty` (human) | 1.5 | Penalizes non-declining trend (log-scale slope) |
| `trend_penalty` (dog) | 1.5 | Same for dogs |
| `spatial_concentration_score` (human) | 1.0 | Penalizes if < 80% of cases in 4 endemic patches |
| `spatial_concentration_score` (dog) | 1.0 | Same for dogs |

Exact match to observed data produces loss = 0.0 (verified by test).

### Model Runner (`run_trial`)

- **4-year burn-in** (no interventions, importation active) to reach quasi-equilibrium
- **6-year data period** (2018–2023) with containment ramping linearly:
  - Human case detection: 20% → `containment_efficacy`
  - Dog tethering: 15% → 70% of `containment_efficacy`
- Uses `IntervenedWaterTransmission` + `IntervenedDualHostTransmission` (properly track `newly_infected` from both within- and cross-species transmission)
- `CalibrationContext` pre-computes parameter-independent data (distance matrix, age pyramids, seasonal profile) once

### Calibration Loop

- **Sobol quasi-random sampling** via `scipy.stats.qmc.Sobol` (falls back to uniform if unavailable)
- Progress table with per-trial timing and all 7 fitness components
- Results sorted by loss; top-5 parameter sets printed
- Outputs: diagnostic plot (`calibration_results.png`) + JSON results

### Usage

```bash
/opt/anaconda3/bin/python3 scripts/guinea_worm_calibration.py --n-trials 64
```
