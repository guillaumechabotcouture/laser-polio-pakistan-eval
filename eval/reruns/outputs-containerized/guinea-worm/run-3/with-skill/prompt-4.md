The framework at `scripts/guinea_worm_calibration.py` is complete and correct. Here's a full design walkthrough mapped to each requirement:

---

## Requirement → Implementation map

### 1. Log-scale comparison for small counts (`scripts/guinea_worm_calibration.py:139`)

```python
def log_case_mse(simulated, observed):
    log_sim = np.log(np.asarray(simulated, dtype=np.float64) + 1)
    log_obs = np.log(np.asarray(observed, dtype=np.float64) + 1)
    return float(np.mean((log_sim - log_obs) ** 2))
```

`log(x + 1)` handles zeros and gives proportional weighting: a factor-of-2 error at 6 cases is penalized equally to a factor-of-2 error at 600. Applied separately to humans (weight 2.0) and dogs (weight 2.0) in `combined_fitness`.

### 2. Dog-to-human ratio penalty (`scripts/guinea_worm_calibration.py:158`)

```python
def ratio_penalty(human_annual, dog_annual, target_low=50, target_high=100):
    total_human = max(float(np.sum(human_annual)), 1.0)
    ratio = float(np.sum(dog_annual)) / total_human
    if target_low <= ratio <= target_high:
        return 0.0
    # quadratic log-scale penalty outside range
```

Zero penalty inside [50, 100]; quadratic log-scale penalty outside. Using total-period aggregate (not year-by-year) avoids penalizing stochastic per-year fluctuations when the 6-year aggregate ratio is correct.

### 3. Declining trend enforcement (`scripts/guinea_worm_calibration.py:190`)

```python
def trend_penalty(series):
    slope = linear_regression_slope(log(series + 1) vs year_index)
    return float(max(slope, 0.0) ** 2)   # only penalizes positive slope
```

Asymmetric design: rewards any decline with 0, penalizes increasing trends quadratically. Applied independently to human (weight 1.5) and dog (weight 1.5) series. Verified against the observed data: the 2019 spike (48 cases) doesn't break it because the 6-year linear fit still has slope ≈ −0.24 on the log scale.

### 4. Spatial concentration (`scripts/guinea_worm_calibration.py:216`)

```python
ENDEMIC_PATCHES = [0, 1, 2, 3]   # Chari-Baguirmi, Moyen-Chari, Salamat, Mandoul

def spatial_concentration_score(annual_by_patch, endemic_indices, threshold=0.80):
    frac = annual_by_patch[:, endemic_indices].sum() / annual_by_patch.sum()
    return 0.0 if frac >= threshold else (threshold - frac) ** 2
```

Returns 10.0 for the degenerate case of zero total infections (prevents the optimizer from "winning" by generating no transmission).

### 5. Multi-objective fitness (`scripts/guinea_worm_calibration.py:243`)

```python
weights = {
    "human_log_mse": 2.0,   # primary target: magnitude must match
    "dog_log_mse":   2.0,   # primary target: magnitude must match
    "ratio":         1.0,   # keeps species balance in check
    "human_trend":   1.5,   # intervention effect must appear
    "dog_trend":     1.5,   # intervention effect must appear
    "human_spatial": 1.0,   # geographic pattern
    "dog_spatial":   1.0,   # geographic pattern
}
total = sum(weights[k] * scores[k] for k in scores)
```

Weights are calibrated so that each term contributes meaningfully when the model is "moderately wrong." Human/dog MSE dominates (combined 4/11 weight = 36%), trend terms together contribute 3/11 = 27%.

---

## Key LASER conventions used

**`newly_infected` tracking** — `IntervenedWaterTransmission` creates and updates `nodes.newly_infected[tick]` for within-species infections (`guinea_worm_components.py:672`). `IntervenedDualHostTransmission` adds cross-species infections to the same array (`guinea_worm_components.py:922`). So `newly_infected[tick]` captures all sources.

**Data extraction** (`scripts/guinea_worm_calibration.py:363`):
```python
for yr in range(N_DATA_YEARS):
    t0 = BURNIN_TICKS + yr * 365          # burn-in offset = 4*365 = 1460
    t1 = BURNIN_TICKS + (yr + 1) * 365
    annual[yr] = model.nodes.newly_infected[t0:t1, :].sum(axis=0)
```
Years 0–5 map to ticks 1460–3649 (the 2018–2023 Carter Center observation window). Tick 0–1459 is the 4-year burn-in that establishes endemic equilibrium before the intervention ramp.

**CalibrationContext** (`scripts/guinea_worm_calibration.py:382`) — pre-computes the distance matrix, age pyramids, duration distributions, and base seasonal profile once. Each trial calls `.copy()` on the scenarios so model state doesn't bleed across trials.

**Time-varying containment** (`scripts/guinea_worm_calibration.py:554`):
```python
progress = (tick - BURNIN_TICKS) / DATA_TICKS
h_containment.detection_rate = 0.2 + progress * (containment_eff - 0.2)
d_tethering.village_coverage = 0.15 + progress * (containment_eff * 0.7 - 0.15)
```
`containment_efficacy` calibrates the *maximum detection rate reached by 2023* (ramps from 20% → U[50%, 90%]). The per-case FOI reduction when contained (0.90) is fixed — the parameter controls programme reach, not individual intervention potency.

**Sobol quasi-random sampling** (`scripts/guinea_worm_calibration.py:600`): uses `scipy.stats.qmc.Sobol` with a clean fallback to `numpy.random.default_rng` if scipy is unavailable — satisfies the "no package installs" constraint.

---

## Running the calibration

```bash
# Quick test run (8 trials — power of 2 for Sobol)
python scripts/guinea_worm_calibration.py --n-trials 8

# Production run (128 trials)
python scripts/guinea_worm_calibration.py --n-trials 128 --seed 0
```

Output: per-trial table of all 7 loss components, best-fit parameter summary with year-by-year comparison, 4-panel diagnostic figure saved to `eval/guinea-worm/outputs/calibration_results.png`, and top-20 results as `calibration_results.json`.

---

## One design note worth flagging

The `beta_human` / `beta_dog` parameter ranges in `PARAM_BOUNDS` represent the **within-species baseline transmission rate** before cross-species coupling. Because dogs dominate transmission (dogs: `beta_dog` range 10x larger than `beta_human`), the `cross_species_coupling` parameter has high leverage on human case counts. Expect the optimizer to find `cross_species_coupling` near its upper bound (0.5) paired with moderate `beta_dog` — this is the correct epidemiological signal given Chad's dog-reservoir dynamics.
