Good news: the complete calibration framework is already implemented in `scripts/guinea_worm_calibration.py` (890 lines). Here's how it maps to every requirement you specified.

---

## Framework Architecture

The file is structured in five sections:

```
Section 1: Scoring functions          (lines 139-286)
Section 2: Model infrastructure       (lines 293-375)
Section 3: Pre-computed context       (lines 382-593)
Section 4: Calibration loop           (lines 600-760)
Section 5: Diagnostics and plotting   (lines 767-869)
```

---

## Requirement Mapping

### 1. Log-scale annual case comparison

`log_case_mse()` at `scripts/guinea_worm_calibration.py:139`:

```python
def log_case_mse(simulated, observed):
    log_sim = np.log(np.asarray(simulated, dtype=np.float64) + 1)
    log_obs = np.log(np.asarray(observed, dtype=np.float64) + 1)
    return float(np.mean((log_sim - log_obs) ** 2))
```

`log(x+1)` guards against zeros and gives equal weight to a factor-of-2 error at 6 cases as at 1935. Both human and dog series get weight `2.0` in the combined loss.

---

### 2. Dog-to-human ratio (50:1–100:1)

`ratio_penalty()` at `scripts/guinea_worm_calibration.py:158`:

```python
def ratio_penalty(human_annual, dog_annual, target_low=50, target_high=100):
    total_human = max(float(np.sum(human_annual)), 1.0)
    ratio = float(np.sum(dog_annual)) / total_human
    if target_low <= ratio <= target_high:
        return 0.0          # no penalty inside range
    # quadratic log-scale penalty outside
    log_ratio = np.log(ratio + 1)
    ...
```

Zero penalty when the 6-year aggregate ratio is in [50, 100]; quadratic log-scale penalty outside. Weight `1.0` in combined loss.

---

### 3. Declining trend via linear regression

`trend_penalty()` at `scripts/guinea_worm_calibration.py:190`:

```python
def trend_penalty(series):
    x = np.arange(n, dtype=np.float64)
    y = np.log(np.asarray(series, dtype=np.float64) + 1)
    slope = np.sum((x - x_mean) * (y - y_mean)) / denom
    return float(max(slope, 0.0) ** 2)   # only penalizes positive slope
```

Fits OLS on log-counts. A declining series returns `0.0`. Applied separately to human and dog series, each at weight `1.5`. This asymmetry (penalize increasing, not declining) handles the 2019 spike in dogs (1,935) without penalizing a run that correctly captures it then falls.

---

### 4. Spatial concentration (3–4 of 8 districts)

`spatial_concentration_score()` at `scripts/guinea_worm_calibration.py:216`:

```python
ENDEMIC_PATCHES = [0, 1, 2, 3]   # Chari-Baguirmi, Moyen-Chari, Salamat, Mandoul

def spatial_concentration_score(annual_by_patch, endemic_indices, threshold=0.80):
    endemic_total = annual_by_patch[:, endemic_indices].sum()
    frac = endemic_total / total
    if frac >= threshold:
        return 0.0
    return float((threshold - frac) ** 2)
```

Requires ≥80% of all cases over the 6 data years to land in the 4 endemic districts. Returns a heavy penalty (`10.0`) if the model produces no cases at all, preventing trivial zero-infection solutions.

---

### 5. Multi-objective fitness

`combined_fitness()` at `scripts/guinea_worm_calibration.py:243`:

```python
weights = {
    "human_log_mse": 2.0,   # primary quantitative target
    "dog_log_mse":   2.0,   # primary quantitative target
    "ratio":         1.0,   # species balance
    "human_trend":   1.5,   # intervention effect, humans
    "dog_trend":     1.5,   # intervention effect, dogs
    "human_spatial": 1.0,   # geographic clustering
    "dog_spatial":   1.0,   # geographic clustering
}
total = sum(weights[k] * scores[k] for k in scores)
```

Total weight = 11.0. Human/dog MSE together dominate (36% of weight each), trend penalties enforce the intervention signal (27%), ratio and spatial act as soft constraints (18%).

---

### 6. Calibration loop with Sobol sampling

`calibrate()` at `scripts/guinea_worm_calibration.py:638`:

```python
ctx = CalibrationContext()            # builds scenario, distances, pyramids once
param_sets = sample_params(n_trials, method="sobol", seed=seed)

for i, params in enumerate(param_sets):
    outputs = run_trial(params, ctx, seed=seed + i)
    loss, scores = combined_fitness(
        outputs["human_annual"], outputs["dog_annual"],
        outputs["human_by_patch"], outputs["dog_by_patch"],
    )
    results.append({...})

results.sort(key=lambda r: r["loss"])
```

`CalibrationContext` (line 382) pre-builds the GeoDataFrames, haversine distance matrix, age pyramids, duration distributions, and base seasonal profile — all parameter-independent — once before the loop.

---

### 7. Intervention ramp in `run_trial()`

`run_trial()` at `scripts/guinea_worm_calibration.py:412`:

```python
BURNIN_YEARS = 4    # no interventions
DATA_YEARS   = 6    # 2018-2023, containment ramps up

for tick in range(TOTAL_TICKS):
    if tick < BURNIN_TICKS:
        h_containment.detection_rate = 0.0
        d_tethering.village_coverage = 0.0
    else:
        progress = (tick - BURNIN_TICKS) / DATA_TICKS
        h_containment.detection_rate = 0.2 + progress * (containment_eff - 0.2)
        d_tethering.village_coverage  = 0.15 + progress * (containment_eff * 0.7 - 0.15)
```

The linear ramp from detection=0.2 → `containment_efficacy` encodes the Carter Center's program scaling up over 2018-2023.

---

## Running the Calibration

```bash
# 64 trials (default, power-of-2 required for Sobol)
python scripts/guinea_worm_calibration.py

# 128 trials, different seed
python scripts/guinea_worm_calibration.py --n-trials 128 --seed 7
```

Output: per-trial table printed live, best-fit comparison table, top-5 parameter sets, diagnostic 4-panel plot saved to `eval/guinea-worm/outputs/calibration_results.png`, and top-20 results as `calibration_results.json`.

---

## Key Design Decisions

| Decision | Rationale |
|---|---|
| Sobol over random | Space-filling; 64 Sobol points cover the 5D box better than 64 uniform random |
| `log(x+1)` not `log(x)` | Chad has years with 6 human cases; log(0) is undefined |
| Trend only penalizes increases | 2019 is a genuine spike in dog data; the model should be able to reproduce it without penalty |
| Asymmetric cross-species coupling | `dog_to_human = coupling`, `human_to_dog = coupling × 0.33` — dogs are the dominant water contamination source |
| Threshold spatial penalty not MSE | District-level surveillance for guinea worm is sparse; matching presence/absence concentration is more reliable than exact counts per patch |
