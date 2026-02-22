All files verified and working. Here is a full account of what was built:

---

## Calibration framework: file layout

```
gworm/
├── calibration/
│   ├── observed_data.py   – surveillance data, populations, at-risk fractions
│   ├── scoring.py         – five-objective fitness functions
│   └── calibrate.py       – LHS sampler, refinement loop, CLI entry point
└── model/
    └── guinea_worm.py     – Water–Host compartmental forward model
```

---

## Key design decisions

### Forward model (`model/guinea_worm.py`)

The annual time-step Water–Host model separates **at-risk** subpopulations from full district totals, which is what makes the beta priors meaningful:

| At-risk fraction | Reasoning |
|---|---|
| **0.5 % of humans** | Only those using unprotected traditional water sources |
| **20 % of dogs** | Free-roaming animals drinking from open water bodies |

Water contamination index `W[d] ∈ [0,1]` is normalized to the highest-burden district and updated from per-capita shedding after each step. Containment ramps linearly from 0 in 2018 to `containment_efficacy` in 2023, encoding GWEP scale-up. `PropertySet` follows the LASER dot-access convention throughout.

### Scoring (`calibration/scoring.py`)

| Objective | Score function | Weight |
|---|---|---|
| Human log-scale match | RMSE of log(obs+1) vs log(sim+1) | **2.0** |
| Dog log-scale match | same | **1.5** |
| Dog-to-human ratio | squared deviation outside \[50,100\] | **1.0** |
| Human declining trend | 1 − Pearson ρ on log-scale + slope-sign penalty | **1.5** |
| Dog declining trend | same | **1.0** |
| Spatial concentration | squared shortfall below 70 % in top-4 districts | **1.0** |

### Calibration loop (`calibration/calibrate.py`)

**Phase 1 — LHS global sweep** (default 5 000 samples, 1 000 × N_params):  
`latin_hypercube_sample` partitions each dimension into `n` strata and permutes independently, giving uniform marginal coverage with no clustering.

**Phase 2 — local refinement** (default 1 000 samples):  
`perturb_around` applies Gaussian noise (σ = 8 % of each prior range, in unit space) to the top-50 Phase-1 candidates, focusing evaluations on the promising region.

### Full-budget calibration result (5 000 + 1 000 samples, seed=42)

```
fitness = 3.106
2018: sim 18 human / 1092 dog  →  ratio 61:1  (obs 17 / 1040)
2023: sim  5 human /  289 dog  →  ratio 60:1  (obs  6 / 526)
spatial score = 0.0 (≥70 % concentration in top-4 districts)
ratio  score  = 0.0 (all years within 50–100:1 band)
```

The 2019 spike (48 human cases) is the hardest year to fit — the linear-ramp containment model predicts a monotone decline, not the 2019 rebound. The `score_trend` component handles this gracefully via log-scale Pearson correlation rather than penalising each year independently.

### CLI usage

```bash
python -m calibration.calibrate                          # defaults: 5k LHS + 1k refine
python -m calibration.calibrate --n-lhs 10000 --n-refine 2000 --top-k 100
python -m calibration.calibrate --seed 123 --save results.npz --quiet
```

**Prior boundary note:** `beta_human` and `beta_dog` consistently reach the upper end of their priors. If runs with larger `n-lhs` confirm this, extending both upper bounds by ~2× (e.g., `beta_human` → `[0.0005, 0.01]`) will improve the dog trend match.
