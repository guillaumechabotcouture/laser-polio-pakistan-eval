All 26 tests pass. Here is a summary of everything produced.

---

## Calibration framework — 4 files, 26 tests passing

```
calibration/
├── __init__.py       clean public API
├── sampler.py        Latin Hypercube Sampler + prior definitions
├── scoring.py        fitness metrics
├── runner.py         LASER model adapter (integration stub)
├── calibrate.py      main sweep loop + CLI + ranking
└── tests.py          26 self-contained tests (numpy only)
```

### Design decisions

**`sampler.py` — LHS over three priors**

| Parameter | Distribution | Notes |
|---|---|---|
| `beta` | U(0.15, 0.30) | linear |
| `gravity_k` | log-uniform 10^U(−3,−1) | sampled in log space, back-transformed |
| `seasonal_amplitude` | U(1.0, 1.5) | linear |

LHS divides each axis into `n_samples` equal strata and draws exactly one point per stratum, then shuffles axes independently. This gives far better prior coverage than pure random Monte Carlo at the same budget.

**`scoring.py` — zero-incidence week proportion**

The primary metric (`score_zero_proportion`) works as follows:

1. Convert `sim_infections → sim_afp` by dividing by `PARALYSIS_RATIO = 200`.
2. Mask simulated weeks wherever the observed AFP record is missing (`NaN`), so the comparison is fair.
3. For each district, compute `p_zero = fraction of weeks with 0 cases`.
4. Districts with fewer than `MIN_OBSERVED_WEEKS = 26` valid records are excluded (NaN).
5. **Score = mean |p_zero_sim − p_zero_obs|** across qualifying districts — lower is better.

Two supplementary metrics are also computed (not used for ranking): relative total burden error and district-rank Spearman correlation.

**`calibrate.py` — sweep loop**

```python
from calibration import run_calibration, save_results

rows = run_calibration(
    n_samples=200,
    obs_afp=afp_array,       # (n_districts, n_weeks), NaN = missing
    base_params=fixed_params, # geography, demography, vaccine schedule
    n_weeks=260,              # 5 years weekly
    seed=42,
)
save_results(rows, "results/calib_sweep.npy")
# rows[0] is the best-fitting run
```

Results are stored as a numpy structured array (`RESULT_DTYPE`) with all parameter values, all scores, and integer ranks. Rows are sorted best-first.

**`runner.py` — integration point**

Replace `_build_model` with your LASER model constructor. The adapter expects:
- `model.params` — mutable namespace LASER processes read for `beta`, `gravity_k`, `seasonal_amplitude`
- `model.run(n_steps)` — advances the simulation
- `model.patches.new_infections` — `(n_patches, n_steps)` array of weekly incidence
