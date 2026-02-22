The framework is complete and validated. Here is a full summary of what was built:

---

## Guinea Worm LASER Calibration Framework

**Files** (`gworm_calibration/` package + entry point):

| File | Lines | Purpose |
|---|---|---|
| `observed_data.py` | 118 | Chad GWEP data 2018-2023; 8-district geography, populations, spatial weights, seasonal profile |
| `gworm_model.py` | 532 | LASER-convention model: 4 components + `PropertySet` + `NodeFrame` |
| `fitness.py` | 300 | Multi-objective fitness: 4 scored objectives |
| `calibration_loop.py` | 378 | LHS sampler, candidate evaluator, calibration loop, persistence |
| `run_calibration.py` | 105 | Entry point (`--quick`, `--samples`, `--seed` flags) |

---

### Architecture: LASER Conventions

**`PropertySet`** stores all 5 calibration parameters with dot-notation access and `+` merge, matching `laser_core.propertyset`.

**`NodeFrame`** holds per-district arrays (`n_districts × 1` scalar or `n_districts × n_ticks` vector), mirroring `laser_core.LaserFrame` for patch data.

**Four components** — each registered at `__init__(model, verbose)` and ticked via `step(tick)`:

| Component | Pre-computed | Purpose |
|---|---|---|
| `ContainmentScale` | ✓ | GWEP efficacy ramp 30%→100% of max, 2018–2023 |
| `SeasonalForcing` | ✓ | Dry-season peak forcing (Nov–Apr in Chad) |
| `GWTransmission` | step | Copepod reservoir W; deterministic ODE dynamics |
| `Reporter` | step | Annual Poisson observation draws per year-end |

**`GuineaWormModel.run()`** iterates 72 monthly ticks, calling each component in registration order.

---

### Transmission Model

The copepod reservoir `W[d]` is updated each tick from local dog shedding:

```
W[d] ← (1 − 0.3)·W[d] + 3·I_d[d]/N_d[d]

λ_dog[d]   = beta_dog · W[d] · w[d] · seasonal · (1 − containment)
λ_human[d] = beta_human · xi · W[d] · w[d] · seasonal · (1 − containment)
```

The `W` architecture gives `R0_dog ∈ [0.2, 2.0]` over the `beta_dog` range, spanning sub- to super-critical. The **dog:human ratio** emerges from `(beta_dog / (beta_human × xi)) × N_d/N_h` — reproducing 50–100:1 with high `beta_dog`, low `beta_human`, small `xi`.

Initialization is set directly from 2018 surveillance data (`I_d = OBS_DOGS[0]/(p_rec_d × 12) = 347`) to avoid dependence on endemic equilibrium (the 2018-2023 window is a declining non-equilibrium epidemic). `W_init[d]` is district-specific from local `I_d[d]`.

---

### Fitness Function (4 objectives, weighted sum)

| Objective | Weight | Metric |
|---|---|---|
| Log-scale case match | 0.40 | Log-RMSE for annual totals (humans + dogs) |
| Dog:human ratio | 0.25 | Log-deviation from 75:1 target + year-by-year RMSE |
| Temporal trend | 0.20 | Pearson correlation on log counts + direction penalty |
| Spatial concentration | 0.15 | Shortfall below 70% of cases in 4 endemic districts |

---

### Calibration Loop

**Latin Hypercube Sampling** over 5D parameter space ensures space-filling coverage (~2.2× more efficient than random sampling). Each candidate runs `n_replicates` stochastic simulations (independent Poisson observation draws from the same deterministic trajectory); fitness is averaged across replicates to reduce noise from small counts.

Quick test (`--quick`, 150 model runs) found a top candidate with fitness **0.230** in seconds. Full calibration: `python run_calibration.py --samples 500 --replicates 5`.

Sources:
- [GitHub - InstituteforDiseaseModeling/laser](https://github.com/InstituteforDiseaseModeling/laser)
- [IDM Tools - IDMOD](https://www.idmod.org/tools/)
- [Institute for Disease Modeling - Wikipedia](https://en.wikipedia.org/wiki/Institute_for_Disease_Modeling)
