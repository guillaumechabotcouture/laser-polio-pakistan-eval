# Difficulty Curve Evaluation Rubric

## Scoring: API Correctness Only (0-3)

This suite scores on a single dimension — API Correctness — since the goal is
to measure how LASER API knowledge degrades as prompt complexity increases.
Structural correctness, domain adaptation, and completeness are deliberately
excluded to keep the measurement clean.

### API Correctness (AC)

| Score | Criteria |
|-------|----------|
| 0 | Abandons LASER entirely (invents framework, uses scipy/mesa/etc.) or imports don't exist |
| 1 | Some real LASER imports but key signatures wrong (wrong arg order, missing required args, invented parameters) |
| 2 | Mostly correct API — imports and signatures right, minor issues (deprecated patterns, wrong defaults, slight misuse) |
| 3 | All LASER API calls valid for v1.0.0 — would pass import and instantiation checks |

## Per-Level Check Items

### D1 -- Minimal Model
- [ ] `from laser.generic import Model` (not invented path)
- [ ] `from laser.generic import PropertySet` (not invented)
- [ ] Scenario is a GeoDataFrame with required columns (population, S, E, I, R, geometry)
- [ ] `Model(scenario, params)` — two positional args, GeoDataFrame + PropertySet
- [ ] `PropertySet(nticks=365, prng_seed=42)` — keyword args
- [ ] `model.run()` called

### D2 -- Add SEIR Components
All D1 checks, plus:
- [ ] `from laser.generic import SEIR` (not `from laser.seir import ...`)
- [ ] `from laser.generic import dists`
- [ ] `dists.gamma(shape=4, scale=1)` — named keyword args, not positional
- [ ] `dists.gamma(shape=5, scale=2)` — for infectious period
- [ ] Components added in order: `SEIR.Susceptible`, `SEIR.Exposed`, `SEIR.Infectious`, `SEIR.Recovered`, `SEIR.Transmission`
- [ ] All five components present (not missing Susceptible or Recovered)

### D3 -- Transmission with expdurdist
All D2 checks, plus:
- [ ] `SEIR.Transmission(model, expdurdist=inf_dist)` — expdurdist wired to the infectious period distribution
- [ ] The distribution object passed to expdurdist is `dists.gamma(shape=5, scale=2)` (not the latent period)
- [ ] expdurdist is a keyword argument (not positional-only)

### D4 -- Births and Deaths
All D3 checks, plus:
- [ ] `from laser.generic import BirthsByCBR` (correct import path)
- [ ] `from laser.generic import MortalityByCDR` (correct import path)
- [ ] `from laser.generic import calc_capacity` (correct import path)
- [ ] CBR passed as per-1000/year (value 30, NOT 0.03 or 0.0000822)
- [ ] CDR passed as per-1000/year (value 10, NOT 0.01 or 0.0000274)
- [ ] `calc_capacity` called with birth rate in per-1000/year and years (not daily rate)
- [ ] `nticks=3650` for 10-year simulation

### D5 -- Gravity Network
All D4 checks, plus:
- [ ] Distance computation present (using LASER utility or manual haversine)
- [ ] `from laser.generic import gravity` (correct import)
- [ ] `gravity(populations, distances, k, a, b, c)` — exactly 6 positional args
- [ ] k=0.01, a=1, b=1, c=1.5 — correct parameter values
- [ ] Network set on model (e.g., `model.network = network` or equivalent)

### D6 -- Row Normalization
All D5 checks, plus:
- [ ] `from laser.generic import row_normalizer` (correct import, not invented name)
- [ ] `row_normalizer(network, max_frac)` — two args: network matrix + scalar
- [ ] max_frac=0.15 (not 15 or 0.015)
- [ ] Applied AFTER gravity(), BEFORE setting on model
- [ ] Result assigned back to network or passed directly to model

### D7 -- Seasonal Forcing via ValuesMap
All D6 checks, plus:
- [ ] `from laser.generic import ValuesMap` (correct import)
- [ ] `ValuesMap.from_timeseries(array, num_patches)` — classmethod with 2 args
- [ ] Input array is 365-length (one year of daily values)
- [ ] num_patches=4 (matches scenario)
- [ ] Passed to `SEIR.Transmission(model, ..., seasonality=values_map)`
- [ ] Seasonal pattern has winter peak ~1.3 and summer trough ~0.7

### D8 -- Custom Importation Component
All D7 checks, plus:
- [ ] Custom class with `__init__(self, model)` signature
- [ ] Custom class with `step(self, tick)` signature
- [ ] Component modifies model state correctly (moves agents S->I or adds infections)
- [ ] Tick-based logic: every 60 ticks, 3 infections in patch 0
- [ ] Component instantiated with model as argument (not added via different pattern)

## Execution Testing

After scoring, attempt to run each output:

```bash
pip install laser-generic  # ensure v1.0.0+ installed
python eval/difficulty-curve/outputs/with-skill/prompt-DN.py 2>&1 | head -50
python eval/difficulty-curve/outputs/without-skill/prompt-DN.py 2>&1 | head -50
```

Record for each:
- **Runs?** Yes / No
- **Error type if no:** ImportError, TypeError, AttributeError, ValueError, RuntimeError
- **How far does it get?** Import / Model init / Component setup / Run start / Completion

## Score Sheet Template

Max score per prompt: 3. Max total per condition: 24 (8 prompts x 3).

| Prompt | Condition | AC (0-3) | Runs? | Error | Notes |
|--------|-----------|----------|-------|-------|-------|
| D1 | with-skill | | | | |
| D1 | without-skill | | | | |
| D2 | with-skill | | | | |
| D2 | without-skill | | | | |
| D3 | with-skill | | | | |
| D3 | without-skill | | | | |
| D4 | with-skill | | | | |
| D4 | without-skill | | | | |
| D5 | with-skill | | | | |
| D5 | without-skill | | | | |
| D6 | with-skill | | | | |
| D6 | without-skill | | | | |
| D7 | with-skill | | | | |
| D7 | without-skill | | | | |
| D8 | with-skill | | | | |
| D8 | without-skill | | | | |

## Aggregate Scoring

| Condition | D1 | D2 | D3 | D4 | D5 | D6 | D7 | D8 | Total (/24) |
|-----------|----|----|----|----|----|----|----|----|-------------|
| WITH skill | | | | | | | | | |
| WITHOUT skill | | | | | | | | | |
| Delta | | | | | | | | | |

## Difficulty Curve Analysis

Plot AC score (y-axis, 0-3) vs. difficulty level (x-axis, D1-D8) for both
conditions. Key metrics to extract:

| Metric | WITH Skill | WITHOUT Skill |
|--------|-----------|--------------|
| **Collapse point** (first level with AC < 2) | D__ | D__ |
| **Skill delta** (collapse shift in levels) | -- | -- |
| **Asymptotic floor** (mean AC at D6-D8) | | |
| **Area under curve** (sum of all 8 AC scores) | /24 | /24 |

## Expected Outcome Hypotheses

| Level Range | Expected Pattern | Rationale |
|-------------|-----------------|-----------|
| D1-D2 | Both conditions score 2-3 | Basic SEIR is well-represented in training data |
| D3 | Minor divergence begins | expdurdist wiring is LASER-specific |
| D4-D5 | Clear divergence | BirthsByCBR, calc_capacity, gravity() are LASER-specific with non-obvious signatures |
| D6-D7 | Without-skill collapses | row_normalizer and ValuesMap are niche APIs unlikely in training data |
| D8 | Both potentially recover slightly | Custom components follow generic Python patterns |

**Predicted collapse points:**
- WITH skill: D7 or D8 (skill provides API reference through D6)
- WITHOUT skill: D4 or D5 (demographic and network APIs are where hallucination begins)

## Framework Abandonment Indicators

Watch for these signs that Claude has abandoned the LASER API:

- **Import fabrication:** `from laser.models import SpatialSEIR` (doesn't exist)
- **Signature invention:** `Model.add_component(Transmission, beta=0.5)` (wrong pattern)
- **Reimplementation:** Writing SEIR ODEs with scipy.integrate instead of using LASER components
- **Library substitution:** Using mesa, epimodels, or custom numpy arrays instead of LASER
- **Partial abandonment:** Correct Model/PropertySet but invented component API
- **Confident hallucination:** Plausible-looking but entirely wrong API calls
