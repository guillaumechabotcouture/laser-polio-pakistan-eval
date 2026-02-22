# Multi-Turn Iterative Evaluation Rubric

## Scoring Dimensions (0-3 each)

### API Correctness (AC)

How well does the code use the actual LASER v1.0.0 API?

| Score | Criteria |
|-------|----------|
| 0 | Imports don't exist, function signatures invented, would crash immediately |
| 1 | Some real imports but key signatures wrong (wrong arg order, missing required args) |
| 2 | Imports and signatures mostly correct, minor issues (deprecated patterns, wrong defaults) |
| 3 | All imports valid, signatures match v1.0.0 API, would pass import/instantiation |

**Check specifically per step:**

Step 1:
- [ ] `from laser.generic import Model` (not invented path)
- [ ] `from laser.generic import SEIR` or direct component imports
- [ ] `Model(scenario, params)` with GeoDataFrame scenario
- [ ] `PropertySet(nticks=365, beta=..., prng_seed=...)`

Step 2:
- [ ] `from laser.generic import gravity, row_normalizer`
- [ ] `gravity(populations, distances, k, a, b, c)` -- 6 positional args
- [ ] `row_normalizer(network, max_fraction)` -- not invented name

Step 3:
- [ ] `from laser.generic import ValuesMap`
- [ ] `ValuesMap.from_timeseries(array, num_patches)` -- correct signature

Step 4:
- [ ] Fix uses correct API (no new hallucinated imports)
- [ ] Error diagnosis addresses actual API misuse if applicable

Step 5:
- [ ] `from laser.generic import BirthsByCBR, MortalityByCDR`
- [ ] CBR/CDR passed as per-1000/year (NOT pre-converted to daily per-capita)
- [ ] `calc_capacity` called with per-1000/year birthrate

### Structural Correctness (SC)

Is the model assembled correctly as a working simulation?

| Score | Criteria |
|-------|----------|
| 0 | Missing critical pieces (no model instantiation, no components, no run call) |
| 1 | Has model + components but wrong assembly (missing compartments, broken flow) |
| 2 | Complete assembly, minor structural issues (component order wrong, missing init) |
| 3 | Correctly ordered components, valid scenario GeoDataFrame, proper initialization, model.run() |

**Check specifically per step:**

Step 1:
- [ ] Scenario is a GeoDataFrame with population, S, E, I, R, geometry columns
- [ ] Components: Susceptible, Exposed, Infectious, Recovered, Transmission (in order)
- [ ] Heterogeneous populations (100k, 200k, 150k, 80k)
- [ ] model.run() called

Step 2:
- [ ] Distance matrix correct for linear arrangement (75, 150, 225 km gaps)
- [ ] Network passed to model via params or scenario
- [ ] Step 1 structure preserved

Step 3:
- [ ] Seasonal array has correct shape (365 values)
- [ ] Seasonality integrated into Transmission component
- [ ] Steps 1-2 structure preserved

Step 4:
- [ ] Fix is structurally sound (doesn't break model assembly)
- [ ] Previous working structure preserved

Step 5:
- [ ] BirthsByCBR and MortalityByCDR added as components
- [ ] calc_capacity used for pre-allocation
- [ ] nticks updated to 3650
- [ ] Population verification logic present

### Domain Adaptation (DA)

How well does it adapt to the specified disease parameters? (Lighter weight for
this suite since we use a generic respiratory disease, not a specific disease.)

| Score | Criteria |
|-------|----------|
| 0 | Ignores specified parameters entirely, uses arbitrary values |
| 1 | Some parameters match specification but key ones wrong (R0, durations) |
| 2 | R0, latent period, and infectious period match specification; initial conditions close |
| 3 | All specified parameters correct, appropriate beta calculation for R0~5, sensible transitions |

**Check specifically:**
- [ ] beta calculated to achieve R0~5 given infectious period of 10 days
- [ ] Latent period ~4 days (gamma-distributed or fixed)
- [ ] Infectious period ~10 days (gamma-distributed or fixed)
- [ ] Initial conditions: 90% S, 1% I, 9% R per patch
- [ ] Gravity parameters match specification: k=0.01, a=1, b=1, c=1.5
- [ ] Seasonal profile matches: winter peak 1.3x, summer trough 0.7x
- [ ] CBR=30, CDR=10 per 1000/year (Step 5)

### Completeness (CO)

Is the output actually runnable code or hand-wavy pseudocode?

| Score | Criteria |
|-------|----------|
| 0 | Pseudocode, prose description, or fragments only |
| 1 | Partial code -- major pieces missing (e.g., scenario creation, parameter definition) |
| 2 | Nearly complete -- could run with minor fixes (missing imports, typos) |
| 3 | Complete, self-contained script that could execute given LASER installed |

---

## Consistency Bonus (CB)

**Scored per step (0 or 1), except Step 1 which has no previous step.**

| Score | Criteria |
|-------|----------|
| 0 | Rewrites from scratch, breaks previous functionality, ignores previous code, or introduces regressions |
| 1 | Correctly builds on previous step's code -- preserves working elements, adds new functionality without regressing |

**Check specifically:**
- Does the code from step N contain all working elements from step N-1?
- Were any previously correct API calls changed to incorrect ones?
- Were any previously working components removed or broken?
- Does the code show incremental modification (not full rewrite)?

---

## Scoring Summary

| Step | AC | SC | DA | CO | CB | Step Max |
|------|----|----|----|----|-----|----------|
| 1 | 0-3 | 0-3 | 0-3 | 0-3 | n/a | **12** |
| 2 | 0-3 | 0-3 | 0-3 | 0-3 | 0-1 | **13** |
| 3 | 0-3 | 0-3 | 0-3 | 0-3 | 0-1 | **13** |
| 4 | 0-3 | 0-3 | 0-3 | 0-3 | 0-1 | **13** |
| 5 | 0-3 | 0-3 | 0-3 | 0-3 | 0-1 | **13** |
| **Total** | | | | | | **64** |

---

## Execution Testing

After scoring, the driver script automatically tests execution at each step.
Record for each:

- **Runs?** Yes / No
- **Error type if no:** ImportError, TypeError, AttributeError, ValueError, RuntimeError
- **How far does it get?** Import / Model init / Component setup / Run start / Completion
- **Step 4 diagnosis:** Did it correctly identify the root cause?

---

## Score Sheet Template

### Per-Step Scores

| Step | Condition | AC | SC | DA | CO | CB | Total | Runs? | Error |
|------|-----------|----|----|----|----|-----|-------|-------|-------|
| 1 | with-skill | | | | | -- | /12 | | |
| 1 | without-skill | | | | | -- | /12 | | |
| 2 | with-skill | | | | | | /13 | | |
| 2 | without-skill | | | | | | /13 | | |
| 3 | with-skill | | | | | | /13 | | |
| 3 | without-skill | | | | | | /13 | | |
| 4 | with-skill | | | | | | /13 | | |
| 4 | without-skill | | | | | | /13 | | |
| 5 | with-skill | | | | | | /13 | | |
| 5 | without-skill | | | | | | /13 | | |

### Aggregate Scores

| Dimension | WITH Skill (/max) | WITHOUT Skill (/max) | Difference |
|-----------|------------------|---------------------|----------:|
| API Correctness | /15 | /15 | |
| Structural Correctness | /15 | /15 | |
| Domain Adaptation | /15 | /15 | |
| Completeness | /15 | /15 | |
| Consistency Bonus | /4 | /4 | |
| **TOTAL** | **/64** | **/64** | |

### Key Metrics

| Metric | WITH Skill | WITHOUT Skill |
|--------|-----------|--------------|
| Steps where code runs | /5 | /5 |
| Final step (5) runs? | | |
| Step 4 root cause found? | | |
| Cumulative regressions (Steps 2-5) | | |
| Code coherence at Step 5 (all features present?) | | |

### Progressive Complexity Tracking

Track how scores change across steps to see if one condition degrades faster:

| Step | WITH total | WITHOUT total | Gap (WITH - WITHOUT) |
|------|-----------|--------------|---------------------|
| 1 | /12 | /12 | |
| 2 | /13 | /13 | |
| 3 | /13 | /13 | |
| 4 | /13 | /13 | |
| 5 | /13 | /13 | |

---

## Expected Outcome Hypotheses

| Dimension | Expected skill advantage | Rationale |
|-----------|------------------------|-----------|
| API Correctness | **Strong** (+1.5 avg) | Consistent LASER API knowledge prevents progressive hallucination |
| Structural Correctness | **Moderate** (+1.0 avg) | Skill's component ordering guidance prevents structural drift |
| Domain Adaptation | **Weak** (+0.3 avg) | Generic disease -- domain knowledge less relevant |
| Completeness | **Moderate** (+0.5 avg) | Skill templates help produce self-contained code |
| Consistency | **Strong** (+0.75 avg) | Skill's discipline layer should prevent code rewrites and regressions |

### Multi-Turn Specific Hypotheses

1. **Cumulative consistency:** WITH-skill should maintain coherent code across
   all 5 steps, while WITHOUT-skill may accumulate technical debt (broken
   imports, lost components, structural drift).

2. **Error recovery (Step 4):** WITH-skill should correctly diagnose API
   misuse; WITHOUT-skill may apply superficial patches or introduce new errors.

3. **Progressive degradation:** WITHOUT-skill scores may decline in later steps
   as code complexity exceeds what can be managed without framework knowledge.

4. **Rewrite tendency:** WITHOUT-skill is more likely to rewrite from scratch
   at each step (losing consistency points), while WITH-skill should show
   incremental modification.

---

## Negative Transfer Risks

Since this suite uses a generic respiratory disease (not polio or guinea worm),
negative transfer from the skill's disease-specific examples is less likely.
However, watch for:

- [ ] Skill biasing toward polio/measles parameters when generic ones are
      specified (R0~5, not R0~6 or R0~12-18)
- [ ] Skill adding vaccination components when not requested
- [ ] Skill imposing monsoon seasonality when winter peak is specified
- [ ] Skill over-engineering the model with components not requested in the
      step (adding demographics before Step 5)
