# Guinea Worm Evaluation Rubric

## Scoring Dimensions (0-3 each, max 12 per prompt, 60 total)

### API Correctness (AC)
How well does the code use the actual LASER v1.0.0 API?

| Score | Criteria |
|-------|----------|
| 0 | Imports don't exist, function signatures invented, would crash immediately |
| 1 | Some real imports but key signatures wrong (wrong arg order, missing required args) |
| 2 | Imports and signatures mostly correct, minor issues (deprecated patterns, wrong defaults) |
| 3 | All imports valid, signatures match v1.0.0 API, would pass import/instantiation |

**Check specifically:**
- [ ] `from laser.generic import Model` (not invented path)
- [ ] `from laser.generic import SEIR` or direct component imports
- [ ] Distribution creation: `dists.gamma(shape=, scale=)` not `dists.gamma(a=, b=)`
- [ ] `Model(scenario, params)` with GeoDataFrame scenario
- [ ] Component signatures: `SEIR.Transmission(model, expdurdist, seasonality=...)`
- [ ] `gravity(populations, distances, k, a, b, c)` — 6 positional args
- [ ] `row_normalizer(network, max_fraction)` — not invented name

### Structural Correctness (SC)
Is the model assembled correctly as a working simulation?

| Score | Criteria |
|-------|----------|
| 0 | Missing critical pieces (no model instantiation, no components, no run call) |
| 1 | Has model + components but wrong assembly, missing host, or wrong compartmental structure (SEIR instead of SEIS) |
| 2 | Both hosts represented with SEIS dynamics, minor structural issues (component order, missing init) |
| 3 | Correct assembly: dual-host tracking, SEIS dynamics (R→S recycling), proper component ordering, model.run() |

**Check specifically:**
- [ ] Scenario is a GeoDataFrame with population, S, E, I, R, geometry columns
- [ ] SEIS dynamics implemented (recovered → susceptible, not absorbing R)
- [ ] Both human and dog populations tracked (dual-host)
- [ ] Separate model instances or property sets per host species
- [ ] Components include S, E, I transitions with R→S recycling
- [ ] Component ordering preserves population invariant
- [ ] Parameters include nticks, beta at minimum
- [ ] model.run() called after setup

### Domain Adaptation (DA)
How well does it adapt to guinea worm biology and Chad's epidemiological context?

| Score | Criteria |
|-------|----------|
| 0 | Generic SEIR with no guinea worm considerations |
| 1 | Correct R0/durations but uses SEIR (with lasting immunity) or misses key biology |
| 2 | Uses SEIS + addresses 2+ of: long pre-patent, dual-host, water-borne transmission, dry-season peak, prevention-only, near-elimination context |
| 3 | Comprehensive adaptation: SEIS dynamics, water-mediated transmission, dual-host with dogs dominant, realistic case counts, appropriate caveats about stochastic dynamics at low numbers |

**Guinea worm elements to look for:**
- [ ] SEIS dynamics (no lasting immunity, R → S recycling)
- [ ] R0 ≈ 1.5-2.5 (not 6+ like polio or 12-18 like measles)
- [ ] Pre-patent period ~365 days (not 3-10 days like most infections)
- [ ] Worm emergence period ~14-21 days
- [ ] Dogs as dominant reservoir in Chad (~95% of detected infections)
- [ ] Water-mediated transmission via copepods (not airborne/fecal-oral)
- [ ] Dry-season peak (Jan-Mar) when stagnant water concentrates copepods
- [ ] Prevention-only interventions (no vaccine exists)
- [ ] Near-elimination case counts (~10-30 human, ~500-1500 dog/year)
- [ ] Acknowledges stochastic dynamics challenge at low case counts
- [ ] No lasting immunity — reinfection occurs
- [ ] Asymmetric cross-species transmission (dogs > humans)

### Completeness (CO)
Is the output actually runnable code or hand-wavy pseudocode?

| Score | Criteria |
|-------|----------|
| 0 | Pseudocode, prose description, or fragments only |
| 1 | Partial code — major pieces missing (e.g., scenario creation, parameter definition) |
| 2 | Nearly complete — could run with minor fixes (missing imports, typos) |
| 3 | Complete, self-contained script that could execute given LASER + data |

## Execution Testing

After scoring, attempt to run each output:

```bash
pip install laser-generic  # ensure v1.0.0+ installed
python eval/guinea-worm/outputs/with-skill/prompt-N.py 2>&1 | head -50
python eval/guinea-worm/outputs/without-skill/prompt-N.py 2>&1 | head -50
```

Record for each:
- **Runs?** Yes / No
- **Error type if no:** ImportError, TypeError, AttributeError, ValueError, RuntimeError
- **How far does it get?** Import / Model init / Component setup / Run start / Completion

## Score Sheet Template

| Prompt | Condition | AC | SC | DA | CO | Total | Runs? | Error |
|--------|-----------|----|----|----|----|-------|-------|-------|
| 1 | with-skill | | | | | /12 | | |
| 1 | without-skill | | | | | /12 | | |
| 2 | with-skill | | | | | /12 | | |
| 2 | without-skill | | | | | /12 | | |
| 3 | with-skill | | | | | /12 | | |
| 3 | without-skill | | | | | /12 | | |
| 4 | with-skill | | | | | /12 | | |
| 4 | without-skill | | | | | /12 | | |
| 5 | with-skill | | | | | /12 | | |
| 5 | without-skill | | | | | /12 | | |

## Aggregate Scoring

| Dimension | WITH Skill (/15) | WITHOUT Skill (/15) | Difference |
|-----------|-----------------|--------------------|-----------:|
| API Correctness | | | |
| Structural Correctness | | | |
| Domain Adaptation | | | |
| Completeness | | | |
| **TOTAL** | **/60** | **/60** | |

## Expected Outcome Hypotheses

| Dimension | Expected skill advantage | Rationale |
|-----------|------------------------|-----------|
| API Correctness | **Strong** (+2.0) | Same framework abandonment pattern expected without skill |
| Structural Correctness | **Moderate** (+1.5) | Dual-host + SEIS are novel; skill helps with component architecture but has no direct example |
| Domain Adaptation | **Weak** (+0.5) | Guinea worm epi knowledge comes from training data and prompts explicitly state SEIS; skill might even bias toward SEIR patterns |
| Completeness | **Moderate** (+1.0) | Skill templates help produce runnable code |

## Negative Transfer Analysis

**Key question:** Does the restructured Layer 1 (discipline) help even when the
disease doesn't match the skill's examples? The validation mindset, assertion
patterns, and "don't reimplement" guidance are disease-agnostic.

**Watch for negative transfer indicators:**
- [ ] SEIR instead of SEIS (ignoring "no lasting immunity" instruction)
- [ ] Vaccination components when no vaccine exists for guinea worm
- [ ] School-term seasonality instead of dry-season peak
- [ ] Gravity migration without water-source adaptation
- [ ] Measles/polio-like R0 (6-18) instead of guinea worm R0 (1.5-2.5)
- [ ] Short incubation periods (days) instead of 365-day pre-patent period
- [ ] Epidemic dynamics instead of near-elimination context
- [ ] Single-host model ignoring the dog reservoir
- [ ] Airborne/fecal-oral transmission framing instead of water-borne

**Positive transfer indicators:**
- [ ] Correct LASER API usage regardless of disease
- [ ] Proper component ordering and model assembly
- [ ] Appropriate use of custom components for SEIS and dual-host
- [ ] Validation assertions for model state
- [ ] Acknowledges model limitations and simplifications
- [ ] Discusses stochastic dynamics challenge at low case counts
- [ ] Correctly implements R→S recycling for SEIS
