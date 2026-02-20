# Evaluation Rubric

## Scoring Dimensions (0-3 each)

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
| 1 | Has model + components but wrong assembly (missing compartments, broken flow) |
| 2 | Complete assembly, minor structural issues (component order wrong, missing init) |
| 3 | Correctly ordered components, valid scenario GeoDataFrame, proper initialization, model.run() |

**Check specifically:**
- [ ] Scenario is a GeoDataFrame with population, S, E, I, R, geometry columns
- [ ] Components include Susceptible, Exposed, Infectious, Recovered, Transmission
- [ ] Component ordering preserves S+E+I+R=N invariant
- [ ] Parameters include nticks, beta at minimum
- [ ] model.run() called after setup

### Domain Adaptation (DA)
How well does it adapt to polio (vs blindly copying measles patterns)?

| Score | Criteria |
|-------|----------|
| 0 | Generic SEIR with no polio-specific considerations |
| 1 | Correct R0/durations for polio but otherwise generic |
| 2 | Addresses 2+ of: paralysis ratio, asymptomatic transmission, OPV dynamics, age structure, monsoon seasonality |
| 3 | Comprehensive polio adaptation with appropriate caveats about model limitations |

**Polio-specific elements to look for:**
- [ ] R0 ≈ 6 (not 12-18 like measles)
- [ ] Latent ~3 days, infectious ~28 days (not 10/8 like measles)
- [ ] Mentions 1:200 paralysis ratio or AFP underascertainment
- [ ] Seasonal forcing aligned to monsoon (Jul-Oct), not school terms
- [ ] OPV mucosal immunity vs IPV humoral immunity distinction
- [ ] Notes that basic SEIR may be insufficient for polio's immunity complexity

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
python eval/outputs/with-skill/prompt-N.py 2>&1 | head -50
python eval/outputs/without-skill/prompt-N.py 2>&1 | head -50
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

## Expected Outcome Hypotheses

| Dimension | Expected skill advantage | Rationale |
|-----------|------------------------|-----------|
| API Correctness | **Strong** (+1.5 avg) | Without skill, Claude likely hallucinates pre-v1.0.0 or nonexistent APIs |
| Structural Correctness | **Strong** (+1.0 avg) | Component ordering, GeoDataFrame setup, network init are non-obvious |
| Domain Adaptation | **Weak** (+0.5 avg) | Polio epi knowledge comes from training data; skill might even bias toward measles patterns |
| Completeness | **Moderate** (+0.5 avg) | Skill's code examples provide template; without it, more pseudocode likely |

**Negative transfer risk:** Prompts 3 and 5 may show the skill causing Claude to shoehorn polio into the measles SEIR pattern (school-term seasonality, CCS framing) rather than adapting appropriately.
