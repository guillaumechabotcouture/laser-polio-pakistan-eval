# LASER Three-Skill Pipeline A/B Evaluation Report

**Date:** 2026-02-20
**Skills tested:** epi-model-parametrization, laser-spatial-disease-modeling, modelops-calabaria
**Model:** Claude Opus 4.6 via `claude --print`

## Methodology

- **WITH skill:** Claude Code in project directory with all 3 skills accessible
- **WITHOUT skill:** Claude Code in a **clean temp directory** with no `.claude/skills/`, no existing scripts, no LASER-specific CLAUDE.md — only a minimal project context establishing the Gates Foundation public health purpose
- **5 prompts** of increasing complexity (basic SEIR → full 20-district model)
- **4 scoring dimensions** (0-3 each): API Correctness, Structural Correctness, Domain Adaptation, Completeness

## Score Sheet

| Prompt | Condition | AC | SC | DA | CO | Total |
|--------|-----------|:--:|:--:|:--:|:--:|:-----:|
| P1: Basic 10-patch SEIR | WITH | 3 | 3 | 3 | 3 | **12** |
| P1: Basic 10-patch SEIR | WITHOUT | 3 | 2 | 2 | 3 | **10** |
| P2: Gravity + Seasonality | WITH | 3 | 3 | 3 | 3 | **12** |
| P2: Gravity + Seasonality | WITHOUT | 0 | 1 | 2 | 2 | **5** |
| P3: Vaccination (RI+SIA) | WITH | 3 | 3 | 3 | 3 | **12** |
| P3: Vaccination (RI+SIA) | WITHOUT | 0 | 1 | 2 | 2 | **5** |
| P4: Calibration framework | WITH | 3 | 3 | 3 | 3 | **12** |
| P4: Calibration framework | WITHOUT | 0 | 0 | 0 | 0 | **0** |
| P5: 20-district model | WITH | 3 | 3 | 3 | 3 | **12** |
| P5: 20-district model | WITHOUT | 0 | 1 | 2 | 2 | **5** |

## Summary Statistics

| Dimension | WITH | WITHOUT | Advantage |
|-----------|:----:|:-------:|:---------:|
| API Correctness | 3.00 | 0.60 | **+2.40** |
| Structural Correctness | 3.00 | 1.00 | **+2.00** |
| Domain Adaptation | 3.00 | 1.60 | **+1.40** |
| Completeness | 3.00 | 1.80 | **+1.20** |
| **Total (per prompt)** | **12.00** | **5.00** | **+7.00** |

| Metric | WITH | WITHOUT |
|--------|:----:|:-------:|
| Sum across 5 prompts | **60/60** | **25/60** |
| Percentage | 100% | 41.7% |

## Key Findings

### 1. Framework Abandonment Pattern

The most revealing finding: without skills, Claude correctly used the LASER API for the simplest case (P1, score 10/12) but **progressively abandoned the framework** as complexity increased. From P2 onward, it:

- Created stand-in `PropertySet` and `LaserFrame` classes instead of importing from `laser.core`
- Reimplemented gravity coupling from scratch instead of using `laser.core.migration.gravity()`
- Built a deterministic ODE metapopulation model instead of using LASER's agent-based engine
- Created custom "Phase" classes referencing nonexistent LASER APIs (`Immunization.RoutineImmunization`, `Immunization.Campaign`)

This suggests base Claude has enough training data exposure to recognize LASER's basic `Model`/`SEIR` pattern but lacks depth for advanced features (gravity migration, `ValuesMap` seasonality, custom components with `on_birth` callbacks, `additional_states`).

### 2. Critical Gotcha Avoidance

The skills documented two critical LASER gotchas that the without-skill condition could not avoid:

- **Vaccination state semantics:** LASER's `TransmissionSE` only checks `state == SUSCEPTIBLE`, ignoring `susceptibility`. The with-skill condition correctly sets `state = RECOVERED` or `state = VACCINATED(4)`. The without-skill condition used "leaky vaccine efficacy" (reducing FOI by 50%) — a valid epidemiological concept but not how LASER works.

- **Birthrate units:** LASER's `BirthsByCBR` expects per-1000/year rates. The without-skill P5 pre-converted to daily per-capita (`25/1000/365`), which would cause silent zero-growth failure in LASER.

### 3. Calibration Collapse (P4)

The without-skill condition produced **empty output** for P4 (calibration). By this point, the accumulated technical debt — a non-LASER aggregated model with no reusable components — made structured calibration infeasible within the session.

### 4. No Measles Negative Transfer

Despite the LASER skill being originally documented with measles examples (England & Wales), **no measles parameter bias** was detected in with-skill outputs:

- R0 = 6 throughout (not measles 12-18)
- Infectious period = 28 days (not measles 7-10)
- OPV at 6 weeks (not MCV1 at 9 months)
- Monsoon seasonal forcing (not school-term)
- 1:200 paralysis ratio correctly applied

### 5. Domain Knowledge vs API Knowledge

The without-skill condition showed **reasonable epidemiological knowledge** (R0=6, monsoon forcing, OPV waning ~3yr, Waziristan corridor) — scoring 1.60/3.00 on Domain Adaptation. But this knowledge couldn't be expressed through the LASER framework without the skill's API guidance, resulting in 0.60/3.00 on API Correctness.

## Comparison to Hypothesized Outcomes

| Dimension | Predicted Advantage | Actual Advantage | Assessment |
|-----------|:-------------------:|:----------------:|:----------:|
| API Correctness | +1.5 | **+2.40** | Exceeded — total framework abandonment not predicted |
| Structural Correctness | +1.0 | **+2.00** | Exceeded — shift to ODE model not predicted |
| Domain Adaptation | +0.5 | **+1.40** | Exceeded — without-skill had reasonable epi knowledge but couldn't apply it |
| Completeness | +0.5 | **+1.20** | Exceeded — P4 complete failure not predicted |

The skill advantage was **stronger than hypothesized** across all dimensions, primarily because the without-skill condition didn't just make API errors — it abandoned the framework entirely.

## Files

### With-Skill Scripts
- `scripts/polio_seir_basic_10patch.py` (P1, 16KB)
- `scripts/polio_gravity_seasonal.py` (P2, 22KB)
- `scripts/polio_seir_10patch.py` (P3, 42KB — with SEIRV vaccination)
- `scripts/calibrate_polio.py` (P4, 31KB)
- `scripts/polio_seir_20district.py` (P5, 24KB)

### Without-Skill Scripts (from clean temp dir)
- `eval/outputs/without-skill/polio_seir_model.py` (P1, 7.5KB)
- `eval/outputs/without-skill/laser_polio_model.py` (P2-P3, 25KB — non-LASER ODE model)
- `eval/outputs/without-skill/polio_laser_model.py` (P4-P5, 18KB — non-LASER tau-leap model)
