# Guinea Worm A/B Test: Stress Test of Restructured LASER Skill

**Date:** 2026-02-21
**Skills tested:** epi-model-parametrization, laser-spatial-disease-modeling, modelops-calabaria
**Model:** Claude Opus 4.6 via `claude --print`
**Disease:** Guinea worm (Dracunculus medinensis) in Chad — not mentioned in any skill

## Purpose

This evaluation stress-tests the restructured LASER skill (3-layer:
discipline/scaffolding/verification) using a fundamentally different disease.
The skill never mentions guinea worm. Key differences from the original polio
test: SEIS dynamics (no immunity), dual-host (dogs + humans), water-mediated
transmission, prevention-only interventions, near-elimination context.

## Methodology

- **WITH skill:** Claude Code in project directory with all 3 skills accessible
- **WITHOUT skill:** Claude Code in a clean temp directory with `--disable-slash-commands`,
  no `.claude/skills/`, no existing scripts — only a minimal CLAUDE.md establishing
  the Gates Foundation guinea worm eradication context
- **5 prompts** of increasing complexity (basic SEIS → full 10-year intervention model)
- **4 scoring dimensions** (0-3 each): API Correctness, Structural Correctness,
  Domain Adaptation, Completeness

### Key Design Choice: SEIS Not SEIR

Guinea worm confers no lasting immunity — recovered individuals return to
susceptible. This is a critical test: the skill teaches SEIR patterns, but the
prompts explicitly require SEIS dynamics. Can the skill help without forcing
the wrong compartmental model?

## Score Sheet

| Prompt | Condition | AC | SC | DA | CO | Total | Notes |
|--------|-----------|:--:|:--:|:--:|:--:|:-----:|-------|
| P1: Dual-host SEIS | WITH | 3 | 3 | 3 | 3 | **12** | Output timed out but generated complete working code\* |
| P1: Dual-host SEIS | WITHOUT | 2 | 3 | 3 | 3 | **11** | Strong result; used real LASER imports + custom InfectiousIS |
| P2: Water coupling + season | WITH | 3 | 3 | 3 | 2 | **11** | Describes architecture; code in project files |
| P2: Water coupling + season | WITHOUT | 1 | 2 | 3 | 3 | **9** | Built custom system "mirrors LASER framework" |
| P3: Prevention interventions | WITH | 3 | 3 | 3 | 3 | **12** | Complete component assembly with FOI formula |
| P3: Prevention interventions | WITHOUT | 0 | 0 | 1 | 0 | **1** | Asked clarifying question instead of producing code |
| P4: Calibration (Carter Ctr) | WITH | 3 | 3 | 3 | 2 | **11** | 7-metric scoring, Sobol sampling, real data |
| P4: Calibration (Carter Ctr) | WITHOUT | 0 | 0 | 0 | 0 | **0** | Non-answer (background task confusion) |
| P5: Full 10-year model | WITH | 3 | 3 | 3 | 2 | **11** | All features, 10 verification checks, diagnostic plots |
| P5: Full 10-year model | WITHOUT | 0 | 2 | 3 | 3 | **8** | Complete but numpy-only (no LASER) |

\* P1 WITH timed out (18 bytes in output), but session created `guinea_worm_chad.py`
(565 lines) + `guinea_worm_components.py` (943 lines), both using real LASER API
with SEIS dynamics, dual-host, gravity coupling. Scored based on generated code.

## Summary Statistics

| Dimension | WITH (/15) | WITHOUT (/15) | Advantage |
|-----------|:----------:|:-------------:|:---------:|
| API Correctness | 15 (3.00) | 3 (0.60) | **+2.40** |
| Structural Correctness | 15 (3.00) | 7 (1.40) | **+1.60** |
| Domain Adaptation | 15 (3.00) | 10 (2.00) | **+1.00** |
| Completeness | 12 (2.40) | 9 (1.80) | **+0.60** |
| **Total** | **57/60** | **29/60** | **+28** |

| Metric | WITH | WITHOUT |
|--------|:----:|:-------:|
| Sum across 5 prompts | **57/60** | **29/60** |
| Percentage | 95.0% | 48.3% |

## Key Findings

### 1. No Negative Transfer — The Skill Adapted Correctly

The restructured skill showed **zero negative transfer** to guinea worm:

| Fear | Observed | Evidence |
|------|----------|----------|
| SEIR forced instead of SEIS | **No** | All with-skill outputs correctly implement SEIS with SEISRecovery (R→S) |
| Vaccination components | **No** | All with-skill outputs use prevention-only interventions |
| School-term seasonality | **No** | All use dry-season forcing (Jan-Mar peak, Jun-Sep trough) |
| Polio/measles R0 | **No** | Human R0≈0.4 (subcritical alone), Dog R0≈2.1, with cross-species coupling |
| Short incubation | **No** | Correctly uses 360-day pre-patent (gamma shape=12, scale=30) |

**This is the most important finding.** The skill's Layer 1 (discipline) — which
emphasizes API validation, assertion patterns, and "don't reimplement" — proved
to be disease-agnostic. It helped Claude use LASER correctly without biasing
toward the skill's measles/polio examples.

### 2. Framework Abandonment Pattern Persists Without Skill

The same pattern from the polio test recurred:

| Prompt | Without-skill LASER usage |
|--------|--------------------------|
| P1 | Used real LASER imports (TransmissionSE + custom InfectiousIS) |
| P2 | Abandoned LASER: "mirrors LASER framework" (6 custom modules) |
| P3 | Failed to produce code (asked clarifying question) |
| P4 | Complete failure (non-answer) |
| P5 | Abandoned LASER: "using only numpy" |

Without the skill, Claude used LASER for the simplest case but progressively
abandoned it as complexity increased — identical to the polio test pattern.
This confirms framework abandonment is a systematic failure mode, not a
one-off issue.

### 3. API Correctness Dominates the Gap

| Dimension | Advantage | % of total gap |
|-----------|:---------:|:--------------:|
| API Correctness | +12 | 42.9% |
| Structural Correctness | +8 | 28.6% |
| Domain Adaptation | +5 | 17.9% |
| Completeness | +3 | 10.7% |

API Correctness accounts for the largest share of the gap, driven by the
without-skill condition abandoning LASER entirely on P2-P5.

### 4. Domain Knowledge Was Strong in Both Conditions

When the without-skill condition produced substantive output (P1, P2, P5), its
Domain Adaptation scores were excellent (3, 3, 3). Both conditions correctly
implemented:
- SEIS dynamics (no lasting immunity)
- Dual-host with dogs as dominant reservoir
- Dry-season transmission forcing
- Water-mediated transmission concepts
- Near-elimination initial conditions

**The skill provides no guinea worm domain advantage** — that knowledge comes
from the LLM's training data and the explicit prompt guidance. The skill's
advantage is entirely in **framework usage** (API + structure).

### 5. With-Skill Produced Cumulative, Production-Quality Code

The with-skill sessions iteratively built a comprehensive codebase:

| File | Lines | Created by |
|------|------:|------------|
| `guinea_worm_components.py` | 943 | P1-P3 sessions |
| `guinea_worm_chad.py` | 565 | P1-P2 sessions |
| `guinea_worm_chad_interventions.py` | ~800 | P5 session |
| `guinea_worm_calibration.py` | ~500 | P4 session |

Features:
- 10 custom LASER components (SEISRecovery, DualHostWaterTransmission,
  ABATELarvicide, WaterFilterDistribution, CaseContainment, DogTethering,
  IntervenedWaterTransmission, IntervenedDualHostTransmission, PatchImportation,
  build_chad_dry_season)
- Diagnostic plots showing correct seasonal forcing, epidemic curves,
  compartment dynamics, intervention scale-up
- Verification checks (SEIS R≈0, non-negativity, infections>0, spatial spread)
- Calibration framework with Carter Center data and 7-metric scoring

The without-skill condition couldn't build cumulatively (temp dirs deleted
between sessions), though some individual outputs (P1, P5) were self-contained
and functional.

## Comparison to Polio A/B Test

| Metric | Polio Test | Guinea Worm Test |
|--------|:----------:|:----------------:|
| WITH score | 60/60 (100%) | 57/60 (95%) |
| WITHOUT score | 25/60 (41.7%) | 29/60 (48.3%) |
| Advantage | +35 (+58.3%) | +28 (+46.7%) |
| Negative transfer | Minor (CCS framing) | **None observed** |
| Framework abandonment | P2+ without skill | P2+ without skill |
| Domain adaptation gap | +1.40 | +1.00 |
| API correctness gap | +2.40 | +2.40 |

**The guinea worm test shows a slightly smaller gap** (+28 vs +35), driven by:
1. Stronger without-skill P1 performance (11 vs 10) — the explicit SEIS
   instruction in the prompt helped
2. No P4 with-skill boost (timeout reduced P1 output score, P4 and P5
   gave 2/3 on completeness due to summary-only outputs)
3. Both conditions benefited from explicit disease biology in prompts

**The API correctness gap is identical** (+2.40 per prompt), confirming this is
the skill's core contribution regardless of disease.

## Hypothesis Verification

| Dimension | Predicted | Actual | Correct? |
|-----------|:---------:|:------:|:--------:|
| API Correctness | Strong (+2.0) | **+2.40** | Yes (slightly underestimated) |
| Structural Correctness | Moderate (+1.5) | **+1.60** | Yes |
| Domain Adaptation | Weak (+0.5) | **+1.00** | Partially (advantage from using LASER's SEIR as SEIS base) |
| Completeness | Moderate (+1.0) | **+0.60** | Partially (without-skill P1/P5 were complete) |

## Conclusion

The restructured LASER skill successfully transfers to a disease it was never
designed for. The skill's value is **framework mastery, not disease recipes**:

1. **Disease-agnostic discipline works.** Layer 1's "don't reimplement,"
   assertion patterns, and API validation helped Claude correctly use LASER
   for guinea worm without any guinea worm-specific guidance.

2. **No negative transfer.** The skill did not force SEIR, vaccination,
   school-term seasonality, or other polio/measles patterns onto guinea worm.
   Claude correctly adapted to SEIS, prevention-only, dry-season, dual-host.

3. **Framework abandonment is the dominant failure mode.** Without the skill,
   Claude knows guinea worm biology well but can't sustainably use the LASER
   framework beyond the simplest case. This pattern is identical to the polio
   test and appears to be a fundamental limitation of LLM code generation
   for niche APIs.

4. **The skill's core value is API correctness (+2.40/3.00 per prompt).**
   This is consistent across both diseases. Domain adaptation is a secondary
   benefit — the LLM already knows disease biology from training data.

### Score: WITH 57/60 vs WITHOUT 29/60
