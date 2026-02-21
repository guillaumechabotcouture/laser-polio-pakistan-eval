# Guinea Worm A/B Test: Stress Test of Restructured LASER Skill

Can a LASER modeling skill designed around measles/polio examples help Claude
Code build a model for a fundamentally different disease it was never trained
on? This evaluation tests the restructured three-layer skill using guinea worm
(Dracunculus medinensis) in Chad.

**Result: WITH skills 57/60 vs WITHOUT skills 29/60. Zero negative transfer.**

## Why Guinea Worm?

The original polio A/B test (60/60 vs 25/60) showed the skill's value, but
polio shares many features with the skill's measles examples (SEIR, vaccination,
respiratory/fecal-oral). Guinea worm is maximally different:

| Aspect | Polio (original test) | Guinea Worm (stress test) |
|--------|:---------------------:|:-------------------------:|
| Compartmental model | SEIR (lasting immunity) | **SEIS (no immunity)** |
| Transmission | Fecal-oral / respiratory | **Water-borne (copepods)** |
| Incubation | ~3 days | **~365 days** |
| Infectious period | ~28 days | ~21 days (worm emergence) |
| R0 | ~6 | **~1.5-2.5** |
| Hosts | Human only | **Human + dogs (dual-host)** |
| Interventions | Vaccination (OPV/IPV) | **Prevention only** |
| Seasonality | Monsoon peak (Jul-Oct) | **Dry season peak (Jan-Mar)** |
| Case numbers | ~36/year paralytic | **~10-30 human, ~500-1500 dog/year** |
| Context | Endemic transmission | **Near-elimination** |

If the skill only works for polio-like diseases, it teaches recipes. If it
works for guinea worm, it teaches **how to use LASER as a framework**.

## Results

| Prompt | WITH | WITHOUT | Gap |
|--------|:----:|:-------:|:---:|
| P1: Dual-host SEIS setup | 12 | 11 | +1 |
| P2: Water coupling + seasonality | 11 | 9 | +2 |
| P3: Prevention interventions | 12 | 1 | +11 |
| P4: Calibration (Carter Center data) | 11 | 0 | +11 |
| P5: Full 10-year integration | 11 | 8 | +3 |
| **Total** | **57** | **29** | **+28** |

### By Dimension

| Dimension | WITH (avg) | WITHOUT (avg) | Advantage |
|-----------|:----------:|:-------------:|:---------:|
| API Correctness | 3.00 | 0.60 | **+2.40** |
| Structural Correctness | 3.00 | 1.40 | **+1.60** |
| Domain Adaptation | 3.00 | 2.00 | **+1.00** |
| Completeness | 2.40 | 1.80 | **+0.60** |

## Key Findings

### 1. Zero Negative Transfer

The skill did not force guinea worm into polio/measles patterns:

- Correctly used **SEIS** (not SEIR) — implemented `SEISRecovery` component
  routing R back to S
- Built **prevention-only** interventions (ABATE, filters, containment,
  tethering) — no vaccination components
- Used **dry-season forcing** (Jan-Mar peak) — not school-term or monsoon
- Modeled **dual-host** dynamics with dogs as dominant reservoir
- Set **R0 ≈ 1.5-2.5** — not polio's 6 or measles' 12-18

### 2. Framework Abandonment Persists Without Skill

Same pattern as the polio test:

| Prompt | Without-skill LASER usage |
|--------|--------------------------|
| P1 | Used real LASER imports (best result) |
| P2 | Abandoned LASER — built 6 custom modules |
| P3 | Failed entirely — asked a question instead of coding |
| P4 | Complete failure — non-answer |
| P5 | Abandoned LASER — numpy-only implementation |

### 3. API Correctness Is the Core Skill Value

The +2.40/3.00 per-prompt API correctness advantage is **identical to the polio
test**. This is the skill's primary contribution regardless of disease.

### 4. Domain Knowledge Is Not the Bottleneck

Both conditions showed strong guinea worm domain knowledge (DA scores: 3.00 vs
2.00). The LLM knows guinea worm biology from training data. The gap is in
**framework usage**, not epidemiological understanding.

## Prompt Design

Prompts were designed for epidemiological realism:

- **SEIS dynamics** explicitly required (no lasting immunity for guinea worm)
- **Real Carter Center surveillance data** for calibration (2018-2023:
  17→6 human cases/year, 1040→526 dog infections/year)
- **Near-elimination initial conditions** (2 infected humans, 50 infected dogs)
- **Water-mediated transmission** biology described in prompts
- **Dogs as dominant reservoir** (~95% of detected infections in Chad)

See [`prompts.md`](prompts.md) for the full prompt documentation.

## Running the Evaluation

```bash
# From the project root:
./eval/guinea-worm/run-eval.sh all     # Run all 5 prompts (both conditions)
./eval/guinea-worm/run-eval.sh 3       # Run single prompt

# Outputs saved to:
#   eval/guinea-worm/outputs/with-skill/prompt-N.md
#   eval/guinea-worm/outputs/without-skill/prompt-N.md
```

Each prompt runs as a one-shot `claude --print` session. The WITH condition
runs in the project directory (skills accessible). The WITHOUT condition runs
in a clean temp directory with `--disable-slash-commands`.

## Scoring

Use [`rubric.md`](rubric.md) to score outputs on four dimensions (0-3 each):

- **API Correctness (AC):** Does the code use real LASER v1.0.0 API?
- **Structural Correctness (SC):** Dual-host tracking, SEIS dynamics,
  proper component ordering?
- **Domain Adaptation (DA):** Guinea worm biology (not generic SEIR)?
- **Completeness (CO):** Runnable code or pseudocode?

## File Structure

```
eval/guinea-worm/
├── README.md                          # This file
├── prompt-{1-5}.txt                   # Plain-text prompts for claude --print
├── prompts.md                         # Consolidated prompt documentation
├── rubric.md                          # 4-dimension scoring rubric (SEIS-adapted)
├── run-eval.sh                        # Automated A/B runner
└── outputs/
    ├── ab-test-report.md              # Full analysis report with scoring
    ├── guinea_worm_dual_host.png      # Diagnostic plot (P2 with-skill)
    ├── guinea_worm_interventions.png  # Intervention plot (P5 with-skill)
    ├── with-skill/prompt-{1-5}.md     # WITH skill outputs
    └── without-skill/prompt-{1-5}.md  # WITHOUT skill outputs
```

The WITH-skill sessions also generated working Python scripts in `scripts/`:

| File | Lines | Description |
|------|------:|-------------|
| `guinea_worm_components.py` | 943 | 10 custom LASER components (SEIS, dual-host, interventions) |
| `guinea_worm_chad.py` | 565 | Main dual-host SEIS simulation |
| `guinea_worm_chad_interventions.py` | ~800 | Full model with all 4 interventions |
| `guinea_worm_calibration.py` | ~500 | Calibration framework with Carter Center data |

## Comparison to Polio Test

| Metric | Polio | Guinea Worm |
|--------|:-----:|:-----------:|
| WITH score | 60/60 (100%) | 57/60 (95%) |
| WITHOUT score | 25/60 (42%) | 29/60 (48%) |
| Advantage | +35 | +28 |
| Negative transfer | Minor | **None** |
| API correctness gap | +2.40 | +2.40 |
| Framework abandonment | P2+ without skill | P2+ without skill |

The consistent +2.40 API correctness advantage across two very different
diseases confirms the skill teaches **framework mastery**, not disease recipes.
