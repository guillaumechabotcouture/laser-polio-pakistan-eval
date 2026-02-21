# Evaluating Claude Code Skills for Spatial Disease Modeling

Can specialized AI coding skills help Claude Code build epidemiological models
using a niche simulation framework? This repository provides an A/B evaluation
measuring the impact of a three-skill pipeline on code generation quality when
building a spatial polio transmission model for Pakistan using the
[LASER framework](https://laser.idmod.org/).

**Result: WITH skills scored 60/60 vs WITHOUT skills 25/60 (41.7%).**

## Why This Matters

Large language models have broad but shallow knowledge of specialized scientific
software. The [LASER framework](https://github.com/laser-base/laser-generic)
(Light Agent Spatial modeling for ERadication) is a niche agent-based disease
modeling toolkit from the Institute for Disease Modeling. Without guidance,
Claude correctly uses LASER's basic API but **abandons the framework entirely**
when complexity increases, falling back to hand-rolled ODE simulations that
cannot answer the underlying public health questions.

This evaluation quantifies that gap and measures whether Claude Code
[skills](https://docs.anthropic.com/en/docs/claude-code/skills) — structured
reference documents loaded into context — can bridge it.

## Results at a Glance

| Dimension | WITH Skills | WITHOUT Skills | Advantage |
|-----------|:----------:|:-------------:|:---------:|
| API Correctness | 3.00 / 3 | 0.60 / 3 | **+2.40** |
| Structural Correctness | 3.00 / 3 | 1.00 / 3 | **+2.00** |
| Domain Adaptation | 3.00 / 3 | 1.60 / 3 | **+1.40** |
| Completeness | 3.00 / 3 | 1.80 / 3 | **+1.20** |
| **Total (per prompt avg)** | **12 / 12** | **5 / 12** | **+7.00** |

### Key Findings

1. **Framework abandonment is the dominant failure mode.** Without skills,
   Claude used the LASER API correctly for the simplest prompt but abandoned it
   from Prompt 2 onward, reimplementing everything as a custom ODE model with
   polyfill classes. By Prompt 4 (calibration), it produced empty output.

2. **Both documented "Critical Gotchas" were validated.** The without-skill
   code hit both silent failure patterns documented in the skill: wrong
   birthrate units (daily per-capita instead of per-1000/year) and wrong
   vaccination mechanism (modifying `susceptibility` instead of `state`). Both
   produce plausible-looking but quantitatively wrong results with no runtime
   errors.

3. **No negative transfer from measles examples.** Despite the skill being
   originally documented with measles examples (England & Wales), the with-skill
   outputs correctly adapted to polio: R0=6 (not measles 12-18), 28-day
   infectious period (not 7-10), OPV at 6 weeks (not MCV1 at 9 months), and
   monsoon seasonal forcing (not school-term).

4. **112 vs 11 correct LASER API calls** (10.2x ratio). The with-skill
   condition made zero API hallucinations across all 5 prompts. The without-skill
   condition made 11 correct calls (all in Prompt 1), then zero from Prompt 2
   onward.

5. **4-8 hours vs 2-6 days to production.** An expert modeler could refine the
   with-skill output into a production-ready calibrated model in hours. The
   without-skill output would require days of rewriting — a 12x time multiplier.

## Evaluation Design

### Three-Skill Pipeline

The evaluation tests a pipeline of three Claude Code skills:

| Skill | Purpose | Trigger |
|-------|---------|---------|
| **epi-model-parametrization** | Research disease parameters, identify calibration targets, guide model structure | "what parameters do I need", "find R0 for" |
| **laser-spatial-disease-modeling** | Build spatial SEIR models using LASER, wrap as BaseModel | "LASER model", "spatial disease model" |
| **modelops-calabaria** | Calibrate, optimize, run scenarios, scale to cloud | "calibrate model", "parameter sweep" |

### Test Conditions

- **WITH skills:** Claude Code session in the project directory with all 3
  skills accessible, plus existing project code and references
- **WITHOUT skills:** Claude Code session in a **clean temporary directory**
  with no `.claude/skills/`, no existing scripts, no LASER-specific knowledge —
  only a minimal CLAUDE.md establishing the public health context

### Five Prompts (Increasing Complexity)

| # | Prompt | Tests | Difficulty |
|---|--------|-------|:----------:|
| 1 | Basic 10-patch SEIR | LASER API, imports, component signatures | Low |
| 2 | Gravity network + monsoon forcing | `gravity()`, `ValuesMap`, spatial coupling | Medium |
| 3 | OPV vaccination (RI + SIA) | Immunization API, custom components, OPV waning | Medium-High |
| 4 | Calibration framework | Loss functions, AFP surveillance, parameter sampling | High |
| 5 | Full 20-district integration | End-to-end model with all features | Highest |

### Scoring Rubric

Each prompt is scored 0-3 on four dimensions (max 12 per prompt, 60 total):

- **API Correctness (AC):** Does the code use real LASER v1.0.0 API?
- **Structural Correctness (SC):** Is the model assembled correctly?
- **Domain Adaptation (DA):** Does it adapt to polio (not blindly copy measles)?
- **Completeness (CO):** Is it runnable code or pseudocode?

Full rubric: [`eval/rubric.md`](eval/rubric.md)

## Running the Evaluation

### Prerequisites

```bash
pip install laser-generic    # LASER framework v1.0.0+
pip install claude-code      # Claude Code CLI (or install via npm)
```

### Run All Prompts

```bash
# From the project directory:
./eval/run-eval.sh all

# Or run a single prompt (both conditions):
./eval/run-eval.sh 3
```

The script runs each prompt through two Claude Code sessions:
1. **WITH skills** — normal session in the project directory
2. **WITHOUT skills** — session with `--disable-slash-commands` in a clean temp dir

Outputs are saved to `eval/outputs/{with,without}-skill/prompt-N.md`.

### Score the Outputs

Use the rubric in [`eval/rubric.md`](eval/rubric.md) to score each output on
the four dimensions. The existing scores are in
[`eval/outputs/ab-test-report.md`](eval/outputs/ab-test-report.md).

## Repository Structure

```
laser-polio-pakistan-eval/
├── README.md                          # This file
├── CLAUDE.md                          # Project context for Claude Code
│
├── .claude/skills/                    # The three-skill pipeline under test
│   ├── laser-spatial-disease-modeling/
│   │   ├── SKILL.md                   #   LASER model building (Steps 1-8)
│   │   ├── references/
│   │   │   ├── laser_api_reference.md #   Complete LASER v1.0.0 API docs
│   │   │   └── wavelet_analysis.md    #   Wavelet phase analysis guide
│   │   └── scripts/
│   │       ├── custom_components.py   #   Importation, VaccinationCampaign
│   │       ├── calibration_metrics.py #   CCS + wavelet + loss bridge
│   │       └── laser_basemodel.py     #   BaseModel wrapper template
│   ├── epi-model-parametrization/
│   │   └── SKILL.md                   #   Parameter research guide
│   └── modelops-calabaria/
│       ├── SKILL.md                   #   Calibration workflow
│       └── references/
│           └── modelops_calabaria_reference.md
│
├── eval/                              # Evaluation framework
│   ├── prompts.md                     #   5 test prompts
│   ├── rubric.md                      #   4-dimension scoring rubric
│   ├── run-eval.sh                    #   Automated A/B runner
│   ├── prompt-{1-5}.txt              #   Extracted plain-text prompts
│   └── outputs/
│       ├── ab-test-report.md          #   Main A/B evaluation report
│       ├── deep-analysis-report.md    #   API audit + failure taxonomy
│       ├── model-calibration-report.md#   13-district calibration write-up
│       ├── with-skill/                #   WITH skill outputs (5 summaries)
│       ├── without-skill/             #   WITHOUT skill outputs + scripts
│       └── *.png, *.csv              #   Diagnostic plots + data
│
├── scripts/                           # WITH-skill generated models
│   ├── polio_seir_basic_10patch.py    #   P1: Basic SEIR
│   ├── polio_gravity_seasonal.py      #   P2: Gravity + monsoon
│   ├── polio_seir_10patch.py          #   P3: Full SEIRV with vaccination
│   ├── calibrate_polio.py             #   P4: LHS calibration framework
│   ├── polio_seir_20district.py       #   P5: 20-district integration
│   ├── custom_components.py           #   Shared custom LASER components
│   └── polio_calibration.py           #   Alternative calibration script
│
└── reference/
    └── polio-model-requirements.md    # Ground truth domain specification
```

## Reports

| Report | Description |
|--------|-------------|
| [`eval/outputs/ab-test-report.md`](eval/outputs/ab-test-report.md) | Main A/B evaluation: score sheet, per-prompt analysis, summary statistics |
| [`eval/outputs/deep-analysis-report.md`](eval/outputs/deep-analysis-report.md) | Deep dive: API usage audit (112 vs 11 calls), failure mode taxonomy, skill-to-code mapping, practical impact assessment |
| [`eval/outputs/model-calibration-report.md`](eval/outputs/model-calibration-report.md) | 13-district Pakistan polio model calibration against MMWR/GPEI surveillance data |

## The LASER Framework

[LASER](https://laser.idmod.org/) (Light Agent Spatial modeling for
ERadication) is an open-source agent-based disease modeling framework from the
Institute for Disease Modeling (IDM). It provides:

- **Agent-based SEIR dynamics** with per-agent state tracking
- **Gravity-model spatial coupling** between geographic patches
- **Seasonal forcing** via `ValuesMap` time series
- **Vaccination components** (`RoutineImmunizationEx`, custom campaigns)
- **Vital dynamics** (births by CBR, mortality by CDR/estimator)

Install: `pip install laser-generic` (includes `laser-core` and `laser-generic`)

Documentation: [laser.idmod.org](https://laser.idmod.org/laser-generic/)

## Project Context

This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation's Global Health division. The work supports polio eradication
efforts by building spatial transmission models to evaluate vaccination
strategies across Pakistan's districts. All modeling uses standard SEIR
compartmental frameworks — purely computational epidemiology for public health
planning.

## License

This evaluation framework and analysis are provided for research and educational
purposes.
