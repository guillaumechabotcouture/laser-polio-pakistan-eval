# Evaluating Claude Code Skills for Spatial Disease Modeling

Can specialized AI coding skills help Claude Code build epidemiological models
using a niche simulation framework? This repository provides A/B evaluations
measuring the impact of a Claude Code
[skill](https://docs.anthropic.com/en/docs/claude-code/skills) on code
generation quality when building spatial disease transmission models using the
[LASER framework](https://laser.idmod.org/).

## Results at a Glance

Four evaluation suites, all run in **Docker-isolated containers** to prevent
the WITHOUT condition from reading LASER source code:

| Suite | WITH Skill | WITHOUT Skill | Gap |
|-------|:----------:|:------------:|:---:|
| **Atomic** (20 API prompts) | 39/40 (97.5%) | 27/40 (67.5%) | +30pp |
| **Difficulty Curve** (8 levels) | 24/24 (100%) | 17/24 (70.8%) | +29pp |
| **Reruns — Polio** (15 sessions) | AC 3.0/3 | AC 0.3/3 | +2.7 |
| **Reruns — Guinea Worm** (15 sessions) | AC 3.0/3 | AC 0.5/3 | +2.5 |
| **Multi-Turn** (5 steps) | 40/64 quality | 24/64 quality | +16 |

**The skill advantage scales with task complexity:**
- Simple API recall: **+30 percentage points**
- Full model building: **+80 percentage points**

### Key Findings

1. **Framework abandonment is deterministic.** Without the skill, Claude
   abandons LASER at a consistent complexity threshold (prompt 2 for polio,
   prompt 3 for guinea worm) and rebuilds in numpy. This reproduces across
   all 3 statistical runs — it's a capability boundary, not random variance.

2. **WITH-skill scores AC=3 on every valid output.** Zero variance across
   48 valid reruns sessions, both diseases. The skill produces perfectly
   consistent API correctness.

3. **The skill enables cross-disease transfer.** Despite being documented
   with measles/polio examples, the skill correctly guides guinea worm
   models (SEIS dynamics, low R0, no vaccination, dry-season forcing).
   WITHOUT the skill, Claude cannot adapt LASER to novel diseases.

4. **WITHOUT-skill Claude is adaptive, not helpless.** Rather than producing
   broken LASER code, it rationally pivots to numpy/scipy. This produces
   running code (3/5 multi-turn steps) but scores zero on API correctness.

5. **The skill is token-efficient.** It produces 40-62% more output (complete
   code vs summaries) in 30% less wall time on complex tasks.

6. **Spatial coupling APIs are the strongest discriminator.** `gravity()`,
   `row_normalizer()`, and `calc_capacity()` are where WITHOUT fails hardest
   — these have non-obvious signatures not in Claude's training data.

## Evaluation Design

### Docker-Isolated A/B Testing

Both conditions run inside Docker containers to ensure clean separation:

- **WITH skill** (`laser-eval-with`): `laser-generic` installed, project
  directory mounted read-only, skill files accessible
- **WITHOUT skill** (`laser-eval-without`): Clean container, no LASER packages,
  no skill files, minimal CLAUDE.md with project context only

```bash
# Run all suites:
./eval/docker/run-all-containerized.sh all -j 4

# Run a single suite:
./eval/docker/run-all-containerized.sh atomic -j 8
./eval/docker/run-all-containerized.sh difficulty
./eval/docker/run-all-containerized.sh reruns -j 6
./eval/docker/run-all-containerized.sh multi-turn
```

### Four Evaluation Suites

| Suite | Design | Purpose | Sessions |
|-------|--------|---------|:--------:|
| **1. Atomic** | 20 isolated API prompts × 2 conditions | Test individual LASER API knowledge | 40 |
| **2. Reruns** | 5 prompts × 3 runs × 2 diseases × 2 conditions | Statistical reproducibility | 60 |
| **3. Difficulty** | 8 prompts of increasing complexity × 2 conditions | Map the skill advantage curve | 16 |
| **4. Multi-Turn** | 5 sequential build steps × 2 conditions × 3 runs | Test iterative model building | 30 |

### Scoring

- **Atomic**: 0-2 per prompt (API correctness)
- **Difficulty**: 0-3 per prompt (API correctness)
- **Reruns**: 0-3 per prompt on API Correctness (AC)
- **Multi-Turn**: 5 dimensions — AC, Structural Correctness, Domain Adaptation,
  Completeness, Consistency Bonus

### Polio Test Prompts (Reruns Suite)

| # | Prompt | Tests | Difficulty |
|---|--------|-------|:----------:|
| 1 | Basic 10-patch SEIR | LASER API, imports, component signatures | Low |
| 2 | Gravity network + monsoon forcing | `gravity()`, `ValuesMap`, spatial coupling | Medium |
| 3 | OPV vaccination (RI + SIA) | Custom components, OPV waning | Medium-High |
| 4 | Calibration framework | Loss functions, AFP surveillance, parameter sampling | High |
| 5 | Full 20-district integration | End-to-end model with all features | Highest |

## Token Efficiency

Output bytes as a proxy for output tokens consumed:

| Suite | WITH avg/prompt | WITHOUT avg/prompt | Ratio |
|-------|:---------------:|:------------------:|:-----:|
| Atomic (simple) | 3,542 bytes | 3,523 bytes | **1.00x** |
| Difficulty D1-D5 | 13,792 bytes | 5,268 bytes | **0.38x** |
| Reruns (polio) | 10,661 bytes | 6,450 bytes | **0.60x** |

On simple prompts, both conditions produce identical output. On complex tasks,
WITHOUT produces 35-62% less output — summaries and descriptions instead of
complete runnable code.

## Repository Structure

```
laser-polio-pakistan-eval/
├── README.md                              # This file
├── CLAUDE.md                              # Project context for Claude Code
│
├── .claude/skills/                        # The skill pipeline under test
│   ├── laser-spatial-disease-modeling/    #   LASER model building
│   ├── epi-model-parametrization/        #   Parameter research guide
│   └── modelops-calabaria/               #   Calibration workflow
│
├── eval/
│   ├── docker/
│   │   ├── Dockerfile                     # Multi-stage: with-laser, without-laser
│   │   └── run-all-containerized.sh       # Unified parallel runner
│   │
│   ├── atomic/                            # Suite 1: 20 isolated API prompts
│   │   ├── prompt-A{01-20}.txt
│   │   ├── rubric.md
│   │   └── outputs-containerized/{with,without}-skill/
│   │
│   ├── difficulty-curve/                  # Suite 3: 8 complexity levels
│   │   ├── prompt-D{1-8}.txt
│   │   ├── rubric.md
│   │   └── outputs-containerized/
│   │       ├── analysis.md                #   Scored results + token analysis
│   │       └── {with,without}-skill/
│   │
│   ├── reruns/                            # Suite 2: statistical re-runs
│   │   ├── analysis.md
│   │   └── outputs-containerized/
│   │       ├── polio/run-{1,2,3}/{with,without}-skill/
│   │       └── guinea-worm/run-{1,2,3}/{with,without}-skill/
│   │
│   ├── multi-turn/                        # Suite 4: iterative building
│   │   ├── run-multi-turn.py              #   5-step sequential driver
│   │   ├── rubric.md
│   │   ├── steps.md
│   │   └── outputs-containerized/
│   │       ├── analysis.md                #   Scored results + methodology notes
│   │       └── run-{1,2,3}/{with,without}-skill/
│   │
│   ├── prompt-{1-5}.txt                   # Polio prompts (used by reruns)
│   ├── rubric.md                          # Main 4-dimension rubric
│   ├── guinea-worm/                       # Guinea worm prompts + outputs
│   │
│   └── outputs/
│       ├── containerized-eval-summary.md  # ★ Full cross-suite analysis
│       ├── ab-test-report.md              # Original polio A/B report
│       └── deep-analysis-report.md        # API audit + failure taxonomy
│
├── scripts/                               # WITH-skill generated models
│   ├── polio_seir_10patch.py              #   Main simulation
│   ├── custom_components.py               #   Vaccination + importation
│   └── ...
│
└── reference/
    └── polio-model-requirements.md        # Ground truth domain specification
```

## Reports

| Report | Description |
|--------|-------------|
| [`eval/outputs/containerized-eval-summary.md`](eval/outputs/containerized-eval-summary.md) | **Full cross-suite analysis** — accuracy, tokens, timing across all 4 suites |
| [`eval/difficulty-curve/outputs-containerized/analysis.md`](eval/difficulty-curve/outputs-containerized/analysis.md) | Difficulty curve: flat +1 advantage, no collapse |
| [`eval/multi-turn/outputs-containerized/analysis.md`](eval/multi-turn/outputs-containerized/analysis.md) | Multi-turn: execution paradox, methodology improvements |
| [`eval/outputs/ab-test-report.md`](eval/outputs/ab-test-report.md) | Original polio A/B: 60/60 vs 25/60 |
| [`eval/guinea-worm/outputs/ab-test-report.md`](eval/guinea-worm/outputs/ab-test-report.md) | Guinea worm stress test: 57/60 vs 29/60 |

## Running the Evaluation

### Prerequisites

- Docker Desktop running
- Claude Code CLI with valid OAuth credentials
- Generate auth token: `claude setup-token`, save to `~/.claude/oauth-token`

### Quick Start

```bash
# Build Docker images and run all 4 suites:
./eval/docker/run-all-containerized.sh all -j 4

# Run individual suites:
./eval/docker/run-all-containerized.sh atomic -j 8
./eval/docker/run-all-containerized.sh reruns -j 6
./eval/docker/run-all-containerized.sh multi-turn

# Run the original (non-containerized) eval:
./eval/run-eval.sh all
```

### Verify Results

```bash
# Check all output files are non-empty:
find eval/*/outputs-containerized -name "*.md" -empty

# Confirm no laser imports in WITHOUT outputs:
grep -r 'laser\.' eval/*/outputs-containerized/without-skill/ || echo 'Clean'

# Score using rubrics in each suite directory
```

## The LASER Framework

[LASER](https://laser.idmod.org/) (Light Agent Spatial modeling for
ERadication) is an open-source agent-based disease modeling framework from the
Institute for Disease Modeling (IDM). It provides:

- **Agent-based SEIR dynamics** with per-agent state tracking
- **Gravity-model spatial coupling** between geographic patches
- **Seasonal forcing** via `ValuesMap` time series
- **Vaccination components** (routine immunization, custom campaigns)
- **Vital dynamics** (births by CBR, mortality by CDR/estimator)

Install: `pip install laser-generic` (includes `laser-core` and `laser-generic`)

Documentation: [laser.idmod.org](https://laser.idmod.org/laser-generic/)

## Project Context

This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation's Global Health division. The work evaluates AI coding skills
for spatial disease transmission modeling using the LASER framework, applied to
polio eradication in Pakistan and guinea worm eradication in Chad. All modeling
uses standard compartmental frameworks (SEIR, SEIS) — purely computational
epidemiology for public health planning.

## License

This evaluation framework and analysis are provided for research and educational
purposes.
