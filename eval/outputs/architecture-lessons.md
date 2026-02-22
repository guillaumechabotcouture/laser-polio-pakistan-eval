# Lessons from Eval for the Disease Modeling Plugin Architecture

This document maps findings from the containerized A/B evaluation (4 suites,
146 sessions) back to the proposed "Disease Modeling Plugin: Technical
Architecture" draft. It identifies what the original architecture got right,
what should change based on evidence, and proposes a revised priority ordering.

---

## What the Architecture Got Right

### 1. "Never build a disease transmission solver from scratch"

The Model Selector Skill's core principle is validated with hard numbers.
Without the skill, Claude abandons LASER and builds numpy solvers in 100% of
cases beyond the simplest prompt. This is not a rare failure — it is the
deterministic, reproducible outcome across all 3 statistical runs.

- Polio: LASER abandoned at prompt 2 (gravity + seasonality) in 3/3 runs
- Guinea worm: LASER abandoned at prompt 3 (custom components) in 3/3 runs
- Multi-turn: LASER abandoned at step 1-2 in 3/3 runs

The architecture's instinct to route users to existing frameworks rather than
letting AI reinvent is exactly correct.

### 2. Skills as the knowledge delivery mechanism

The architecture bets on Skills (SKILL.md + references/) as the right vehicle
for domain knowledge. The eval confirms this decisively:

| Metric | WITH Skill | WITHOUT Skill |
|--------|:----------:|:------------:|
| Atomic API correctness | 97.5% | 67.5% |
| Full model building (reruns) | 100% | 11-20% |
| Multi-turn quality (/64) | 62.7 ± 0.5 | 28.3 ± 10.4 |
| Run-to-run variance | Near zero | 22x higher |

Skills are the single most impactful component in the entire proposed
architecture.

### 3. Cross-framework skill applicability

The architecture proposes one skill per framework, not one per disease. The
guinea worm eval validates this — the LASER skill (documented with
measles/polio examples) correctly guided guinea worm models with SEIS dynamics,
low R0, no vaccination, and dry-season forcing. Zero negative transfer.

This means framework-specific skills can serve the full disease portfolio
without per-disease customization.

### 4. Progressive disclosure model

The architecture describes Skills loading only name/description at startup,
full content on demand. The eval confirms this doesn't waste tokens on simple
tasks (identical output volume on atomic prompts) but enables dramatically
more complete output on complex ones (+40-62% more code). The progressive
disclosure hypothesis is confirmed.

### 5. The ToolUniverse pattern

The principle of "unified discovery layer over heterogeneous backends" is
sound. The Skills layer maps to ToolUniverse's "Tool Finder" — lightweight
routing to the right framework, with deep domain knowledge loaded on demand.

---

## What Should Change Based on Evidence

### 1. Skills are the highest-ROI component — reorder build priorities

The architecture puts the LASER Guide Skill at Phase 2 ("Low-Medium effort")
and the Modeling Gateway at Phase 1 ("Medium-High effort, Very High impact").
The eval shows this is backwards:

- The LASER skill alone (no gateway, no MCP, no agents) produces 97-100% API
  correctness on full model building tasks
- The gateway provides execution infrastructure, which is valuable but additive
- Skills provide the knowledge foundation that everything else depends on

**Revised priority matrix:**

| Component | Architecture Priority | Recommended Priority | Evidence |
|-----------|:---:|:---:|---|
| Framework Skills (LASER, Starsim) | Phase 1-2 | **Phase 1, highest** | +80pp on model building |
| Evaluation Framework | Phase 2 | **Phase 1** | Built and proven (4 suites) |
| Model Selector Skill | Phase 1, "linchpin" | Phase 1, lower | Modelers know their framework |
| Calibration Skill + calabaria | Phase 2-3 | **Phase 1-2** | WHERE=0 on calibration without skill |
| Modeling Gateway MCP | Phase 1 | **Phase 2** | Skills alone provide most value |
| Commands / Agents | Phase 1-2 | Phase 2 | Additive, not foundational |

### 2. Skill internal structure matters enormously

The architecture treats Skills as "SKILL.md + references/ + scripts/" without
specifying internal architecture. The eval shows that structure drives quality:

- The tested skill uses a 3-layer architecture:
  - **Layer 1 — Discipline**: What to do, what NOT to do, critical gotchas
  - **Layer 2 — Scaffolding**: Step-by-step workflow with code templates
  - **Layer 3 — Verification**: How to validate outputs, common failure modes
- This produces near-zero variance (std=0.5 out of 64 across 3 multi-turn runs)
- A poorly structured skill would likely produce much weaker results

**Recommendation:** Add a "Skill Architecture Guide" to the architecture
document specifying:

```
framework-guide/
├── SKILL.md              # Layer 1: Discipline + workflow overview
│                         #   - "ALWAYS use X, NEVER do Y"
│                         #   - Critical gotchas (units, state semantics)
│                         #   - Step-by-step modeling workflow
├── references/
│   ├── api-reference.md  # Layer 2: Complete API documentation
│   └── patterns.md       #   Common patterns, templates, examples
└── scripts/
    ├── components.py     # Layer 3: Reusable verified components
    └── validation.py     #   Verification/validation utilities
```

Key design principles from the eval:
- **Gotchas prevent silent failures.** The skill documents that `BirthsByCBR`
  expects per-1000/year rates (not daily per-capita) and that vaccination must
  set `state = RECOVERED` (not `susceptibility`). These prevent errors that
  produce plausible but wrong results with no runtime errors.
- **Code templates > prose descriptions.** The skill provides working code
  snippets for each API. WITHOUT the skill, Claude describes APIs in prose
  but doesn't produce complete runnable code (62% less output on complex tasks).
- **Reference files should be exhaustive.** The API reference covers every
  class, method, and parameter. This is what enables 97-100% API correctness.

### 3. The Model Selector Skill is less critical than assumed

The architecture calls it "the linchpin component." The eval suggests
framework-specific Skills are the actual linchpin.

**Why:** The bottleneck is framework-specific API knowledge, not framework
selection. A disease modeler already knows they want LASER or Starsim — they
need help using the API correctly, not choosing which framework.

**Where Model Selector IS valuable:** Non-expert users (grantees, new team
members) who don't know the ecosystem. But for the core workflow, framework
Skills are the linchpin.

### 4. The eval methodology is a solved problem

Section 6.4 of the architecture asks Anthropic for "guidance on how to evaluate
whether Skills improve model quality." This is now a solved problem:

**What we built:**
- Docker-isolated A/B testing (prevents WITHOUT reading source code)
- 4 complementary suites:
  - Atomic (20 isolated API prompts) — tests individual knowledge
  - Difficulty Curve (8 levels) — maps the advantage curve
  - Statistical Re-Runs (3 runs × 2 diseases) — tests reproducibility
  - Multi-Turn (5 sequential steps × 3 runs) — tests iterative building
- Rubrics with 1-5 scoring dimensions
- Statistical significance testing (paired t-tests, Cohen's d)
- Resume support, parallel execution, automated scoring

**Recommendation:** Share this with Anthropic as a deliverable rather than
requesting their methodology. It demonstrates rigorous skill evaluation at a
level they haven't publicly documented.

### 5. The "token hungry MCP" concern is backwards

The architecture worries about multiple MCP servers consuming context ("clunky
and token hungry"). The eval shows Skills make Claude MORE token-efficient:

| Suite | WITH output | WITHOUT output | Ratio |
|-------|:----------:|:-------------:|:-----:|
| Atomic (simple) | 3,542 bytes | 3,523 bytes | 1.00x |
| Difficulty D1-D5 | 13,792 bytes | 5,268 bytes | **2.62x** |
| Reruns (polio) | 10,661 bytes | 6,450 bytes | **1.65x** |

WITH-skill produces MORE output (complete code) in LESS time (30% faster on
multi-turn). The skill enables confident code generation rather than hedging
with summaries. The gateway unification for token efficiency is solving a
non-problem.

### 6. The architecture underestimates the execution environment problem

The multi-turn eval revealed that correct LASER code can fail to run due to
environment issues:
- `plt.show()` blocking in non-interactive mode
- Path issues in containerized environments
- Missing data files or dependencies

The architecture's Compute Dispatch Skill is described as "Low effort, 1 day."
But execution environment setup is more critical than that. The skill needs to
guide not just what code to write, but:
- How to configure matplotlib for non-interactive use
- How to handle data file paths across environments
- How to verify the environment has required packages

### 7. Calibration needs Skills + MCP together

The architecture has "Calibration Skill" and "calabaria integration" as
separate Phase 2-3 items. The eval shows that calibration (reruns prompt 4) is
where WITHOUT completely falls apart — AC=0 across all runs for both diseases.

Calibration is where Skills + MCP together provide the most value:
- Skills for the calibration pattern (parameter ranges, loss functions, sampling)
- MCP for execution (submitting sweeps, retrieving results)

These should be developed together, not sequentially.

### 8. Quantify the "What Claude Already Knows" baseline

The eval establishes exactly what Claude knows about LASER without any help:

| API Area | Claude's Baseline | Skill Needed? |
|----------|:---:|:---:|
| Basic SEIR concepts | Good (AC=2) | Minimal |
| Import paths | Unreliable (`laser_core` vs `laser.core`) | **Yes** |
| `gravity()` 6-arg signature | Poor (reimplements manually) | **Critical** |
| `row_normalizer()` | Partial (wrong param name) | **Yes** |
| `calc_capacity()` | Poor (wrong signature) | **Critical** |
| Vital dynamics units | Partial (values right, signatures wrong) | **Yes** |
| `ValuesMap.from_timeseries` | Moderate | Yes |
| Custom component protocol | Good (`__init__`/`step` is intuitive) | Minimal |
| State enum values | Good (knows S=0, E=1, I=2, R=3) | Minimal |
| Vaccination gotcha | Good (when prompted about it) | Yes (unprompted) |

This baseline should inform which APIs get the most coverage in each
framework-specific skill. Focus depth on the weakest areas (spatial coupling,
vital dynamics signatures) rather than treating all APIs equally.

---

## Revised Architecture Recommendations

### Phase 1: Skills + Eval (8-12 weeks)

1. **Framework-specific Skills** — Build for LASER and Starsim using the
   3-layer architecture (discipline/scaffolding/verification). These provide
   the largest measured impact (+80pp on model building).

2. **Evaluation Framework** — Already built. Formalize as a reusable tool:
   Docker-isolated A/B, 4 suite types, automated scoring. Use to validate
   each new skill before release.

3. **Model Selector Skill** — Lightweight routing skill. Less critical than
   framework skills but easy to build.

4. **Calibration Skill** — High-impact for the full modeling workflow.
   Document calabaria patterns, parameter sampling, loss function design.

### Phase 2: Execution + Integration (8-12 weeks)

5. **Modeling Gateway MCP** — Now that Skills provide the knowledge layer,
   add execution infrastructure. Start with single-run dispatch, add sweep
   support.

6. **Compute Dispatch Skill** — Include execution environment guidance
   (matplotlib backend, path handling, dependency verification).

7. **Commands** — `/new-model`, `/calibrate`, `/check-job` built on top of
   working Skills + MCP.

### Phase 3: Ecosystem + Scale

8. **Agents** — model-architect, data-scout, calibration-reviewer. Additive
   value on top of the Skills + MCP foundation.

9. **Additional framework skills** — CMS, calabaria advanced workflows.

10. **ToolUniverse integration** — Expose IDM tools to the broader ecosystem.

---

## What to Share with Anthropic

The architecture document proposes several asks of Anthropic (Sections 6.1-6.6).
Based on the eval results, revise these:

| Original Ask | Status | Revised Approach |
|---|---|---|
| 6.1 Architecture review | Still valuable | Share eval results to demonstrate our approach |
| 6.2 CC Web environment support | Still needed | Technical integration question |
| 6.3 Skills design feedback | **We have data now** | Share the 3-layer architecture + eval proving it works |
| 6.4 Evaluation methodology | **Solved** | Share our 4-suite Docker-isolated framework as a deliverable |
| 6.5 Case study visibility | Stronger with data | "97-100% API correctness, +80pp on model building" is a compelling story |
| 6.6 Co-development partnership | Stronger with data | We have a proven, generalizable pattern for scientific domain skills |

The eval transforms the Anthropic conversation from "we need help designing
this" to "we built this, here's what works, let's scale it together."
