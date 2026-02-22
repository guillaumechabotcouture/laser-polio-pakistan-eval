# Suite 2: Statistical Re-Runs — Scoring Protocol and Analysis

## Purpose

Establish score distributions and confidence intervals for the LASER skill A/B
tests by re-running each evaluation 3 times. The original A/B tests (N=1)
produced:

- **Polio:** WITH 60/60 vs WITHOUT 25/60
- **Guinea worm:** WITH 57/60 vs WITHOUT 29/60

This suite answers: Are these results reproducible? What is the variance?

## Design

- 2 diseases (polio, guinea worm)
- 5 prompts per disease (same as originals)
- 2 conditions (with-skill, without-skill)
- 3 independent runs
- **Total: 60 scored sessions**

## Scoring Protocol

### Rubrics

Use the same rubrics as the original evaluations:

- **Polio:** `eval/rubric.md`
- **Guinea worm:** `eval/guinea-worm/rubric.md`

### Dimensions (0-3 each, max 12 per session)

| Dimension | Abbrev | What it measures |
|-----------|--------|------------------|
| API Correctness | AC | Correct use of LASER v1.0.0 imports, signatures, patterns |
| Structural Correctness | SC | Working model assembly (GeoDataFrame, components, run call) |
| Domain Adaptation | DA | Disease-specific parameter and biology adaptation |
| Completeness | CO | Runnable code vs pseudocode/fragments |

### Scoring Rules

1. Score each of the 60 outputs independently using the appropriate rubric.
2. Do NOT let knowledge of which run a session belongs to bias scoring.
3. Score all sessions for one disease before moving to the other.
4. Recommended order: score all with-skill first, then all without-skill
   (reduces anchoring from condition comparison within a prompt).

---

## Score Sheets

### Polio — Run 1

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

**Run 1 totals:** WITH ___/60 | WITHOUT ___/60

### Polio — Run 2

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

**Run 2 totals:** WITH ___/60 | WITHOUT ___/60

### Polio — Run 3

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

**Run 3 totals:** WITH ___/60 | WITHOUT ___/60

### Guinea Worm — Run 1

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

**Run 1 totals:** WITH ___/60 | WITHOUT ___/60

### Guinea Worm — Run 2

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

**Run 2 totals:** WITH ___/60 | WITHOUT ___/60

### Guinea Worm — Run 3

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

**Run 3 totals:** WITH ___/60 | WITHOUT ___/60

---

## Analysis Templates

### 1. Headline Scores Across Runs

#### Polio

| Run | WITH (/60) | WITHOUT (/60) | Difference | Delta % |
|-----|-----------|--------------|------------|---------|
| Original | 60 | 25 | +35 | +140% |
| Run 1 | | | | |
| Run 2 | | | | |
| Run 3 | | | | |
| **Mean** | | | | |
| **Std Dev** | | | | |
| **95% CI** | | | | |

#### Guinea Worm

| Run | WITH (/60) | WITHOUT (/60) | Difference | Delta % |
|-----|-----------|--------------|------------|---------|
| Original | 57 | 29 | +28 | +97% |
| Run 1 | | | | |
| Run 2 | | | | |
| Run 3 | | | | |
| **Mean** | | | | |
| **Std Dev** | | | | |
| **95% CI** | | | | |

### 2. Per-Prompt Score Distribution (WITH skill)

#### Polio WITH skill

| Prompt | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Std |
|--------|-------|-------|-------|-----|------|-----|-----|
| 1 | /12 | /12 | /12 | | | | |
| 2 | /12 | /12 | /12 | | | | |
| 3 | /12 | /12 | /12 | | | | |
| 4 | /12 | /12 | /12 | | | | |
| 5 | /12 | /12 | /12 | | | | |

#### Polio WITHOUT skill

| Prompt | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Std |
|--------|-------|-------|-------|-----|------|-----|-----|
| 1 | /12 | /12 | /12 | | | | |
| 2 | /12 | /12 | /12 | | | | |
| 3 | /12 | /12 | /12 | | | | |
| 4 | /12 | /12 | /12 | | | | |
| 5 | /12 | /12 | /12 | | | | |

#### Guinea Worm WITH skill

| Prompt | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Std |
|--------|-------|-------|-------|-----|------|-----|-----|
| 1 | /12 | /12 | /12 | | | | |
| 2 | /12 | /12 | /12 | | | | |
| 3 | /12 | /12 | /12 | | | | |
| 4 | /12 | /12 | /12 | | | | |
| 5 | /12 | /12 | /12 | | | | |

#### Guinea Worm WITHOUT skill

| Prompt | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Std |
|--------|-------|-------|-------|-----|------|-----|-----|
| 1 | /12 | /12 | /12 | | | | |
| 2 | /12 | /12 | /12 | | | | |
| 3 | /12 | /12 | /12 | | | | |
| 4 | /12 | /12 | /12 | | | | |
| 5 | /12 | /12 | /12 | | | | |

### 3. Per-Dimension Variance Analysis

Are certain dimensions stable across runs or highly variable?

#### Polio — WITH skill (sum across 5 prompts, /15 per dimension)

| Dimension | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Range |
|-----------|-------|-------|-------|-----|------|-----|-------|
| AC | /15 | /15 | /15 | | | | |
| SC | /15 | /15 | /15 | | | | |
| DA | /15 | /15 | /15 | | | | |
| CO | /15 | /15 | /15 | | | | |

#### Polio — WITHOUT skill (sum across 5 prompts, /15 per dimension)

| Dimension | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Range |
|-----------|-------|-------|-------|-----|------|-----|-------|
| AC | /15 | /15 | /15 | | | | |
| SC | /15 | /15 | /15 | | | | |
| DA | /15 | /15 | /15 | | | | |
| CO | /15 | /15 | /15 | | | | |

#### Guinea Worm — WITH skill (sum across 5 prompts, /15 per dimension)

| Dimension | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Range |
|-----------|-------|-------|-------|-----|------|-----|-------|
| AC | /15 | /15 | /15 | | | | |
| SC | /15 | /15 | /15 | | | | |
| DA | /15 | /15 | /15 | | | | |
| CO | /15 | /15 | /15 | | | | |

#### Guinea Worm — WITHOUT skill (sum across 5 prompts, /15 per dimension)

| Dimension | Run 1 | Run 2 | Run 3 | Min | Mean | Max | Range |
|-----------|-------|-------|-------|-----|------|-----|-------|
| AC | /15 | /15 | /15 | | | | |
| SC | /15 | /15 | /15 | | | | |
| DA | /15 | /15 | /15 | | | | |
| CO | /15 | /15 | /15 | | | | |

### 4. Framework Abandonment Consistency

In the original tests, WITHOUT-skill sessions tended to abandon LASER at later
prompts (which require more specific API knowledge). Does this pattern hold
across re-runs?

**Definition:** "Framework abandonment" = the session either (a) uses an
entirely different framework (e.g., custom Python, Mesa, NetworkX), or (b)
writes pseudocode that references LASER conceptually but does not use actual
LASER API calls.

#### Polio — WITHOUT skill: Framework abandonment by prompt

| Prompt | Run 1 | Run 2 | Run 3 | Abandoned in N/3 runs |
|--------|-------|-------|-------|----------------------|
| 1 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 2 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 3 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 4 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 5 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |

**Threshold prompt:** The first prompt where abandonment occurs consistently
(2+ of 3 runs): ___

#### Guinea Worm — WITHOUT skill: Framework abandonment by prompt

| Prompt | Run 1 | Run 2 | Run 3 | Abandoned in N/3 runs |
|--------|-------|-------|-------|----------------------|
| 1 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 2 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 3 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 4 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |
| 5 | LASER / Abandoned | LASER / Abandoned | LASER / Abandoned | /3 |

**Threshold prompt:** ___

**WITH skill:** Record any framework abandonment in with-skill sessions.

| Disease | Run 1 | Run 2 | Run 3 |
|---------|-------|-------|-------|
| Polio | None / Prompt ___ | None / Prompt ___ | None / Prompt ___ |
| Guinea worm | None / Prompt ___ | None / Prompt ___ | None / Prompt ___ |

### 5. Cross-Run Correlation

Does the same prompt always produce similar scores across runs, or is there
high variance for specific prompts?

#### Polio: Per-prompt score variance (Total /12, across 3 runs)

| Prompt | Condition | Run 1 | Run 2 | Run 3 | Variance | Stable? |
|--------|-----------|-------|-------|-------|----------|---------|
| 1 | with | | | | | |
| 1 | without | | | | | |
| 2 | with | | | | | |
| 2 | without | | | | | |
| 3 | with | | | | | |
| 3 | without | | | | | |
| 4 | with | | | | | |
| 4 | without | | | | | |
| 5 | with | | | | | |
| 5 | without | | | | | |

#### Guinea Worm: Per-prompt score variance (Total /12, across 3 runs)

| Prompt | Condition | Run 1 | Run 2 | Run 3 | Variance | Stable? |
|--------|-----------|-------|-------|-------|----------|---------|
| 1 | with | | | | | |
| 1 | without | | | | | |
| 2 | with | | | | | |
| 2 | without | | | | | |
| 3 | with | | | | | |
| 3 | without | | | | | |
| 4 | with | | | | | |
| 4 | without | | | | | |
| 5 | with | | | | | |
| 5 | without | | | | | |

**Stability threshold:** variance <= 1.0 = Stable, variance > 2.0 = Unstable

### 6. Execution Testing Summary

| Disease | Condition | Run 1 pass | Run 2 pass | Run 3 pass | Pass rate |
|---------|-----------|-----------|-----------|-----------|-----------|
| Polio | with | /5 | /5 | /5 | /15 |
| Polio | without | /5 | /5 | /5 | /15 |
| Guinea worm | with | /5 | /5 | /5 | /15 |
| Guinea worm | without | /5 | /5 | /5 | /15 |

---

## Expected Outputs

### Reproducibility of N=1 Findings

The original single-run results were:

| Disease | WITH | WITHOUT | Difference |
|---------|------|---------|-----------|
| Polio | 60/60 | 25/60 | +35 |
| Guinea worm | 57/60 | 29/60 | +28 |

**Key questions this suite answers:**

1. **Overall mean and std:** What is the mean WITH and WITHOUT score across 3
   runs? Is the standard deviation small (< 3 points) or large (> 5 points)?

2. **Per-dimension stability:** Is AC always the strongest skill advantage, or
   does it vary? Example: AC might be stable at +2.4 across runs while DA
   might range from +0.0 to +1.5.

3. **Framework abandonment threshold:** In the original tests, WITHOUT
   condition tended to abandon LASER at specific prompt complexity levels.
   Does this happen at the same prompt number each time, or is it variable?

4. **N=1 reproducibility:** Can we say with 95% confidence that WITH > WITHOUT
   by at least X points? With N=3 runs (each /60), a paired t-test or
   bootstrap CI will provide this.

### Statistical Methods

For N=3 paired observations (WITH_i, WITHOUT_i) for i = 1, 2, 3:

1. **Paired differences:** d_i = WITH_i - WITHOUT_i
2. **Mean difference:** d_bar = mean(d_1, d_2, d_3)
3. **Std of differences:** s_d = std(d_1, d_2, d_3)
4. **95% CI:** d_bar +/- t_{0.025, df=2} * s_d / sqrt(3)
   where t_{0.025, df=2} = 4.303
5. **Effect size:** Cohen's d = d_bar / s_d

If CI lower bound > 0, the skill advantage is statistically significant at
p < 0.05 even with only N=3.

### Confidence Interval Computation

#### Polio

| Metric | Value |
|--------|-------|
| d_1 (Run 1: WITH - WITHOUT) | |
| d_2 (Run 2: WITH - WITHOUT) | |
| d_3 (Run 3: WITH - WITHOUT) | |
| d_bar (mean difference) | |
| s_d (std of differences) | |
| SE (s_d / sqrt(3)) | |
| 95% CI lower | |
| 95% CI upper | |
| Cohen's d | |
| Significant? (CI_lower > 0) | |

#### Guinea Worm

| Metric | Value |
|--------|-------|
| d_1 (Run 1: WITH - WITHOUT) | |
| d_2 (Run 2: WITH - WITHOUT) | |
| d_3 (Run 3: WITH - WITHOUT) | |
| d_bar (mean difference) | |
| s_d (std of differences) | |
| SE (s_d / sqrt(3)) | |
| 95% CI lower | |
| 95% CI upper | |
| Cohen's d | |
| Significant? (CI_lower > 0) | |

#### Combined (pooled across diseases)

| Metric | Value |
|--------|-------|
| N observations | 6 |
| d_bar (mean difference) | |
| s_d (std of differences) | |
| SE (s_d / sqrt(6)) | |
| 95% CI lower | |
| 95% CI upper | |
| Cohen's d | |
| Significant? | |

---

## Summary Findings

*(Fill in after scoring all 60 sessions)*

### Headline Result

> The LASER skill advantage of +___ to +___ points (/60) was observed across
> 3 independent runs for [polio/guinea worm/both], with a 95% CI of
> [___, ___]. The original N=1 findings [are/are not] reproducible.

### Dimension Stability

| Dimension | Most stable or most variable? | Notes |
|-----------|-------------------------------|-------|
| AC | | |
| SC | | |
| DA | | |
| CO | | |

### Framework Abandonment Pattern

> WITHOUT-skill sessions consistently abandoned LASER at prompt ___ for polio
> and prompt ___ for guinea worm. The abandonment threshold was
> [consistent/variable] across runs.

### Limitations

- N=3 runs gives limited statistical power (wide CIs expected).
- Scorer subjectivity: a single scorer introduces potential bias. Consider
  having a second scorer independently score a subset for inter-rater
  reliability.
- LLM temperature/randomness: sessions on different days may encounter
  different model versions or server-side configuration changes.
- The "resume support" feature means partially completed runs may span
  different model snapshots if run across multiple days.
