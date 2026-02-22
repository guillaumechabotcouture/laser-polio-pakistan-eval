# Multi-Turn Analysis — Containerized Results (3 runs × 2 conditions)

## Methodology

- 5 sequential steps per condition, each building on the previous
- 2 conditions: WITH skill (laser-generic installed) vs WITHOUT skill (clean container)
- 3 independent runs per condition (30 total steps)
- Improved methodology vs initial run: 1200s timeout, retry on timeout, tighter
  step 1 prompt ("do NOT add gravity yet"), inline code instruction, resume support
- Scored on 5 dimensions: AC (0-3), SC (0-3), DA (0-3), CO (0-3), CB (0-1)

## Execution Summary

| | Run 1 | Run 2 | Run 3 | Total |
|---|:---:|:---:|:---:|:---:|
| **WITH running** | 1/5 | 2/5 | 3/5 | **6/15** |
| **WITHOUT running** | 0/5 | 0/5 | 3/5 | **3/15** |

Note: Most WITH "failures" are `plt.show()` timeouts (120s code execution
limit), not broken code. Run 2 step 4 discovered `matplotlib.use("Agg")`
fix, which resolved the issue for steps 4-5.

---

## Per-Step Scoring

### Run 1 — WITH Skill (63/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 3 | 3 | 3 | 3 | -- | 12/12 | timeout (plt.show) | All LASER imports correct; proper component order |
| 2 | 2 | 3 | 3 | 3 | 1 | 12/13 | timeout | Manual gravity (no `gravity()` import); correct network |
| 3 | 3 | 3 | 3 | 3 | 1 | 13/13 | YES | `ValuesMap.from_timeseries` correct |
| 4 | 3 | 3 | 3 | 3 | 1 | 13/13 | timeout | Correctly diagnoses symmetric I.C. → seeds patch 0 only |
| 5 | 3 | 3 | 3 | 3 | 1 | 13/13 | timeout | BirthsByCBR/MortalityByCDR correct; calc_capacity used |
| **Sum** | **14** | **15** | **15** | **15** | **4** | **63/64** | **1/5** | |

### Run 1 — WITHOUT Skill (43/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 0 | 2 | 2 | 3 | -- | 7/12 | FAILED | `from laser_core import LaserFrame` — wrong module; reimplements SEIR |
| 2 | 0 | 2 | 3 | 3 | 1 | 9/13 | FAILED | Same wrong import; adds gravity manually (correct math) |
| 3 | 0 | 2 | 3 | 3 | 1 | 9/13 | FAILED | Reinvents ValuesMap class; correct seasonal profile |
| 4 | 0 | 2 | 3 | 3 | 1 | 9/13 | timeout | Stubs PropertySet/LaserFrame to fix ImportError |
| 5 | 0 | 2 | 3 | 3 | 1 | 9/13 | timeout | Reinvents BirthsByCBR/MortalityByCDR; CBR=30/CDR=10 correct |
| **Sum** | **0** | **10** | **14** | **15** | **4** | **43/64** | **0/5** | |

### Run 2 — WITH Skill (63/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 3 | 3 | 3 | 3 | -- | 12/12 | timeout | All LASER imports correct |
| 2 | 2 | 3 | 3 | 3 | 1 | 12/13 | timeout | Manual gravity; correct network |
| 3 | 3 | 3 | 3 | 3 | 1 | 13/13 | timeout | ValuesMap.from_timeseries correct |
| 4 | 3 | 3 | 3 | 3 | 1 | 13/13 | YES | `matplotlib.use("Agg")` fix; correct diagnosis |
| 5 | 3 | 3 | 3 | 3 | 1 | 13/13 | YES | BirthsByCBR/MortalityByCDR correct; growth verified 1.97% |
| **Sum** | **14** | **15** | **15** | **15** | **4** | **63/64** | **2/5** | |

### Run 2 — WITHOUT Skill (23/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 2 | 2 | 2 | 3 | -- | 9/12 | timeout | `laser.core.PropertySet` + hallucinated `grid()`, `initialize_population()` |
| 2 | 0 | 0 | 0 | 0 | 0 | 0/13 | FAILED | Code extraction failure — only got a visualization fragment |
| 3 | 0 | 1 | 1 | 2 | 0 | 4/13 | timeout | Rewrites from scratch as numpy Euler; wrong pops, beta, gravity |
| 4 | 0 | 1 | 1 | 2 | 1 | 5/13 | timeout | Builds on step 3's wrong params |
| 5 | 0 | 1 | 1 | 2 | 1 | 5/13 | timeout | Reinvents vital dynamics; correct units but wrong base params |
| **Sum** | **2** | **5** | **5** | **9** | **2** | **23/64** | **0/5** | |

### Run 3 — WITH Skill (62/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 3 | 3 | 2 | 3 | -- | 11/12 | YES | `dists.constant_int(4)` instead of gamma (minor) |
| 2 | 2 | 3 | 3 | 3 | 1 | 12/13 | timeout | Manual gravity; correct network |
| 3 | 3 | 3 | 3 | 3 | 1 | 13/13 | YES | ValuesMap.from_timeseries; smooth cos transitions |
| 4 | 3 | 3 | 3 | 3 | 1 | 13/13 | YES | Seeds patch 0; confirms spatial propagation |
| 5 | 3 | 3 | 3 | 3 | 1 | 13/13 | timeout | BirthsByCBR/MortalityByCDR correct; calc_capacity |
| **Sum** | **14** | **15** | **14** | **15** | **4** | **62/64** | **3/5** | |

### Run 3 — WITHOUT Skill (19/64)

| Step | AC | SC | DA | CO | CB | Total | Runs? | Notes |
|------|:--:|:--:|:--:|:--:|:---:|:-----:|:-----:|-------|
| 1 | 0 | 0 | 0 | 0 | -- | 0/12 | NO CODE | No Python code extracted from response |
| 2 | 0 | 1 | 1 | 2 | 0 | 4/13 | YES | Rewrites as numpy SIR (no E!); wrong pops/beta |
| 3 | 0 | 1 | 1 | 2 | 1 | 5/13 | FAILED | Builds on step 2; syntax error |
| 4 | 0 | 1 | 1 | 2 | 1 | 5/13 | YES | Fixes syntax error; seasonal forcing works |
| 5 | 0 | 1 | 1 | 2 | 1 | 5/13 | YES | Vital dynamics stubs; CBR/CDR correct; wrong base params |
| **Sum** | **0** | **4** | **4** | **8** | **3** | **19/64** | **3/5** | |

---

## Cross-Run Aggregates

### Totals per condition

| Run | WITH | WITHOUT | Delta |
|-----|:----:|:-------:|:-----:|
| 1 | 63/64 | 43/64 | +20 |
| 2 | 63/64 | 23/64 | +40 |
| 3 | 62/64 | 19/64 | +43 |
| **Mean** | **62.7** | **28.3** | **+34.3** |

### By dimension (mean ± std)

| Dimension | Max | WITH | WITHOUT | Delta |
|-----------|:---:|:----:|:-------:|:-----:|
| API Correctness | 15 | 14.0 ± 0.0 | 0.7 ± 0.9 | **+13.3** |
| Structural Correctness | 15 | 15.0 ± 0.0 | 6.3 ± 2.6 | **+8.7** |
| Domain Adaptation | 15 | 14.7 ± 0.5 | 7.7 ± 4.5 | **+7.0** |
| Completeness | 15 | 15.0 ± 0.0 | 10.7 ± 3.1 | **+4.3** |
| Consistency Bonus | 4 | 4.0 ± 0.0 | 3.0 ± 0.8 | **+1.0** |
| **Total** | **64** | **62.7 ± 0.5** | **28.3 ± 10.4** | **+34.3** |
| Steps running | 5 | 2.0 | 1.0 | +1.0 |

### Variance comparison

- **WITH std = 0.5** (scores: 63, 63, 62) — near-zero variance
- **WITHOUT std = 10.4** (scores: 43, 23, 19) — **22x higher variance**

The skill produces rock-solid consistency. WITHOUT is a lottery.

---

## Key Findings

### 1. The methodology improvements resolved the execution paradox

The initial run showed WITH 1/5 running vs WITHOUT 3/5 — suggesting the skill
*hurt* executability. The 3-run replication shows this was a fluke:

| | Initial run | 3-run mean |
|---|:---:|:---:|
| WITH running | 1/5 | 2.0/5 |
| WITHOUT running | 3/5 | 1.0/5 |

The tighter step 1 prompt and inline code instruction fixed WITH's
over-engineering problem. The `plt.show()` timeout is the remaining issue —
Run 2 discovered the `matplotlib.use("Agg")` fix mid-run.

### 2. API Correctness is the decisive dimension

AC accounts for 13.3 of the 34.3-point gap (39%). WITH scores AC=14/15 in
every run (losing 1 point for manual gravity instead of `gravity()` import).
WITHOUT scores AC=0 in 2/3 runs, AC=2 in 1 run. The skill's value is
overwhelmingly about LASER API knowledge.

### 3. WITHOUT-skill variance is 22x higher

WITH totals: {63, 63, 62}, std=0.5. WITHOUT totals: {43, 23, 19}, std=10.4.
The skill eliminates run-to-run variance. WITHOUT outcomes depend on whether
the initial LASER import guess succeeds (Run 1 vs Run 2-3).

### 4. WITHOUT Run 1 was an outlier — a coherent reimplementation

Run 1 WITHOUT scored 43/64 because it maintained a consistent, well-structured
reimplementation across all 5 steps. Although it never used real LASER APIs
(AC=0), it got domain parameters right (DA=14/15), produced complete code
(CO=15/15), and maintained perfect consistency (CB=4/4). This was the best
possible WITHOUT outcome — and it still scored 20 points below WITH.

### 5. Cascading failures in WITHOUT runs 2-3

- **Run 2**: Code extraction failure at step 2 (only got a visualization
  fragment) forced a complete rewrite at step 3 with wrong parameters that
  persisted through steps 4-5.
- **Run 3**: Step 1 produced no extractable code. Step 2 started from scratch
  with a simplified SIR (no E compartment!) and wrong populations.

### 6. WITH consistently diagnoses the "identical curves" problem correctly

All 3 WITH runs correctly identified that symmetric initial conditions cause
identical curves across patches, and fixed it by seeding infection in only
patch 0. This is a sophisticated epidemiological insight that demonstrates
domain understanding, not just API knowledge.

### 7. The `gravity()` import is the one consistent WITH weakness

None of the 3 WITH runs used `from laser.core.migration import gravity`. All
computed the gravity network manually using numpy. This is mathematically
correct but loses 1 AC point per run. The skill's gravity documentation may
not strongly enough push toward the specific API function.

---

## Comparison: Initial Run vs 3-Run Replication

| Metric | Initial (n=1) | 3-Run (n=3 mean) | Change |
|--------|:---:|:---:|---|
| WITH total | 40/64 | 62.7/64 | **+22.7** — methodology fixes helped enormously |
| WITHOUT total | 24/64 | 28.3/64 | +4.3 — modest improvement |
| WITH running | 1/5 | 2.0/5 | +1.0 |
| WITHOUT running | 3/5 | 1.0/5 | -2.0 — initial run was the fluke |
| Delta | +16 | +34.3 | Gap widened with better methodology |

The initial run's "execution paradox" (WITHOUT running more code) was an
artifact of:
1. WITH over-engineering step 1 (fixed by tighter prompt)
2. Short timeouts killing correct WITH code (fixed by 1200s timeout)
3. `plt.show()` blocking (identified in run 2, fixable with `Agg` backend)
4. n=1 variance (WITHOUT Run 1 was an outlier best-case)

---

## Statistical Significance

Paired t-test on total scores (WITH minus WITHOUT):

| Run | WITH | WITHOUT | d_i |
|-----|------|---------|-----|
| 1 | 63 | 43 | +20 |
| 2 | 63 | 23 | +40 |
| 3 | 62 | 19 | +43 |

- d_bar = 34.3, s_d = 10.2, SE = 5.9
- t(2) = 5.83, p = 0.028
- 95% CI: [8.9, 59.8]
- Cohen's d = 3.4

**p < 0.05 with n=3.** The skill advantage on multi-turn model building is
statistically significant even with small sample size, driven by the massive
and consistent AC gap.
