# Multi-Turn Analysis — Containerized Results

## Execution Summary

| Step | WITH Skill | | WITHOUT Skill | |
|------|-----------|---|--------------|---|
| | Time | Runs? | Time | Runs? |
| 1 | 415s | **NO** (14,597 chars) | 600s | **NO** (timeout, 0 chars) |
| 2 | 146s | **NO** (271 chars) | 243s | **YES** (9,177 chars) |
| 3 | 347s | **YES** (17,149 chars) | 543s | **YES** (13,608 chars) |
| 4 | 600s | **NO** (timeout) | 593s | **YES** (14,671 chars) |
| 5 | 240s | **NO** (23,283 chars) | 600s | **NO** (timeout, 0 chars) |
| **Running** | | **1/5** | | **3/5** |

## What Each Condition Actually Did

**WITH-skill** used correct LASER APIs throughout:
- **Step 1**: Full SEIR model with LASER imports — but **over-engineered** (added
  gravity, row_normalizer, ValuesMap not requested). Failed to run (likely
  container path issue with `--print` mode writing to file)
- **Step 2**: Only described parameter tweaks (271 chars) since gravity was
  already in step 1. Not a complete script
- **Step 3**: Full rebuild with seasonal forcing via `ValuesMap.from_timeseries`.
  Ran successfully
- **Step 4**: Timeout (600s)
- **Step 5**: Full 10-year model with `BirthsByCBR`, `MortalityByCDR`, correct
  units. Failed to run

**WITHOUT-skill** abandoned LASER entirely and built a pure numpy SIR model:
- **Step 1**: Timeout (600s) — couldn't figure out LASER
- **Step 2**: Gave up on LASER, built a numpy SIR from scratch with manual
  gravity. Ran
- **Step 3**: Added custom `ValuesMap` class (reimplemented, not imported from
  LASER). Ran
- **Step 4**: Correctly diagnosed and fixed a spatial coupling bug
  (`foi_home = visit_frac @ foi`). Ran
- **Step 5**: Timeout (600s)

## Per-Step Scoring

### Rubric: 5 dimensions
- **AC** (API Correctness, 0-3): Correct LASER v1.0.0 API usage
- **SC** (Structural Correctness, 0-3): Working simulation assembly
- **DA** (Domain Adaptation, 0-3): Correct disease parameters
- **CO** (Completeness, 0-3): Runnable code vs pseudocode
- **CB** (Consistency Bonus, 0-1): Builds on previous step (steps 2-5 only)

### WITH Skill

| Step | AC | SC | DA | CO | CB | Total | Runs? |
|------|----|----|----|----|-----|-------|-------|
| 1 | 3 | 2 | 3 | 3 | -- | 11/12 | No |
| 2 | 2 | 1 | 2 | 0 | 0 | 5/13 | No |
| 3 | 3 | 3 | 3 | 3 | 1 | 13/13 | Yes |
| 4 | 0 | 0 | 0 | 0 | 0 | 0/13 | No (timeout) |
| 5 | 3 | 2 | 3 | 3 | 0 | 11/13 | No |
| **Sum** | **11** | **8** | **11** | **9** | **1** | **40/64** | **1/5** |

### WITHOUT Skill

| Step | AC | SC | DA | CO | CB | Total | Runs? |
|------|----|----|----|----|-----|-------|-------|
| 1 | 0 | 0 | 0 | 0 | -- | 0/12 | No (timeout) |
| 2 | 0 | 2 | 2 | 3 | 0 | 7/13 | Yes |
| 3 | 0 | 2 | 2 | 3 | 1 | 8/13 | Yes |
| 4 | 0 | 3 | 2 | 3 | 1 | 9/13 | Yes |
| 5 | 0 | 0 | 0 | 0 | 0 | 0/13 | No (timeout) |
| **Sum** | **0** | **7** | **6** | **9** | **2** | **24/64** | **3/5** |

## Aggregate by Dimension

| Dimension | WITH (/max) | WITHOUT (/max) | Delta |
|-----------|------------|---------------|-------|
| API Correctness | 11/15 | 0/15 | **+11** |
| Structural Correctness | 8/15 | 7/15 | +1 |
| Domain Adaptation | 11/15 | 6/15 | +5 |
| Completeness | 9/15 | 9/15 | 0 |
| Consistency Bonus | 1/4 | 2/4 | -1 |
| **TOTAL** | **40/64** | **24/64** | **+16** |
| Steps running | 1/5 | 3/5 | **-2** |

## Key Findings

1. **WITH-skill wins on quality (40 vs 24) but loses on execution (1/5 vs 3/5).**
   Correct LASER code that doesn't run is less useful than wrong-framework code
   that does run.

2. **WITHOUT abandoned LASER entirely.** After timing out on step 1, it pivoted
   to a pure numpy SIR model at step 2. Scored 0/15 on API Correctness but
   3/5 on running code. A rational adaptive strategy.

3. **WITH over-engineered step 1**, adding gravity, row_normalizer, and
   ValuesMap when only a basic SEIR was requested. This caused step 2 ("add
   gravity") to produce a 271-char parameter tweak instead of a complete script.
   **Negative transfer from the skill.**

4. **Both conditions lost steps to 600s timeouts** — WITH lost step 4, WITHOUT
   lost steps 1 and 5. The `--print` containerized mode appears slower than
   interactive use.

5. **WITHOUT showed better consistency** (CB: 2/4 vs 1/4). Steps 2→3→4 built
   incrementally on each other with a coherent numpy SIR that accumulated
   features. WITH had two gaps (step 2 incomplete, step 4 timeout) that broke
   the chain.

6. **WITHOUT's step 4 debugging was excellent** — correctly identified the
   spatial coupling bug and fixed it cleanly. Strong generic programming skill
   even without LASER knowledge.

## Comparison to Rubric Hypotheses

| Hypothesis | Actual |
|------------|--------|
| WITH +1.5 avg AC advantage | **+2.2** — even larger than predicted |
| WITHOUT shows progressive degradation | **Wrong** — WITHOUT improved over steps 2-4 |
| WITH maintains code coherence | **Wrong** — step 2 gap and step 4 timeout broke coherence |
| WITHOUT rewrites from scratch | **Partially right** — step 2 was a fresh start, but steps 3-4 were incremental |

## The Execution Paradox

The multi-turn suite exposes a fundamental tension: **API correctness and
executability are inversely correlated in containerized evaluation.** WITH wrote
correct LASER code that failed to run (container environment issues, `--print`
mode limitations, runtime errors). WITHOUT wrote non-LASER code that ran
perfectly because numpy has no environment issues.

---

# Methodology Improvements for Future Runs

## Problem 1: Timeouts (3/10 steps lost)

**Issue:** 600s timeout killed WITH step 4, WITHOUT steps 1 and 5. The Docker +
`--print` mode is significantly slower than interactive use. Each timeout wastes
10 minutes AND cascades (step 5 can't build on a timed-out step 4).

**Fixes:**
- Increase timeout from 600s to 900s or 1200s in `_call_claude_docker()`
- Add retry logic: if a step times out, retry once before recording failure
- Log the timeout explicitly in the response file (currently just "ERROR:
  Claude call timed out after 300 seconds" — the 300s message is from the
  non-Docker path; Docker uses 600s but the error message is wrong)

## Problem 2: n=1 per condition (no statistical power)

**Issue:** With a single run per condition, any result could be random variance.
The WITHOUT condition outperformed WITH on execution (3/5 vs 1/5) — is this
reproducible or a fluke?

**Fixes:**
- Run 3x per condition (like the reruns suite). This means 30 total steps
  instead of 10, but gives variance estimates
- At minimum run 2x to check if the pattern is stable
- Could parallelize: run 3 WITH + 3 WITHOUT simultaneously since they're
  independent conversations

## Problem 3: WITH over-engineers step 1

**Issue:** The skill causes Claude to add gravity, row_normalizer, and
ValuesMap in step 1 when only a basic SEIR is asked for. This means step 2
("add gravity") becomes a parameter tweak instead of a structural addition,
producing incomplete code (271 chars).

**Fixes:**
- Make step 1 prompt more explicit: "Build a BASIC model with NO spatial
  coupling, NO seasonal forcing, NO demographics. Just the 5 SEIR components
  and model.run(). We will add features in subsequent steps."
- Or restructure the step sequence to avoid overlap (e.g., step 1 = SEIR +
  demographics, step 2 = gravity, step 3 = seasonality)

## Problem 4: `--print` doesn't capture file writes

**Issue:** In Docker, Claude sometimes writes code to a file inside the
container rather than outputting it as text. The `--print` flag only captures
text output. This systematically disadvantages the WITH condition since the
skill encourages tool-use workflows (write file → verify → report).

**Fixes:**
- Mount a writable output volume: `-v /tmp/output:/output` and instruct
  Claude to write code there
- After the Docker run, check both stdout AND the output volume for code
- Or: add to the prompt "Output the complete code in a ```python code block
  in your response. Do not write it to a file."
- Alternative: use `--output-format json` if available to capture structured
  output

## Problem 5: Code execution runs on HOST, not in container

**Issue:** `try_run_code()` runs extracted code on the host via
`/opt/anaconda3/bin/python3`. This means:
- WITH code (uses LASER) succeeds if host has LASER installed
- WITHOUT code (uses numpy only) always succeeds
- The execution test doesn't reflect what would happen in the container

**Fix:** Run the extracted code inside the same Docker container where it was
generated. This would make the execution test consistent with the generation
environment. For WITHOUT, this means numpy-only code should still run. For
WITH, this means LASER code runs in the container where LASER is installed.

## Problem 6: Step 4 dynamic prompt creates asymmetric conditions

**Issue:** Step 4's prompt depends on whether step 3 ran. If both conditions'
step 3 ran (as happened here), both get the "identical curves" debugging
prompt. But if one condition's step 3 fails, it gets a different prompt (fix
the error), making the two conditions incomparable at step 4.

**Fixes:**
- Use a fixed step 4 prompt regardless of step 3 outcome. E.g., always use
  the "identical curves" scenario, or always inject a known bug
- Or: keep the dynamic prompt but note in analysis when the two conditions
  received different step 4 prompts

## Problem 7: No resume support for multi-turn

**Issue:** Unlike the other suites, `run-multi-turn.py` has no skip-if-exists
logic. Re-running overwrites all existing good outputs — this destroyed 8 good
outputs in a previous run.

**Fixes:**
- Add file-existence check: if all 5 step response files exist and are >100
  bytes, skip the condition
- Or: add a `--force` flag to overwrite, defaulting to resume
- At minimum, add a `--condition-only` flag to re-run just one condition

## Summary of Recommended Changes

| Priority | Change | Impact |
|----------|--------|--------|
| **High** | Increase timeout to 900-1200s | Eliminates ~30% of failures |
| **High** | Run 3x per condition | Enables statistical comparison |
| **High** | Add resume support to driver | Prevents data loss on re-runs |
| **Medium** | Make step 1 prompt more restrictive | Prevents over-engineering cascade |
| **Medium** | Add "output code in response" instruction | Fixes `--print` capture issue |
| **Medium** | Run extracted code in Docker container | Consistent execution environment |
| **Low** | Fix step 4 asymmetry | Cleaner comparison |
| **Low** | Fix timeout error message (says 300s, actual 600s) | Cosmetic |
