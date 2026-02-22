# Suites 1 & 3 Analysis Report: Atomic Competency Map + Difficulty Curve

## Executive Summary

**Both suites are INVALID as measures of skill-provided LASER knowledge.**

The WITHOUT-skill condition was designed to test "base LLM knowledge of LASER
API" but actually tests "can Claude Code reverse-engineer an installed Python
package by reading its source code." The `--dangerously-skip-permissions` flag
grants full filesystem and code execution access, allowing Claude to read the
globally installed `laser-generic` source at
`/opt/anaconda3/lib/python3.12/site-packages/laser/` and execute code against
it — even from a clean tmpdir with `--disable-slash-commands`.

The skill's value cannot be measured until this confound is fixed.

---

## The Confound: Source Code Contamination

### Mechanism

```
WITHOUT condition setup:
  1. Create clean tmpdir (no .claude/skills/)
  2. Write minimal CLAUDE.md
  3. Run: claude --print --dangerously-skip-permissions --disable-slash-commands

What --disable-slash-commands blocks:
  ✗ Skill loading from .claude/skills/
  ✗ Slash command invocation

What --dangerously-skip-permissions ALLOWS:
  ✓ Read any file on the filesystem
  ✓ Execute any bash command
  ✓ Write files to disk
  ✓ Run Python against installed packages
```

Claude Code in the WITHOUT condition can:
1. Read `/opt/anaconda3/lib/python3.12/site-packages/laser/**/*.py`
2. Write `.py` files to the tmpdir
3. Execute them with `/opt/anaconda3/bin/python3`
4. Iterate on errors by reading tracebacks

The tmpdir is destroyed after the session (`rm -rf "$tmpdir"`), erasing
all evidence of what Claude explored and wrote.

### Direct Evidence of Source Reading

WITHOUT-skill outputs contain **explicit references to LASER source line
numbers and internal module paths** that could only come from reading the
installed package:

| Prompt | Evidence |
|--------|---------|
| A07 | "The key line in the implementation (`migration.py:137`)" |
| A12 | "Internally (see `vitaldynamics.py:186-190`)" |
| A17 | "The underlying enum (from `laser/generic/shared.py`)" |
| A04 | Preamble: "Now I have the details." |
| A07 | Preamble: "Now I have the full source code." |
| A12 | Preamble: "Now I have the full picture." |

WITHOUT outputs consistently use **internal import paths** rather than the
convenience re-exports, proving source tree exploration:

```python
# WITHOUT-skill uses internal paths (from reading source):
from laser.core.migration import gravity, row_normalizer, distance
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.propertyset import PropertySet
from laser.core.utils import calc_capacity

# WITH-skill uses public API paths (from skill documentation):
from laser.generic import gravity, row_normalizer
from laser.generic import BirthsByCBR, MortalityByCDR
from laser.generic import Model, PropertySet
```

WITHOUT outputs even use **parameter names from the actual source** that
differ from the rubric:
- `max_rowsum` (actual) vs `max_frac` (rubric) for `row_normalizer`
- `nnodes` (actual) vs `num_patches` (rubric) for `ValuesMap.from_timeseries`

### Evidence of Code Execution

All difficulty curve WITHOUT outputs are prose summaries with execution timing
and stochastic numerical results:

| Prompt | Timing | Stochastic Result |
|--------|--------|-------------------|
| D1 | — | "468,697 recovered by day 365" |
| D2 | — | "peaked at 177,751 infectious on day 32" |
| D4 | "~5.6 seconds" | "173,620 births, 59,474 deaths" |
| D7 | "~5 seconds" | "475,836 total infections" |
| D8 | "~6 seconds" | "Peak infectious: 211,588 at tick 28" |

These stochastic numbers differ across prompts (different seeds/parameters)
and match realistic epidemiological dynamics. They could only come from
actual execution of working code against the installed `laser-generic`.

---

## What the Eval Actually Measures

| Intended Measurement | Actual Measurement |
|---------------------|-------------------|
| WITH: Skill + LLM knowledge | WITH: Skill + LLM + source reading + execution |
| WITHOUT: LLM knowledge only | WITHOUT: LLM + source reading + execution |
| Delta: Value of the skill | Delta: Value of skill minus value of source reading |

The skill's apparent value is **deflated** because the WITHOUT condition has
an alternative knowledge source (the installed package) that partially
substitutes for the skill. The true skill value is higher than measured.

---

## Raw Scores (Unreliable Due to Confound)

### Suite 1: Atomic Competency Map

**WITH: 40/40 (100%) | WITHOUT: 28-36/40 (70-90%)**

The WITHOUT score depends on scoring strictness:
- **Lenient (initial scoring):** 28/40 — accepted summaries at face value
- **Strict (re-audit):** 36/40 — acknowledged that code actually ran correctly

| Prompt | WITH | WITHOUT | Delta | WITHOUT Knowledge Source |
|--------|:----:|:-------:|:-----:|------------------------|
| A01 | 2 | 1-2 | 0-1 | Source reading (internal import paths) |
| A02 | 2 | 2 | 0 | Standard geopandas knowledge |
| A03 | 2 | 2 | 0 | Source reading (actual constructor API) |
| A04 | 2 | 2 | 0 | Source reading ("Now I have the details") |
| A05 | 2 | 2 | 0 | Source reading (internal path) |
| A06 | 2 | 1-2 | 0-1 | Source reading + manual implementation |
| A07 | 2 | 2 | 0 | Source reading (cites migration.py:137) |
| A08 | 2 | 2 | 0 | Source reading (uses actual param `nnodes`) |
| A09 | 2 | 1-2 | 0-1 | Source reading (internal paths, extra params) |
| A10 | 2 | 0-1 | 1-2 | Source reading (but got ordering wrong) |
| A11 | 2 | 2 | 0 | Source reading (cites vitaldynamics.py) |
| A12 | 2 | 2 | 0 | Source reading (cites lines 186-190) |
| A13 | 2 | 2 | 0 | Source reading (internal path) |
| A14 | 2 | 1 | 1 | Ambiguous (wrong state value) |
| A15 | 2 | 1-2 | 0-1 | Source reading (demographics module) |
| A16 | 2 | 2 | 0 | Source reading (internal path) |
| A17 | 2 | 2 | 0 | Source reading (cites shared.py) |
| A18 | 2 | 0-2 | 0-2 | Source reading (explored measles subpackage) |
| A19 | 2 | 2 | 0 | Source reading (understood TransmissionSE) |
| A20 | 2 | 2 | 0 | Standard numpy knowledge |

### Suite 3: Difficulty Curve

**WITH: 24/24 (100%) | WITHOUT: 18-22/24 (75-92%)**

| Level | WITH | WITHOUT | Delta | Notes |
|-------|:----:|:-------:|:-----:|-------|
| D1 | 3 | 2-3 | 0-1 | SIR framing but code ran |
| D2 | 3 | 2-3 | 0-1 | Component order may differ |
| D3 | 3 | 2-3 | 0-1 | expdurdist wiring unclear |
| D4 | 3 | 3 | 0 | Demographics correct, code ran |
| D5 | 3 | 3 | 0 | Gravity correct, code ran |
| D6 | 3 | 2-3 | 0-1 | row_normalizer correct |
| D7 | 3 | 2-3 | 0-1 | Seasonal forcing works |
| D8 | 3 | 2-3 | 0-1 | Custom component, code ran |

### Category Analysis (Atomic)

| Category | WITH | WITHOUT | Apparent Delta | Real Delta |
|----------|:----:|:-------:|:-:|:-:|
| Core objects | 6/6 | 4-6/6 | 0-2 | Unknown |
| Distributions | 6/6 | 5-6/6 | 0-1 | Unknown |
| Spatial coupling | 8/8 | 5-8/8 | 0-3 | Unknown |
| SEIR components | 6/6 | 2-4/6 | 2-4 | Unknown |
| Demographics | 10/10 | 10/10 | **0** | Unknown |
| State semantics | 4/4 | 4/4 | **0** | Unknown |

---

## What IS Reliable From These Results

Despite the confound, some findings are trustworthy:

### 1. WITH-skill produces perfect scores (40/40, 24/24)

The WITH condition is valid: the skill provides correct API knowledge, Claude
follows it, and the resulting code is correct. This confirms the skill works.

### 2. Component ordering is a real knowledge gap

Even with source code access, the WITHOUT condition got component ordering
wrong (A10: score 0-1). This is the ONE thing that source reading didn't
easily reveal — because the correct ordering is a semantic convention, not
something obvious from reading individual component implementations.

### 3. The skill provides faster, more reliable results

The WITH condition uses the skill's documented public API paths. The WITHOUT
condition had to explore the source tree, read internal modules, and
reverse-engineer the API. The skill substitutes instant documentation for
expensive runtime exploration.

### 4. Rubric has inaccuracies

The WITHOUT condition (reading actual source) sometimes produced MORE correct
answers than the rubric expected:
- `max_rowsum` (actual) vs `max_frac` (rubric)
- `nnodes` (actual) vs `num_patches` (rubric)
- `PropertySet(*dicts)` (actual) vs `PropertySet(keyword=args)` (rubric)

These rubric errors should be corrected regardless of the eval redesign.

---

## Recommended Fix

The WITHOUT condition must be isolated from the installed `laser-generic`
package. Options:

### Option A: Virtual Environment Without LASER (Recommended)

```bash
# Create a clean venv without laser-generic
python -m venv /tmp/laser-eval-clean-env
source /tmp/laser-eval-clean-env/bin/activate
# Do NOT install laser-generic
# Run WITHOUT condition using this Python
```

### Option B: Docker Container

```bash
docker run --rm -v $(pwd):/workspace python:3.12-slim \
  bash -c "pip install claude-code && cd /workspace && \
  claude --print --dangerously-skip-permissions < prompt.txt"
```

### Option C: Block Filesystem Access

Add `--disallowedTools` to prevent reading the package:
```bash
claude --print --dangerously-skip-permissions \
  --disable-slash-commands \
  --disallowedTools "Read(/opt/*)" \
  < prompt.txt
```

### Option D: Uninstall Before Running

```bash
pip uninstall laser-generic -y
# Run WITHOUT condition
pip install laser-generic  # Reinstall for WITH condition
```

**Option A is recommended** because it's clean, reversible, and doesn't
affect the system Python. Options B-D have various drawbacks (Docker
complexity, incomplete blocking, risk of forgetting to reinstall).

However, note that even Option A doesn't prevent Claude from reading
documentation via web search (if available). The cleanest test would also
use `--disallowedTools "WebSearch WebFetch"` in the WITHOUT condition.

---

## Impact on Existing Polio/Guinea Worm Evals

**The same confound likely affects the original A/B tests** (polio 60/60 vs
25/60, guinea worm 57/60 vs 29/60). However, those tests showed MUCH larger
deltas (58% vs 42%, 95% vs 48%), suggesting the confound was less impactful
there — possibly because the bundled prompts required more integrated
knowledge that source reading alone couldn't easily provide.

The atomic tests' small delta (0-30%) makes them more vulnerable to this
confound than the bundled tests' larger delta (53-58%).

---

## Conclusion

**These eval results cannot be used to measure the skill's value.** The
WITHOUT condition's high scores reflect Claude Code's ability to
reverse-engineer installed packages, not base LLM knowledge of LASER.

The eval infrastructure (prompts, rubrics, runners) is sound. Only the
isolation mechanism needs fixing. Once the WITHOUT condition is properly
sandboxed from the installed package, re-running both suites should produce
the expected larger deltas that measure actual skill contribution.
