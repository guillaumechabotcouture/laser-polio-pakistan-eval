# Multi-Turn Iterative Evaluation: 5-Step Conversation Protocol

## Purpose

This suite tests whether the LASER skill helps with **iterative development**
-- building a model step by step with feedback between turns. Unlike the polio
and guinea worm suites (one-shot prompts), this sends 5 sequential prompts
where each builds on the previous output. The driver script extracts code from
each response, attempts to run it, and feeds the result back as context for the
next step.

This evaluates a fundamentally different capability: maintaining coherent,
accumulating code across multiple turns while progressively adding complexity.

## Key Differences from One-Shot Suites

| Aspect | One-Shot (Suites 1-3) | Multi-Turn (Suite 4) |
|--------|----------------------|---------------------|
| Prompts | Independent, no context | Sequential, builds on previous |
| Disease | Polio/guinea worm (specific) | Generic respiratory (lighter domain) |
| Evaluation focus | API + domain knowledge | Consistency + progressive complexity |
| Code continuity | Each prompt standalone | Each step must preserve previous work |
| Error handling | Not tested | Step 4 explicitly tests recovery |
| Scoring | 4 dimensions | 4 dimensions + consistency bonus |

## Protocol

For each condition (WITH-skill, WITHOUT-skill), run all 5 steps sequentially.
Each step's prompt includes the code from the previous step as context. The
automated driver (`run-multi-turn.py`) handles this orchestration.

### Condition A: WITH skills (project directory)

The driver runs Claude from the project root directory where all skills
(`.claude/skills/`) are available.

### Condition B: WITHOUT skills (clean temp directory)

The driver runs Claude from a temporary directory with only a minimal
`CLAUDE.md` and `--disable-slash-commands` to prevent skill loading.

---

## Step 1 -- Base Model

**Tests:** LASER API fundamentals, correct imports, component assembly, basic SEIR

### Prompt

```
Using the LASER framework (laser-generic package), build a basic 4-patch SEIR
model for a respiratory disease. Patches have populations 100k, 200k, 150k,
80k. Disease: R0~5, latent period 4 days, infectious period 10 days. Initialize
with 90% susceptible, 1% infectious, 9% recovered. Run for 1 year. Write
complete Python code. Do not install any packages.
```

### What to look for

- Correct imports: `from laser.generic import Model, SEIR, PropertySet, dists`
- Scenario as GeoDataFrame with population, S, E, I, R, geometry columns
- Components in correct order: Susceptible, Exposed, Infectious, Recovered, Transmission
- `PropertySet(nticks=365, beta=..., prng_seed=...)` with appropriate beta for R0~5
- `model.run()` called
- Heterogeneous populations (not all identical)

---

## Step 2 -- Add Gravity Coupling

**Tests:** Spatial coupling, gravity model, row normalization, code integration

### Prompt

```
Add gravity-model spatial coupling to this model. The 4 patches are arranged
in a line 75km apart. Use gravity parameters k=0.01, a=1, b=1, c=1.5.
Row-normalize so no patch exports more than 15%. Update the code and show the
complete modified script.
```

### What to look for

- `from laser.generic import gravity, row_normalizer`
- `gravity(populations, distances, k, a, b, c)` with 6 positional args
- Distance matrix correctly computed for linear arrangement (75, 150, 225 km)
- `row_normalizer(network, 0.15)` applied
- Network passed to model via `params` or `scenario`
- All Step 1 code preserved and functional

---

## Step 3 -- Add Seasonal Forcing

**Tests:** ValuesMap usage, seasonal transmission, code accumulation

### Prompt

```
Add seasonal forcing to the transmission. Winter peak (days 0-90) at 1.3x
baseline, summer trough (days 150-240) at 0.7x baseline. Use LASER's ValuesMap
to create the seasonal profile. Update the code and show the complete modified
script.
```

### What to look for

- `from laser.generic import ValuesMap`
- `ValuesMap.from_timeseries(array, num_patches)` with correct signature
- Seasonal array of length 365 with smooth or piecewise transitions
- Winter peak (days 0-90) at 1.3x, summer trough (days 150-240) at 0.7x
- Seasonality integrated into Transmission component
- Steps 1-2 code preserved and functional

---

## Step 4 -- Error Recovery (DYNAMIC)

**Tests:** Debugging, error diagnosis, fix quality, preservation of working code

This step's prompt is **generated dynamically** by the driver script based on
whether Step 3's code actually ran successfully.

### If Step 3 code ERRORS:

```
The model crashes with the following error:

[actual error message from Step 3 execution]

Fix the issue and show the complete corrected script. Do not install any
packages.
```

### If Step 3 code RUNS:

```
The model runs but all patches show identical infection curves, suggesting
spatial coupling isn't working. Debug and fix the issue. Show the complete
corrected script. Do not install any packages.
```

### What to look for

- Correctly identifies the root cause (not a superficial fix)
- Preserves all working functionality from Steps 1-3
- Does not introduce new bugs while fixing the reported issue
- If error recovery: fix addresses the actual traceback, not a different issue
- If debugging: investigates coupling parameters, network initialization, or
  transmission integration

---

## Step 5 -- Add Demographics + Verify

**Tests:** BirthsByCBR/MortalityByCDR, calc_capacity, long simulation, self-verification

### Prompt

```
Add births (CBR=30 per 1000/year) and deaths (CDR=10 per 1000/year) using
LASER's BirthsByCBR and MortalityByCDR. Use calc_capacity to pre-allocate for
10 years. Extend the simulation to 10 years. After running, print the total
population at the start and end -- does the ~2% annual growth rate look
correct? Show the complete script.
```

### What to look for

- `from laser.generic import BirthsByCBR, MortalityByCDR`
- `from laser.generic import calc_capacity` (or equivalent)
- CBR=30 and CDR=10 passed as per-1000/year (NOT pre-converted to daily rates)
- `calc_capacity` called with per-1000/year birthrate
- `nticks` updated to 3650 (10 years)
- Population verification: prints start and end population, checks ~2% growth
- All Steps 1-4 code preserved and functional
- No regression in spatial coupling or seasonal forcing

---

## Execution Notes

The driver script (`run-multi-turn.py`) handles:

1. Sending Step 1 prompt to Claude
2. Extracting Python code from the response
3. Attempting to execute the extracted code
4. Including the extracted code as context in the next step's prompt
5. Dynamically generating Step 4 based on Step 3 execution results
6. Saving all prompts, responses, extracted code, and execution logs

Each step uses `/opt/anaconda3/bin/python3` for code execution with a 60-second
timeout.
