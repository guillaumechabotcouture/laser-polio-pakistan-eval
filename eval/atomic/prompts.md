# Suite 1: Atomic Competency Map — LASER Skill Evaluation

## Overview

This suite tests **individual API patterns in isolation**. Each of the 20
micro-prompts targets exactly one LASER framework concept, making it possible
to identify precisely which API elements the skill teaches effectively versus
which ones Claude already knows (or consistently gets wrong).

Unlike the integration prompts in `eval/prompts.md` which test end-to-end model
assembly, these atomic prompts decompose the LASER API surface into testable
units. This enables fine-grained diagnosis: if an integration test fails, the
atomic results reveal whether the failure stems from a wrong import path, wrong
function signature, wrong units, or wrong conceptual understanding.

## Protocol

For each prompt below, run two conditions:

### Condition A: WITH skills (project directory)

```bash
cd laser-polio-pakistan-eval
claude --print --dangerously-skip-permissions < eval/atomic/prompt-A01.txt \
  > eval/atomic/outputs/with-skill/prompt-A01.md
```

All skills in `.claude/skills/` are available, along with existing project code
and references.

### Condition B: WITHOUT skills (clean temp directory)

```bash
TMPDIR=$(mktemp -d)
cp eval/atomic/prompt-A01.txt "$TMPDIR/"
cat > "$TMPDIR/CLAUDE.md" << 'EOF'
# Project Context
This is a public health epidemiological modeling project using the open-source
LASER framework (https://laser.idmod.org/) from the Institute for Disease
Modeling. The project builds spatial disease transmission models for public
health planning.
EOF

cd "$TMPDIR"
unset CLAUDECODE
claude --print --dangerously-skip-permissions --disable-slash-commands \
  < prompt-A01.txt > eval/atomic/outputs/without-skill/prompt-A01.md
```

The `--disable-slash-commands` flag prevents skill loading. The clean temp
directory ensures no access to `.claude/skills/`, existing scripts, or
project-specific CLAUDE.md content beyond basic context.

### Automated Runner

```bash
./eval/atomic/run-eval.sh A01      # Run one prompt, both conditions
./eval/atomic/run-eval.sh all      # Run all 20 prompts
```

### General Rules

1. Each prompt runs as a single one-shot generation (no follow-up questions)
2. Prompts use generic "spatial disease model" framing (not polio-specific)
3. Each prompt tests exactly ONE API pattern
4. Scoring is binary per-concept: 0 (wrong/missing), 1 (partial), 2 (correct)

## Context

Both conditions receive minimal project context establishing the public health
purpose. The WITH condition additionally has access to all skill files and
existing project code. Prompts are deliberately generic (not polio-specific) to
test API knowledge in isolation rather than domain adaptation.

## Prompts

### A01 — Model(scenario, params)
**Tests:** Core model instantiation with GeoDataFrame + PropertySet

```
Using the LASER framework (laser-generic package), create a Model instance with
a 4-patch GeoDataFrame scenario (each patch 50,000 population, all susceptible)
and a PropertySet with nticks=365 and prng_seed=20250101. Show the complete
Python code. Do not install any packages.
```

### A02 — GeoDataFrame scenario construction
**Tests:** Correct scenario columns (population, S, E, I, R, geometry)

```
Using the LASER framework, build a GeoDataFrame scenario for a 6-patch spatial
model. Each patch should have columns for population (varying from 10k to 100k),
initial S, E, I, R counts (98% S, 1% E, 1% I, 0% R), and geometry (Point
coordinates). Show the complete Python code. Do not install any packages.
```

### A03 — PropertySet
**Tests:** PropertySet creation with keyword args, parameter access

```
Using the LASER framework, create a PropertySet containing: nticks=3650,
beta=0.25, gamma=1/14, sigma=1/5, prng_seed=42, and a custom parameter
seasonal_amplitude=1.3. Show how to access individual parameters. Do not install
any packages.
```

### A04 — dists.gamma
**Tests:** Gamma distribution creation with shape/scale parameterization

```
Using the LASER framework's distribution utilities, create a gamma distribution
with mean 10 days and shape parameter 5 (so scale=2). Show how to create it and
explain what it would be used for in an SEIR model (e.g., infectious period
duration). Do not install any packages.
```

### A05 — dists.normal
**Tests:** Normal distribution creation with loc/scale parameterization

```
Using the LASER framework's distribution utilities, create a normal distribution
with mean 28 days and standard deviation 5 days. Show how to create it using
LASER's dists module. Do not install any packages.
```

### A06 — gravity network
**Tests:** gravity() function with correct 6-arg signature (pops, dists, k, a, b, c)

```
Using the LASER framework, build a gravity-model migration network for 4 patches
with populations [50000, 100000, 75000, 200000] arranged in a line 100km apart.
Use parameters k=0.01, a=1, b=1, c=2. Show the complete code to compute the
distance matrix and gravity network. Do not install any packages.
```

### A07 — row_normalizer
**Tests:** Row normalization function with (network, max_frac) signature

```
Using the LASER framework, given a gravity network matrix for 4 patches, apply
row normalization so that no patch exports more than 20% of its force of
infection. Show the complete code using LASER's row_normalizer function. Do not
install any packages.
```

### A08 — ValuesMap.from_timeseries
**Tests:** ValuesMap creation from time-series array for multiple patches

```
Using the LASER framework, create a seasonal transmission profile using
ValuesMap. Build a 365-day array where transmission peaks at 1.3x in summer
(days 180-270) and troughs at 0.7x in winter (days 0-90, 330-365), then create
a ValuesMap for 4 patches. Do not install any packages.
```

### A09 — SEIR.Transmission setup
**Tests:** Transmission component with expdurdist + seasonality kwargs

```
Using the LASER framework, set up the SEIR.Transmission component for a model.
Use a gamma-distributed infectious period with shape=5, scale=5.6 (mean 28
days). Include a seasonal ValuesMap. Show the complete code for creating the
distribution and initializing the Transmission component. Do not install any
packages.
```

### A10 — SEIR component ordering
**Tests:** Correct addition order (S, E, I, R, Transmission) and rationale

```
Using the LASER framework, add all four SEIR state-tracking components
(Susceptible, Exposed, Infectious, Recovered) plus Transmission to a Model
instance. Show the correct execution order for these components and explain why
the order matters. Do not install any packages.
```

### A11 — BirthsByCBR (units)
**Tests:** CBR unit convention (per-1000/year, NOT pre-converted daily per-capita)

```
Using the LASER framework, add a birth process to a model using BirthsByCBR.
The crude birth rate is 44 per 1000 population per year. Show the complete code.
Be explicit about what units BirthsByCBR expects for the birth rate parameter.
Do not install any packages.
```

### A12 — MortalityByCDR (units)
**Tests:** CDR unit convention (per-1000/year, NOT pre-converted daily per-capita)

```
Using the LASER framework, add a mortality process to a model using
MortalityByCDR. The crude death rate is 13 per 1000 population per year. Show
the complete code. Be explicit about what units MortalityByCDR expects. Do not
install any packages.
```

### A13 — calc_capacity (units)
**Tests:** Pre-allocation calculation with per-1000/year birthrate

```
Using the LASER framework, use calc_capacity to pre-allocate agent storage for a
model. Initial population is 500,000 with a crude birth rate of 44 per 1000 per
year over a 10-year simulation. Show the complete code and explain what units
calc_capacity expects. Do not install any packages.
```

### A14 — Custom component pattern
**Tests:** __init__(self, model) + step(self, tick) protocol

```
Using the LASER framework, write a custom component class called PulseImportation
that introduces 5 new infections into patch 0 every 90 ticks. The component must
follow LASER's component protocol with __init__(self, model) and step(self, tick)
methods. Do not install any packages.
```

### A15 — AliasedDistribution
**Tests:** Discrete distribution sampling for age pyramids

```
Using the LASER framework, create an AliasedDistribution representing an age
pyramid: 30% aged 0-4, 25% aged 5-14, 20% aged 15-29, 15% aged 30-49, 10%
aged 50+. Show how to create it and sample 1000 ages from it. Do not install any
packages.
```

### A16 — distance (haversine)
**Tests:** LASER's haversine distance function import and usage

```
Using the LASER framework, compute the haversine distance in kilometers between
Lahore (31.55N, 74.35E) and Karachi (24.86N, 67.01E) using LASER's distance
function. Show the complete code. Do not install any packages.
```

### A17 — State enum values
**Tests:** SUSCEPTIBLE=0, EXPOSED=1, INFECTIOUS=2, RECOVERED=3

```
Using the LASER framework, write code that explicitly shows the integer values of
all four SEIR states (SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED). Show where to
import them from and what numeric value each state has. Do not install any
packages.
```

### A18 — model.network manual setup
**Tests:** Assigning gravity network to model after construction

```
Using the LASER framework, after creating a Model instance, manually set its
gravity-based migration network. Show how to compute a gravity network for 4
patches and assign it to the model's network attribute. Do not install any
packages.
```

### A19 — Vaccination gotcha
**Tests:** Must set state=RECOVERED (not susceptibility) for TransmissionSE

```
Using the LASER framework, write code to vaccinate 80% of susceptible agents in a
model so they become immune to transmission. IMPORTANT: LASER's TransmissionSE
component checks the agent state field, not susceptibility. Show the correct way
to make vaccinated agents immune. Do not install any packages.
```

### A20 — Seasonal normalization
**Tests:** Building a mean-1.0 seasonal array for ValuesMap

```
Using the LASER framework, create a 365-day seasonal forcing array for disease
transmission where the peak is 1.5x the baseline and the trough is 0.5x the
baseline, with the annual mean equal to approximately 1.0. Show the complete code
and verify the mean. Do not install any packages.
```

## Competency Categories

The 20 prompts group into six conceptual categories:

| Category | Prompts | What it tests |
|----------|---------|---------------|
| **Core objects** | A01, A02, A03 | Model, scenario GeoDataFrame, PropertySet |
| **Distributions** | A04, A05, A15 | dists.gamma, dists.normal, AliasedDistribution |
| **Spatial coupling** | A06, A07, A16, A18 | gravity(), row_normalizer(), distance(), model.network |
| **SEIR components** | A09, A10, A14 | Transmission setup, component ordering, custom components |
| **Demographics** | A08, A11, A12, A13, A20 | ValuesMap, BirthsByCBR, MortalityByCDR, calc_capacity, seasonality |
| **State semantics** | A17, A19 | State enum values, vaccination gotcha (state vs susceptibility) |
