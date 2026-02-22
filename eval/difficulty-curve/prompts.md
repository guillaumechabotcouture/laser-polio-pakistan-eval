# Suite 3: Difficulty Curve — Finding the Complexity Threshold for Framework Abandonment

## Goal

Find the exact complexity threshold where Claude abandons the LASER framework
API and begins hallucinating or reimplementing from scratch. By adding one LASER
concept at a time with monotonically increasing difficulty, we can identify the
"cliff" where API correctness collapses — and whether the LASER skill shifts
that cliff to a higher complexity level.

## Design Rationale

### Monotonically increasing complexity
Each level D1-D8 adds exactly one new LASER concept on top of all previous
requirements. This produces a clean difficulty curve rather than testing
independent skills in isolation (Suite 2: Atomic) or end-to-end integration
(Suite 1: Polio).

### Self-contained prompts
Every prompt is fully self-contained — not "extend the previous model." This
eliminates confounds from context accumulation and ensures each level is an
independent measurement point on the difficulty curve.

### Generic disease avoids knowledge confounds
All prompts use a fictional "Republic of Eastland" with a generic respiratory
disease (R0~5, latent 4 days, infectious 10 days). This avoids confounding
LASER API knowledge with disease-specific epidemiological knowledge, which is
what Suites 1 and 2 already test.

### Single scoring dimension
We score on API Correctness only (0-3), not structural correctness, domain
adaptation, or completeness. This keeps the measurement clean: we are
specifically measuring "does the model know how to call the LASER API?" as a
function of prompt complexity.

## Difficulty Levels

| Level | New Concept | Cumulative API Calls | Key New APIs |
|-------|------------|---------------------|--------------|
| D1 | Model + Scenario + PropertySet | ~2 | `Model`, `PropertySet`, GeoDataFrame |
| D2 | SEIR components + distributions | ~6 | `SEIR.Susceptible/Exposed/Infectious/Recovered`, `SEIR.Transmission`, `dists.gamma` |
| D3 | Transmission with expdurdist | ~8 | `SEIR.Transmission(model, expdurdist=...)` |
| D4 | Births, deaths, capacity | ~11 | `BirthsByCBR`, `MortalityByCDR`, `calc_capacity` |
| D5 | Gravity network + distances | ~14 | `gravity()`, distance computation, `model.network` |
| D6 | Row normalization | ~16 | `row_normalizer(network, max_frac)` |
| D7 | Seasonal forcing via ValuesMap | ~19 | `ValuesMap.from_timeseries()`, `seasonality=` param |
| D8 | Custom component (importation) | ~22 | `__init__(self, model)` + `step(self, tick)` pattern |

## Protocol

For each prompt, run two conditions:

### Condition A: WITH skills (project directory)

```bash
cd laser-polio-pakistan-eval
claude --print --dangerously-skip-permissions < eval/difficulty-curve/prompt-DN.txt \
  > eval/difficulty-curve/outputs/with-skill/prompt-DN.md
```

All skills (`epi-model-parametrization`, `laser-spatial-disease-modeling`,
`modelops-calabaria`) are available via the `.claude/skills/` directory, along
with existing project code and references.

### Condition B: WITHOUT skills (clean temp directory)

```bash
# Create isolated environment with no LASER-specific knowledge
TMPDIR=$(mktemp -d)
cp eval/difficulty-curve/prompt-DN.txt "$TMPDIR/"
cat > "$TMPDIR/CLAUDE.md" << 'EOF'
# Project Context
This is a public health epidemiological modeling project using the open-source
LASER framework (https://laser.idmod.org/) from the Institute for Disease
Modeling. The project builds spatial disease transmission models for public
health planning.
EOF

cd "$TMPDIR"
unset CLAUDECODE  # Required for nested claude sessions
claude --print --dangerously-skip-permissions --disable-slash-commands \
  < prompt-DN.txt > eval/difficulty-curve/outputs/without-skill/prompt-DN.md
```

The `--disable-slash-commands` flag prevents skill loading. The clean temp
directory ensures no access to `.claude/skills/`, existing scripts, or
project-specific CLAUDE.md content beyond basic context.

### General Rules

1. Each prompt runs as a single one-shot generation (no follow-up questions)
2. Pre-extracted prompts are in `eval/difficulty-curve/prompt-D{1-8}.txt`
3. The automated runner `eval/difficulty-curve/run-eval.sh` handles both conditions

## Prompts

### D1 -- Minimal Model (~2 API calls)
**New concept:** Model + Scenario GeoDataFrame + PropertySet
**Tests:** Can Claude instantiate the most basic LASER objects?

```
Using the LASER framework (laser-generic package), create a basic spatial
disease model for a respiratory illness in the Republic of Eastland. Set up a
Model with a 4-patch GeoDataFrame scenario (populations: 100k, 200k, 150k, 80k;
all susceptible) and a PropertySet with nticks=365 and prng_seed=42. Call
model.run() to execute a 1-year simulation. Write complete Python code. Do not
install any packages.
```

### D2 -- Add SEIR Components (~6 API calls)
**New concept:** SEIR state components + Transmission + dists.gamma
**Tests:** Can Claude add the standard five SEIR components in correct order?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k). The disease has R0~5, a latent period of
4 days (gamma distributed, shape=4, scale=1), and an infectious period of 10
days (gamma distributed, shape=5, scale=2). Initialize with 90% susceptible, 0%
exposed, 1% infectious, 9% recovered. Add all SEIR state components
(Susceptible, Exposed, Infectious, Recovered) and Transmission in the correct
order. Run for 1 year. Write complete Python code. Do not install any packages.
```

### D3 -- Transmission with Distributions (~8 API calls)
**New concept:** Wiring expdurdist parameter to Transmission
**Tests:** Can Claude connect distribution objects to Transmission's expdurdist?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k). Disease parameters: R0~5, latent period
gamma(shape=4, scale=1) mean 4 days, infectious period gamma(shape=5, scale=2)
mean 10 days. Initialize: 90% S, 0% E, 1% I, 9% R. Add all SEIR components in
correct order. The Transmission component should use the gamma-distributed
infectious period as its expdurdist parameter. Run for 1 year. Write complete
Python code. Do not install any packages.
```

### D4 -- Add Births and Deaths (~11 API calls)
**New concept:** BirthsByCBR + MortalityByCDR + calc_capacity
**Tests:** Can Claude use demographic components with per-1000/year rates?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k). Disease: R0~5, latent gamma(shape=4,
scale=1), infectious gamma(shape=5, scale=2). Init: 90% S, 0% E, 1% I, 9% R.
Add all SEIR components + Transmission in correct order. Add demographics:
births at crude birth rate 30 per 1000/year using BirthsByCBR, deaths at crude
death rate 10 per 1000/year using MortalityByCDR. Use calc_capacity to
pre-allocate agent storage for the birth rate over 10 years. Run for 10 years
(nticks=3650). Write complete Python code. Do not install any packages.
```

### D5 -- Add Gravity Network (~14 API calls)
**New concept:** Distance computation + gravity() + model.network
**Tests:** Can Claude build and attach a gravity-model migration network?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k) at coordinates (30N, 50E), (30.5N,
50.5E), (31N, 51E), (31.5N, 51.5E). Disease: R0~5, latent gamma(shape=4,
scale=1), infectious gamma(shape=5, scale=2). Init: 90% S, 0% E, 1% I, 9% R.
Add all SEIR components + Transmission. Add BirthsByCBR(CBR=30) and
MortalityByCDR(CDR=10). Use calc_capacity for 10 years. Compute pairwise
distances using LASER's distance function. Build a gravity network with k=0.01,
a=1, b=1, c=1.5 and set it on the model. Run for 10 years. Write complete
Python code. Do not install any packages.
```

### D6 -- Add Row Normalization (~16 API calls)
**New concept:** row_normalizer(network, max_fraction)
**Tests:** Can Claude apply row normalization to constrain network flows?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k) at coordinates (30N, 50E), (30.5N,
50.5E), (31N, 51E), (31.5N, 51.5E). Disease: R0~5, latent gamma(shape=4,
scale=1), infectious gamma(shape=5, scale=2). Init: 90% S, 0% E, 1% I, 9% R.
Add all SEIR components + Transmission. Add BirthsByCBR(CBR=30) and
MortalityByCDR(CDR=10). Use calc_capacity for 10 years. Compute distances,
build gravity network (k=0.01, a=1, b=1, c=1.5), then apply row_normalizer so
no patch exports more than 15% of its force of infection. Set the normalized
network on the model. Run for 10 years. Write complete Python code. Do not
install any packages.
```

### D7 -- Add Seasonal Forcing via ValuesMap (~19 API calls)
**New concept:** ValuesMap.from_timeseries() + seasonality parameter
**Tests:** Can Claude create and wire a ValuesMap for seasonal transmission?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k) at coordinates (30N, 50E), (30.5N,
50.5E), (31N, 51E), (31.5N, 51.5E). Disease: R0~5, latent gamma(shape=4,
scale=1), infectious gamma(shape=5, scale=2). Init: 90% S, 0% E, 1% I, 9% R.
Add all SEIR components + Transmission. Add BirthsByCBR(CBR=30) and
MortalityByCDR(CDR=10). Use calc_capacity for 10 years. Build gravity network
(k=0.01, a=1, b=1, c=1.5) with row_normalizer(0.15). Create a 365-day seasonal
forcing array peaking at 1.3x in winter (days 0-90) and troughing at 0.7x in
summer (days 150-240), with mean approximately 1.0. Use
ValuesMap.from_timeseries to create a seasonal profile for 4 patches and pass it
to Transmission's seasonality parameter. Run for 10 years. Write complete Python
code. Do not install any packages.
```

### D8 -- Add Custom Importation Component (~22 API calls)
**New concept:** Custom component with __init__(model) + step(tick) pattern
**Tests:** Can Claude write a custom LASER component that integrates with the model?

```
Using the LASER framework (laser-generic package), build a spatial SEIR model
for a respiratory illness in the Republic of Eastland with 4 patches
(populations: 100k, 200k, 150k, 80k) at coordinates (30N, 50E), (30.5N,
50.5E), (31N, 51E), (31.5N, 51.5E). Disease: R0~5, latent gamma(shape=4,
scale=1), infectious gamma(shape=5, scale=2). Init: 90% S, 0% E, 1% I, 9% R.
Add all SEIR components + Transmission. Add BirthsByCBR(CBR=30) and
MortalityByCDR(CDR=10). Use calc_capacity for 10 years. Build gravity network
(k=0.01, a=1, b=1, c=1.5) with row_normalizer(0.15). Add seasonal forcing via
ValuesMap (winter peak 1.3x, summer trough 0.7x). Write a custom
PulseImportation component (following LASER's __init__(model) + step(tick)
pattern) that introduces 3 infections into patch 0 every 60 ticks. Run for 10
years. Write complete Python code. Do not install any packages.
```

## Expected Difficulty Curve Shape

We expect an S-shaped degradation in API correctness as complexity increases:

- **D1-D3 (plateau):** Both conditions likely score well — basic SEIR is in
  training data
- **D4-D5 (divergence):** BirthsByCBR/MortalityByCDR and gravity() are
  LASER-specific; without-skill condition begins to hallucinate
- **D6-D8 (collapse):** row_normalizer, ValuesMap, and custom components
  require specific API knowledge; without-skill condition abandons framework

The skill should shift the collapse point rightward (higher complexity before
abandonment).

## Analysis

After scoring, plot API Correctness (0-3) vs. difficulty level (D1-D8) for
both conditions. The key metrics are:

1. **Collapse point:** First level where score drops below 2 (mostly correct)
2. **Skill delta:** How many levels further the with-skill condition sustains
   correctness
3. **Asymptotic floor:** What score the without-skill condition stabilizes at
   after collapse (0 = total abandonment, 1 = partial knowledge)
