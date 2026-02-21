# A/B Test Prompts: LASER Skill Evaluation

## Protocol

For each prompt below, run two conditions:

### Condition A: WITH skills (project directory)

```bash
cd laser-polio-pakistan-eval
claude --print --dangerously-skip-permissions < eval/prompt-N.txt \
  > eval/outputs/with-skill/prompt-N.md
```

All three skills (`epi-model-parametrization`, `laser-spatial-disease-modeling`,
`modelops-calabaria`) are available via the `.claude/skills/` directory, along
with existing project code and references.

### Condition B: WITHOUT skills (clean temp directory)

```bash
# Create isolated environment with no LASER-specific knowledge
TMPDIR=$(mktemp -d)
cp eval/prompt-N.txt "$TMPDIR/"
cat > "$TMPDIR/CLAUDE.md" << 'EOF'
# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation. The work supports polio eradication efforts by building
spatial transmission models using the LASER framework (laser-generic package).
EOF

cd "$TMPDIR"
unset CLAUDECODE  # Required for nested claude sessions
claude --print --dangerously-skip-permissions --disable-slash-commands \
  < prompt-N.txt > eval/outputs/without-skill/prompt-N.md
```

The `--disable-slash-commands` flag prevents skill loading. The clean temp
directory ensures no access to `.claude/skills/`, existing scripts, or
project-specific CLAUDE.md content beyond basic context.

### General Rules

1. Each prompt runs as a single one-shot generation (no follow-up questions)
2. Pre-extracted prompts are in `eval/prompt-{1-5}.txt`
3. The automated runner `eval/run-eval.sh` handles both conditions

## Context

Both conditions receive minimal project context establishing the Gates
Foundation public health purpose. The WITH condition additionally has access to
all skill files and existing project code.

## Prompts

### Prompt 1: Basic Model Setup
**Tests:** LASER API knowledge, correct imports, component signatures

```
For our Gates Foundation polio eradication modeling project, use the LASER
framework (laser-generic package) to set up a basic spatial SEIR model for
polio transmission dynamics across 10 patches representing districts. Each
patch has a population of 100,000. Use a transmission rate corresponding to
R0≈6, a latent period of 3 days, and an infectious period of 28 days.
Initialize with 95% recovered (immune) and 5 infectious per patch. Write
the complete Python code to configure and run a 10-year simulation.
Do not install any packages.
```

### Prompt 2: Gravity Network + Seasonal Forcing
**Tests:** Migration model setup, seasonal transmission, spatial coupling

```
Extend the basic LASER polio model to include:
1. A gravity-model migration network between patches, with distance decay
   exponent c=1.5 and destination population exponent b=0.5
2. Seasonal forcing that peaks during July-October (monsoon season) with
   amplitude 1.3x baseline and troughs during December-March at 0.7x baseline
3. Row-normalize the network so no patch exports more than 15% of its
   force of infection

Write complete Python code using the LASER framework. Assume patches are
arranged in a line 50km apart. Do not install any packages.
```

### Prompt 3: Vaccination Components
**Tests:** Immunization API, domain adaptation (OPV/IPV), custom component design

```
Add vaccination to the LASER polio model:
1. Routine immunization at 6 weeks of age with 80% coverage (OPV)
2. Supplementary Immunization Activities (SIAs) every 6 months targeting
   children aged 0-5 years with 90% coverage
3. The model should track vaccinated individuals separately from naturally
   recovered, since OPV provides mucosal immunity while natural infection
   provides stronger, longer-lasting immunity

Show the complete component setup using LASER's immunization classes where
possible, and custom components where needed. Do not install any packages.
```

### Prompt 4: Calibration Framework
**Tests:** Calibration metric design, adaptation from measles CCS to polio AFP

```
Design a calibration framework for the LASER polio model that:
1. Compares simulated incidence to AFP (Acute Flaccid Paralysis) surveillance
   data, accounting for the ~1:200 paralysis-to-infection ratio
2. Uses a fitness metric based on the proportion of zero-incidence weeks per
   district (similar to critical community size analysis)
3. Samples from these parameter ranges:
   - beta: U(0.15, 0.30)
   - gravity_k: 10^U(-3, -1)
   - seasonal_amplitude: U(1.0, 1.5)
4. Ranks simulations by goodness of fit

Write the calibration loop and scoring functions. Use LASER framework
conventions. Do not install any packages.
```

### Prompt 5: Full Integration (Hardest)
**Tests:** End-to-end model assembly, all components working together

```
Build a complete spatial polio transmission model for Pakistan using the LASER
framework with these specifications:
- 20 districts with heterogeneous populations (range 50k-2M)
- SEIR dynamics with R0≈6, 3-day latent period, 28-day infectious period
- Gravity-model spatial coupling
- Monsoon-season transmission forcing (peak Jul-Oct)
- OPV routine immunization at 80% coverage
- SIA campaigns every 6 months
- Importation of 3 infections every 60 days for the first 5 years
- 20-year simulation with 10-year burn-in
- Output: weekly incidence per district, proportion of zero-incidence weeks

Write complete, runnable Python code. Do not install any packages.
```
