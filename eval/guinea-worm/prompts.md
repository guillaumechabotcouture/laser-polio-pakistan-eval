# A/B Test Prompts: Guinea Worm Stress Test of Restructured LASER Skill

## Purpose

This evaluation stress-tests the restructured LASER skill (3-layer:
discipline/scaffolding/verification) using a fundamentally different disease.
The skill never mentions guinea worm — success means the skill teaches **how to
use LASER as a framework** rather than providing disease-specific recipes.

## Key Differences from Polio Test

| Aspect | Polio (original test) | Guinea Worm (stress test) |
|--------|----------------------|--------------------------|
| Compartmental model | SEIR (lasting immunity) | SEIS (no immunity, reinfection) |
| Transmission route | Fecal-oral / respiratory | Water-borne (copepods in stagnant water) |
| Incubation | ~3 days | ~365 days (pre-patent period) |
| Infectious period | ~28 days | ~21 days (worm emergence) |
| R0 | ~6 | ~1.5-2.5 |
| Hosts | Human only | Human + dogs (dual-host, dogs dominant) |
| Interventions | Vaccination (OPV/IPV) | Prevention only (ABATE, filters, tethering) |
| Seasonality | Monsoon peak (Jul-Oct) | Dry season peak (Jan-Mar) |
| Case numbers | ~36/year paralytic (Pakistan) | ~10-30 human, ~500-1500 dog/year (Chad) |
| Context | Endemic transmission | Near-elimination |
| Immunity | Lifelong after infection | None — reinfection occurs |

## Protocol

For each prompt below, run two conditions:

### Condition A: WITH skills (project directory)

```bash
cd laser-polio-pakistan-eval
claude --print --dangerously-skip-permissions < eval/guinea-worm/prompt-N.txt \
  > eval/guinea-worm/outputs/with-skill/prompt-N.md
```

All skills (`epi-model-parametrization`, `laser-spatial-disease-modeling`,
`modelops-calabaria`) are available via the `.claude/skills/` directory, along
with existing project code and references.

### Condition B: WITHOUT skills (clean temp directory)

```bash
# Create isolated environment with no LASER-specific knowledge
TMPDIR=$(mktemp -d)
cp eval/guinea-worm/prompt-N.txt "$TMPDIR/"
cat > "$TMPDIR/CLAUDE.md" << 'EOF'
# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation's Global Health division. The work supports guinea worm
eradication efforts by building spatial transmission models to evaluate
intervention strategies across Chad's endemic districts. All modeling uses the
open-source LASER framework (https://laser.idmod.org/) from the Institute for
Disease Modeling. This is purely computational epidemiology for public health
planning.
EOF

cd "$TMPDIR"
unset CLAUDECODE  # Required for nested claude sessions
claude --print --dangerously-skip-permissions --disable-slash-commands \
  < prompt-N.txt > eval/guinea-worm/outputs/without-skill/prompt-N.md
```

The `--disable-slash-commands` flag prevents skill loading. The clean temp
directory ensures no access to `.claude/skills/`, existing scripts, or
project-specific CLAUDE.md content beyond basic context.

### General Rules

1. Each prompt runs as a single one-shot generation (no follow-up questions)
2. Pre-extracted prompts are in `eval/guinea-worm/prompt-{1-5}.txt`
3. The automated runner `eval/guinea-worm/run-eval.sh` handles both conditions

## Prompts

### Prompt 1: Basic Dual-Host SEIS Setup
**Tests:** LASER API knowledge, creative adaptation for multi-host, SEIS (not SEIR), correct disease parameters

```
For our Gates Foundation guinea worm eradication modeling project, use the LASER
framework (laser-generic package) to set up a spatial disease model for guinea
worm (Dracunculus medinensis) transmission in Chad. Key disease biology:

- Guinea worm confers NO lasting immunity — recovered individuals return to
  susceptible (use SEIS dynamics, not SEIR)
- Pre-patent period: 10-14 months (from copepod ingestion to worm emergence)
- Infectious period: 2-4 weeks (during worm emergence through skin)
- R0 ≈ 1.5-2 (low, water-mediated transmission through copepods)
- Transmission is water-mediated: infected hosts enter water sources, release
  larvae, copepods ingest larvae, new hosts drink contaminated water

Model 8 patches representing endemic districts in southern Chad. Track two
host populations: humans (total ~2M across patches, range 100k-500k per
district) and dogs (~200K total, range 10k-50k per district), since dogs are
the dominant reservoir in Chad (~95% of detected infections are in dogs).

Initialize in a near-elimination state reflecting current conditions: seed
2 infected humans and 50 infected dogs nationally, distributed across the
3 most endemic patches. All others susceptible.

Write complete Python code to run a 5-year simulation. Output annual case
counts for both species. Do not install any packages.
```

### Prompt 2: Water-Source Transmission + Seasonality
**Tests:** Adapting spatial coupling for water-borne transmission, cross-host dynamics, seasonal forcing

```
Extend the guinea worm LASER model to include realistic transmission dynamics:

1. Water-source mediated transmission: Guinea worm spreads when infected hosts
   enter stagnant water to soothe emerging worms, releasing larvae that are
   ingested by copepods, which are then consumed by new hosts drinking
   unfiltered water. Model spatial coupling based on shared water sources —
   neighboring districts share seasonal ponds and rivers, with coupling
   strength decaying with distance (exponent c=2.0, destination population
   exponent b=0.3)

2. Dry-season transmission forcing: Stagnant pools during Chad's dry season
   (Oct-May) concentrate copepods and increase transmission. Peak transmission
   Jan-Mar at 1.8x baseline, declining through Apr-May. Rainy season (Jun-Sep)
   at 0.4x baseline as flowing water dilutes copepod density and disperses
   larvae.

3. Cross-host transmission via shared water sources: Dogs and humans
   contaminate the same water bodies. In Chad, dog infections drive the
   epidemic (~95% of cases), so dog-to-water-to-human transmission is
   critical. Use asymmetric cross-species coupling reflecting dogs' higher
   contribution.

Maintain SEIS dynamics (no lasting immunity — recovered hosts return to
susceptible). Write complete Python code using the LASER framework.
Do not install any packages.
```

### Prompt 3: Prevention-Based Interventions
**Tests:** Implementing prevention-based (non-vaccine) interventions using LASER components with SEIS dynamics

```
Add prevention-based interventions to the guinea worm LASER model. No vaccine
exists for guinea worm — all control relies on breaking the water transmission
cycle:

1. ABATE (temephos) larvicide: Applied monthly to known water sources during
   transmission season, kills copepods, reduces water-source transmission by
   80% in treated areas. Coverage varies by district (20-80%).

2. Cloth/pipe water filters: Distributed to households, removes copepods from
   drinking water. Current adoption ~60% nationally, each filter reduces
   individual ingestion risk by 95%.

3. Case containment: When a worm emerges, the case is detected and the patient
   is prevented from entering water sources for the ~30-day emergence period.
   70% of human cases are successfully contained, reducing their contribution
   to water contamination by 90%.

4. Dog tethering/management: 50% of villages with known dog infections tether
   infected dogs to prevent water access during worm emergence, reducing
   dog-to-water transmission by 85%.

Show the complete component setup using LASER's framework. These must work
with SEIS dynamics (recovered individuals return to susceptible, since there
is no lasting immunity to guinea worm). Do not install any packages.
```

### Prompt 4: Calibration Against Real Carter Center Data
**Tests:** Multi-objective calibration with two host species using real surveillance data

```
Design a calibration framework for the guinea worm LASER model using real
surveillance data from the Carter Center's Guinea Worm Eradication Program
in Chad:

Observed annual data (human cases / dog infections):
- 2018: 17 / 1,040
- 2019: 48 / 1,935
- 2020: 12 / 1,507
- 2021: 14 / 830
- 2022: 13 / 688
- 2023: 6 / 526

The calibration should:
1. Match total annual human cases and dog infections using log-scale
   comparison (appropriate for small counts)
2. Reproduce the ~50:1 to 100:1 dog-to-human infection ratio
3. Capture the declining trend in both species reflecting intervention
   scale-up over this period
4. Account for spatial concentration: cases cluster in 3-4 of 8 districts
   (Chari-Baguirmi, Moyen-Chari, Salamat, Mandoul)

Sample from parameter ranges:
- beta_human: U(0.0005, 0.005)
- beta_dog: U(0.005, 0.05)
- cross_species_coupling: U(0.05, 0.5)
- seasonal_amplitude: U(0.3, 0.8)
- containment_efficacy: U(0.5, 0.9)

Use a multi-objective fitness function with appropriate weighting for the
different targets. Write the calibration loop and scoring functions using
LASER framework conventions. Do not install any packages.
```

### Prompt 5: Full Integration (Hardest)
**Tests:** End-to-end assembly with SEIS dynamics, dual-host, real geography, all interventions

```
Build a complete spatial guinea worm transmission model for Chad using the
LASER framework with these specifications:

Disease dynamics:
- SEIS (no lasting immunity — guinea worm confers no acquired immunity,
  recovered individuals return to susceptible)
- Humans: R0≈1.5, 365-day pre-patent period, 21-day worm emergence
- Dogs: R0≈2.5, 365-day pre-patent period, 21-day worm emergence
- Water-source mediated transmission (not direct contact — larvae released
  into water, ingested via copepods)
- Cross-species coupling via shared water sources (asymmetric: dogs are
  the dominant reservoir, ~95% of infections)

Spatial structure:
- 8 districts in the endemic corridor with heterogeneous populations
  (human/dog): Chari-Baguirmi (450k/40k), Moyen-Chari (400k/35k),
  Salamat (350k/30k), Mandoul (300k/25k), Logone Oriental (250k/20k),
  Tandjile (200k/15k), Mayo-Kebbi Est (150k/10k), Lac (100k/5k)
- Gravity-model spatial coupling representing shared water source networks

Seasonality:
- Dry-season peak (Jan-Mar at 1.8x, rainy season Jun-Sep at 0.4x)

Interventions (all prevention-based, no vaccine exists):
- ABATE larvicide: starts at 30% coverage year 1, increases 10%/year to
  max 90%
- Filter distribution: starts at 40%, increases 10%/year to max 90%
- Case containment: 70% efficacy, improving 3%/year
- Dog tethering: 40% coverage, increasing 5%/year

Initial conditions (near-elimination context):
- Seed 5 infected humans and 200 infected dogs in the 3 most endemic
  districts (Chari-Baguirmi, Moyen-Chari, Salamat)
- Importation: 1 infected human traveler every 90 days for first 3 years

Simulation:
- 10-year simulation
- Output: annual human cases and dog infections per district
- Track whether elimination (zero cases for 3 consecutive years) is
  achieved and when

Write complete, runnable Python code. Do not install any packages.
```

## Negative Transfer Risks

Watch for the skill causing guinea worm to be forced into standard patterns:
- **SEIR instead of SEIS** — guinea worm has no lasting immunity
- **Vaccination components** when no vaccine exists
- **School-term seasonality** instead of dry-season peak
- **Gravity migration** without water-source adaptation
- **Measles/polio-like R0** (6-18) instead of guinea worm R0 (1.5-2.5)
- **Short incubation** (days) instead of 365-day pre-patent period
- **Epidemic dynamics** instead of near-elimination context
- **Single-host model** ignoring the dog reservoir
