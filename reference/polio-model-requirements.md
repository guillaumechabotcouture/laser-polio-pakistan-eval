# Polio Model Domain Requirements

Reference specification extracted from existing `polio-seir-model/` and domain knowledge.
Used as ground truth for evaluating generated models.

## Disease Biology

- **Pathogen**: Poliovirus (WPV1 focus for Pakistan)
- **Infectious period**: ~28 days (gamma = 1/28)
- **Latent period**: ~3 days (sigma = 1/3)
- **R0**: ~6 (beta_scale ≈ 0.2143 for daily timestep)
- **Paralysis rate**: ~1 in 200 infections are paralytic (AFP)
- **Asymptomatic transmission**: ~99% of infections are subclinical but still shed virus

## Compartmental Structure

A basic SEIR is insufficient for polio. A realistic model needs:

### Immunity states (at minimum)
- Naive (never infected/vaccinated)
- Recently immune (strong protection, post-infection or post-vaccination)
- Waned (partial protection, reduced susceptibility)
- Partially susceptible (can be re-infected but with lower viral shedding)

### Age structure
- 0-5 years (primary target population, highest susceptibility)
- 5-15 years (school-age, lower susceptibility)
- 15+ years (adults, potential silent transmission reservoir)

### Vaccination
- **OPV (oral polio vaccine)**: Live attenuated, induces mucosal immunity, causes secondary transmission (vaccine virus shedding), risk of cVDPV
- **IPV (inactivated polio vaccine)**: Induces humoral immunity only, no mucosal immunity, no shedding
- **Routine immunization**: Multiple doses at scheduled ages (6, 10, 14 weeks + boosters)
- **SIAs (Supplementary Immunization Activities)**: Mass campaigns targeting 0-5 year olds
- **Reactive vaccination**: Triggered by AFP case detection

## Spatial Structure (Pakistan)

- **154 districts** across 4 provinces + federal territories
- **Core endemic reservoirs**: Quetta block (Balochistan), Peshawar corridor (KP), Karachi
- **Cross-border transmission**: Afghanistan (southern corridor)
- **Population**: ~230 million, high birth rate (~3.5%)
- **Seasonal forcing**: Transmission peaks in monsoon season (Jul-Oct)

## Surveillance Data

- **AFP (Acute Flaccid Paralysis)**: Primary surveillance — reports paralytic cases by district/week
- **Environmental surveillance**: Sewage sampling detects poliovirus in ~50 sites
- **Genetic sequencing**: Tracks transmission chains via VP1 sequencing

## Key Modeling Challenges

1. Asymptomatic bulk of transmission (AFP sees only ~0.5% of infections)
2. OPV dynamics (vaccine-derived immunity AND vaccine-derived transmission)
3. Population movement patterns (nomadic populations, internally displaced persons)
4. Spatial heterogeneity in vaccination coverage
5. Waning immunity requiring booster dynamics
