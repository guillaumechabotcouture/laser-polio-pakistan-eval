The model runs successfully. Here's a summary of `eastland_seir.py`:

**Architecture:**
- 4 patches of Eastland (100k, 200k, 150k, 80k population) defined as a `GeoDataFrame` with rectangular geometries
- Gravity-based spatial coupling between patches computed automatically by the `Model`

**Disease parameters:**
- `beta = R0 / mean_infectious = 5 / 10 = 0.5`
- Latent period: `dists.gamma(4, 1)` — mean 4 days
- Infectious period: `dists.gamma(5, 2)` — mean 10 days

**SEIR components in correct order:**
1. `Susceptible` — initializes agents and S compartment
2. `Transmission` (TransmissionSE) — computes force of infection, S→E transitions; uses `inf_dist` as its `expdurdist` per your specification
3. `Exposed` — decrements incubation timers, E→I transitions
4. `Infectious` (InfectiousIR) — decrements infectious timers, I→R transitions
5. `Recovered` — tracks the R compartment

**Results (365-day run):**
- Peak infectious: **120,021** at day 48
- Total infections: **468,602** out of 530,000
- Final state: ~8,400 S / 0 E / 0 I / ~521,600 R — epidemic burns through the susceptible population and resolves within a year, consistent with R0=5
