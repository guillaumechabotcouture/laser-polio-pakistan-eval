The simulation ran successfully in about 6 seconds. Here's the complete file at `eastland_seir.py`. Key results from the 10-year run:

**Model configuration:**
- 4 patches (530K total population) with gravity-coupled spatial transmission
- R0=5, beta=0.5, latent ~4 days (Gamma(4,1)), infectious ~10 days (Gamma(5,2))
- Seasonal forcing oscillating between 1.3x (winter) and 0.7x (summer)
- CBR=30, CDR=10 vital dynamics; pulse importation of 3 cases every 60 ticks

**Results:**
- **Peak infectious**: 211,588 at tick 28 (~day 28) — a rapid initial epidemic given 1% seeded as infectious with R0=5
- **Final population**: 643,491 (up from 530,000 due to net births exceeding deaths)
- **Final state**: ~25% susceptible, ~75% recovered, near-zero E/I — the population is largely in a post-epidemic steady state with periodic re-seeding from pulse importation
- **Total births**: 172,930 | **Total deaths**: 59,439 over the 10-year period

The custom `PulseImportation` component follows the standard LASER `__init__(model)` + `step(tick)` pattern and injects 3 infections into patch 0 every 60 ticks, maintaining ongoing transmission pressure throughout the simulation.
