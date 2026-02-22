The model runs successfully. Here's a summary of what the code does:

**`eastland_seir.py`** — Spatial SEIR model for the Republic of Eastland

**Setup:**
- **4 patches** with populations 100k, 200k, 150k, and 80k (530k total)
- **Coordinates** in EPSG:4326 with gravity-based spatial coupling between patches
- **Initial state:** 90% Susceptible, 0% Exposed, 1% Infectious, 9% Recovered

**Disease parameters:**
- **beta = 0.5** (derived from R0/infectious_period = 5/10)
- **Latent period:** gamma(shape=4, scale=1) → mean 4 days
- **Infectious period:** gamma(shape=5, scale=2) → mean 10 days

**Components** (in correct SEIR order):
1. `Susceptible` — initializes S compartment
2. `Exposed` — handles E→I transitions via incubation timer
3. `Infectious` — handles I→R transitions via infectious timer
4. `Recovered` — absorbing state (permanent immunity)
5. `Transmission` — computes force of infection with spatial coupling, handles S→E

**Results** (365-day run):
- The epidemic peaked at **177,751 infectious** on **day 32**
- By year's end, 98.4% of the population had been infected (521,779 recovered), with only 8,221 still susceptible — consistent with R0=5 dynamics
