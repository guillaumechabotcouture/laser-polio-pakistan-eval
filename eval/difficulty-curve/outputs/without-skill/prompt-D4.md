The model ran successfully in ~5.6 seconds. Here's the file: `eastland_seir.py`

**Key design decisions based on the actual LASER API:**

**Model setup:**
- Created a `GeoDataFrame` scenario with 4 patches (100k, 200k, 150k, 80k population) and simple polygon geometries
- Used `initialize_population()` from `laser.core.utils` with fractional initial states: 90% S, 0% E, 1% I, 9% R
- `calc_capacity` is called internally by `Model.__init__()` using the `birthrates` argument to pre-allocate agent storage for births over the 10-year simulation

**Disease parameters:**
- `beta = R0 / mean_infectious_period = 5.0 / 10 = 0.5`
- Latent duration: `dists.gamma(4, 1)` — Numba-wrapped gamma distribution (mean=4 days)
- Infectious duration: `dists.gamma(5, 2)` — Numba-wrapped gamma distribution (mean=10 days)

**Components in correct SEIR order:**
1. `SEIR.Susceptible` — initializes agents/node counts for S
2. `SEIR.Exposed` — E→I transitions (decrement etimer, set itimer)
3. `SEIR.Infectious` — I→R transitions (decrement itimer)
4. `SEIR.Recovered` — R bookkeeping (permanent immunity)
5. `SEIR.Transmission` — S→E based on force of infection with spatial coupling
6. `BirthsByCBR` — births at CBR=30/1000/yr (newborns enter S)
7. `MortalityByCDR` — deaths at CDR=10/1000/yr (removes from all compartments)

**Results (10-year simulation):**
- 470,504 total infections
- 173,620 births, 59,474 deaths
- Final state: 170,872 S / 0 E / 0 I / 473,274 R
