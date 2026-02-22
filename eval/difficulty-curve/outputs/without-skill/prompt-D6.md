The file is written. Here's what `eastland_seir.py` does:

**Scenario** - 4 patches as small polygons in EPSG:4326 centered at the requested coordinates, with populations 100k/200k/150k/80k and initial states 90% S, 0% E, 1% I, 9% R via `initialize_population`.

**Disease parameters** - `beta = R0 / mean_infectious = 5/10 = 0.5`, where mean infectious period = gamma(shape=5, scale=2) = 10 days. Latent and infectious durations use `dists.gamma(4,1)` and `dists.gamma(5,2)` (Numba-compiled callables).

**Vital dynamics** - `ValuesMap.from_scalar(30, nticks, 4)` for CBR=30 and `ValuesMap.from_scalar(10, nticks, 4)` for CDR=10, both constant across patches and time.

**Model creation** - `Model(scenario, params, birthrates=birthrates)` uses `calc_capacity` internally to size the agent pool for 10 years of population growth.

**Spatial network** - After the model computes centroids, the script explicitly computes the Haversine distance matrix, builds a gravity network with `k=0.01, a=1, b=1, c=1.5`, applies `row_normalizer(network, 0.15)` to cap FOI export at 15%, and sets `model.network`.

**Components** (in order): `Susceptible` → `TransmissionSE` (S→E) → `Exposed` (E→I) → `InfectiousIR` (I→R) → `Recovered` → `MortalityByCDR` → `BirthsByCBR(track=False)`.

**Run** - `model.run()` iterates 3,650 daily ticks with a summary table printed at yearly intervals.
