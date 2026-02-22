Here is the complete code. The `/eval-work` directory is read-only in this environment, so the file is at `/tmp/eastland_seir.py`, but here is the full listing:

```python
"""
Spatial SEIR model – Republic of Eastland
==========================================
4-patch respiratory illness model built with the LASER framework
(laser-generic package, https://laser.idmod.org/).

Patches
-------
  Patch 0:  100,000 people  (30.0°N, 50.0°E)
  Patch 1:  200,000 people  (30.5°N, 50.5°E)
  Patch 2:  150,000 people  (31.0°N, 51.0°E)
  Patch 3:   80,000 people  (31.5°N, 51.5°E)

Disease (respiratory)
---------------------
  R0  ~ 5
  Latent period     : Gamma(shape=4, scale=1)  → mean  4 days
  Infectious period : Gamma(shape=5, scale=2)  → mean 10 days
  β   = R0 / mean_infectious_period = 5 / 10  = 0.50 /day

Initial conditions
------------------
  S=90 %  |  E=0 %  |  I=1 %  |  R=9 %

Vital dynamics
--------------
  CBR = 30 /1000/year  (BirthsByCBR)
  CDR = 10 /1000/year  (MortalityByCDR)

Spatial coupling
----------------
  Gravity model : k=0.01, a=1, b=1, c=1.5
  row_normalizer: max_rowsum=0.15 (≤15 % of FOI exported per patch)

Simulation
----------
  Duration : 10 years (3 650 daily ticks)
  Capacity : sized via calc_capacity (laser-core)
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point

from laser.core import PropertySet
from laser.core.migration import distance, gravity, row_normalizer
from laser.core.utils import calc_capacity
import laser.core.distributions as dists

from laser.generic.model import Model
from laser.generic import SEIR
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR

# ──────────────────────────────────────────────────────────────────────────────
# 1.  Simulation constants
# ──────────────────────────────────────────────────────────────────────────────

NUM_YEARS = 10
NTICKS    = NUM_YEARS * 365          # 3 650 daily time-steps

# Disease: β = R0 / mean_infectious_period
R0_VALUE      = 5.0
MEAN_INF_DAYS = 5 * 2               # gamma(shape=5, scale=2) → mean = 10 days
BETA          = R0_VALUE / MEAN_INF_DAYS   # 0.50 /day

CBR = 30.0   # crude birth rate  [per 1 000 per year]
CDR = 10.0   # crude death rate  [per 1 000 per year]

SEED = 42

# ──────────────────────────────────────────────────────────────────────────────
# 2.  Patch definitions
# ──────────────────────────────────────────────────────────────────────────────

populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int64)
lats        = np.array([30.0, 30.5, 31.0, 31.5])   # degrees North
lons        = np.array([50.0, 50.5, 51.0, 51.5])   # degrees East
num_nodes   = len(populations)

patch_names = [
    "Eastland-A (30°N 50°E)",
    "Eastland-B (30.5°N 50.5°E)",
    "Eastland-C (31°N 51°E)",
    "Eastland-D (31.5°N 51.5°E)",
]

# ──────────────────────────────────────────────────────────────────────────────
# 3.  Initial compartment sizes:  S=90 %, E=0 %, I=1 %, R=9 %
# ──────────────────────────────────────────────────────────────────────────────

I0        = np.round(0.01 * populations).astype(np.int64)
E0        = np.zeros(num_nodes, dtype=np.int64)
R0_counts = np.round(0.09 * populations).astype(np.int64)
S0        = populations - E0 - I0 - R0_counts      # remainder ≈ 90 %

assert np.all(S0 + E0 + I0 + R0_counts == populations), (
    "Compartment counts must sum exactly to population"
)

print("Initial compartment sizes:")
for i, name in enumerate(patch_names):
    print(f"  {name}: N={populations[i]:,}  "
          f"S={S0[i]:,}  E={E0[i]:,}  I={I0[i]:,}  R={R0_counts[i]:,}")

# ──────────────────────────────────────────────────────────────────────────────
# 4.  Scenario GeoDataFrame
#     Required columns: nodeid, population, S, E, I, R, geometry
#     Point(lon, lat) – Shapely uses (x=lon, y=lat) convention
# ──────────────────────────────────────────────────────────────────────────────

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(num_nodes, dtype=np.int32),
        "population": populations,
        "S":          S0,
        "E":          E0,
        "I":          I0,
        "R":          R0_counts,
    },
    geometry=[Point(lon, lat) for lat, lon in zip(lats, lons)],
    crs="EPSG:4326",
)

# ──────────────────────────────────────────────────────────────────────────────
# 5.  Vital-dynamics rate arrays   shape: (nticks, num_nodes)
# ──────────────────────────────────────────────────────────────────────────────

birthrates     = np.full((NTICKS, num_nodes), CBR, dtype=np.float64)
mortalityrates = np.full((NTICKS, num_nodes), CDR, dtype=np.float64)

# ──────────────────────────────────────────────────────────────────────────────
# 6.  Explicit capacity calculation via calc_capacity (10-year horizon)
#     calc_capacity(birthrates, initial_pop, safety_factor) → per-node int array
#     The Model constructor also calls this internally when birthrates is given.
# ──────────────────────────────────────────────────────────────────────────────

capacities = calc_capacity(birthrates, populations, safety_factor=1.0)

print("\nAgent capacities (calc_capacity, 10 years, CBR=30):")
for i, name in enumerate(patch_names):
    print(f"  {name}: {capacities[i]:,}")
print(f"  Total agents allocated: {capacities.sum():,}")

# ──────────────────────────────────────────────────────────────────────────────
# 7.  Model parameters (PropertySet)
# ──────────────────────────────────────────────────────────────────────────────

params = PropertySet({"nticks": NTICKS, "beta": BETA, "seed": SEED})

# ──────────────────────────────────────────────────────────────────────────────
# 8.  Initialise the Model
#     Passing birthrates causes the constructor to call calc_capacity internally
#     to size the people LaserFrame for the full 10-year duration.
# ──────────────────────────────────────────────────────────────────────────────

model = Model(
    scenario,
    params,
    birthrates=birthrates,
    name="eastland_seir",
)

print(f"\nModel initialised: {model.people.count:,} agents, "
      f"{model.nodes.count} patches, {NTICKS:,} ticks ({NUM_YEARS} years)")

# ──────────────────────────────────────────────────────────────────────────────
# 9.  Duration distributions  (Numba-JIT callables: f(tick, node) → float32)
#     Latent period     : Gamma(shape=4, scale=1)   mean = 4 days
#     Infectious period : Gamma(shape=5, scale=2)   mean = 10 days
# ──────────────────────────────────────────────────────────────────────────────

exp_dur = dists.gamma(shape=4, scale=1)   # latent-period distribution
inf_dur = dists.gamma(shape=5, scale=2)   # infectious-period distribution

# ──────────────────────────────────────────────────────────────────────────────
# 10.  Pairwise distances and gravity network
# ──────────────────────────────────────────────────────────────────────────────

# Haversine pairwise distances (km)  → (4 × 4) symmetric matrix
dist_matrix = distance(lats, lons, lats, lons)

print("\nPairwise Haversine distances (km):")
print(np.round(dist_matrix, 1))

# Gravity network:  network[i,j] = k × pop_i^a × pop_j^b / dist_ij^c
#   k=0.01, a=1, b=1, c=1.5  (diagonal forced to 0 by gravity())
net_raw = gravity(
    populations.astype(np.float64),
    dist_matrix,
    k=0.01,
    a=1.0,
    b=1.0,
    c=1.5,
)

print("\nGravity network (raw):")
print(np.round(net_raw, 4))

# Row-normalise: scale rows so no patch exports more than 15 % of its FOI
net_norm = row_normalizer(net_raw, max_rowsum=0.15)

print("\nNormalised gravity network:")
print(np.round(net_norm, 6))
print("Row sums (≤ 0.15):", np.round(net_norm.sum(axis=1), 6))

# Override the Model's default gravity network with the custom normalised one
model.network = net_norm

# ──────────────────────────────────────────────────────────────────────────────
# 11.  Assemble SEIR components
#
#  Initialisation order (dependency-driven):
#    Susceptible  → assigns all agents to patches; default state = SUSCEPTIBLE
#    Exposed      → promotes scenario.E agents per patch: S → E (etimer set)
#    Infectious   → promotes scenario.I agents per patch: S → I (itimer set)
#    Recovered    → promotes scenario.R agents per patch: S → R (permanent)
#    Transmission → adds FOI tracking arrays; no state changes at init
#    Deaths/Births → vital-dynamics tracking arrays
#
#  Step order (per-tick, called by model.run()):
#    _initialize_flows copies state[t] → state[t+1] for all compartments
#    1. susceptible  – validates S census; no transitions
#    2. transmission – FOI from I[t]; draws S → E stochastically (etimer set)
#    3. exposed      – etimer countdown; E → I on expiry (itimer set)
#    4. infectious   – itimer countdown; I → R on expiry
#    5. recovered    – validates R census; permanent immunity
#    6. deaths       – stochastic CDR mortality; marks agents DECEASED
#    7. births       – Poisson CBR births; new agents appended as S
# ──────────────────────────────────────────────────────────────────────────────

# Initialise in dependency order (Susceptible must be first) ------------------

# Adds people.nodeid, people.state (all=SUSCEPTIBLE), nodes.S
susceptible  = SEIR.Susceptible(model)

# Exposed.__init__: converts scenario.E agents S→E; adds people.etimer, nodes.E
exposed      = SEIR.Exposed(model, exp_dur, inf_dur)

# InfectiousIR.__init__: converts scenario.I agents S→I; adds people.itimer, nodes.I
infectious   = SEIR.Infectious(model, inf_dur)

# Recovered.__init__: converts scenario.R agents S→R; adds nodes.R (no timers)
recovered    = SEIR.Recovered(model)

# TransmissionSE.__init__: adds nodes.forces, nodes.newly_infected; no state changes
transmission = SEIR.Transmission(model, exp_dur)

# Vital dynamics
deaths = MortalityByCDR(model, mortalityrates)

# track=False: no age-pyramid needed; newborns are still added as S each tick
births = BirthsByCBR(model, birthrates, pyramid=None, track=False)

# Set component step order ----------------------------------------------------
model.components = [
    susceptible,    # 1. validates S census
    transmission,   # 2. S → E  (FOI from I[t] + spatial network coupling)
    exposed,        # 3. E → I  (etimer countdown)
    infectious,     # 4. I → R  (itimer countdown)
    recovered,      # 5. validates R census; permanent immunity
    deaths,         # 6. stochastic CDR mortality
    births,         # 7. Poisson CBR births (new S agents)
]

# ──────────────────────────────────────────────────────────────────────────────
# 12.  Run the simulation for 10 years
# ──────────────────────────────────────────────────────────────────────────────

print(f"\nStarting SEIR simulation: {NTICKS:,} ticks ({NUM_YEARS} years)")
print(f"  β={BETA:.3f}/day  R0≈{R0_VALUE}  "
      f"latent~Gamma(4,1)  infectious~Gamma(5,2)")
print(f"  CBR={CBR}/1000/yr  CDR={CDR}/1000/yr")

model.run()

print("Simulation complete.\n")

# ──────────────────────────────────────────────────────────────────────────────
# 13.  Summary output
# ──────────────────────────────────────────────────────────────────────────────

print("=" * 70)
print("RESULTS – End of year 10 (final tick)")
print("=" * 70)

for state, label in [("S", "Susceptible"), ("E", "Exposed"),
                     ("I", "Infectious"),  ("R", "Recovered")]:
    if hasattr(model.nodes, state):
        vals = getattr(model.nodes, state)[-1]
        print(f"\n  {label}:")
        for i, name in enumerate(patch_names):
            print(f"    {name:30s}: {vals[i]:>10,}")
        print(f"    {'Total':30s}: {vals.sum():>10,}")

print("\n  Cumulative vital dynamics (10 years):")
if hasattr(model.nodes, "deaths"):
    cd = model.nodes.deaths.sum(axis=0)
    print(f"    Deaths  by patch: {list(cd)}  → total {cd.sum():,}")
if hasattr(model.nodes, "births"):
    cb = model.nodes.births.sum(axis=0)
    print(f"    Births  by patch: {list(cb)}  → total {cb.sum():,}")

print("\n  Epidemic incidence summary:")
if hasattr(model.nodes, "newly_infected"):
    inc = model.nodes.newly_infected          # shape: (nticks+1, num_nodes)
    peak_inc  = inc.max(axis=0)
    peak_day  = inc.argmax(axis=0)
    total_inc = inc.sum(axis=0)
    print(f"    {'Patch':30s}  {'Peak/day':>10}  {'Peak day':>9}  "
          f"{'Total infected':>14}")
    for i, name in enumerate(patch_names):
        print(f"    {name:30s}  {peak_inc[i]:>10,}  "
              f"{peak_day[i]:>9}  {total_inc[i]:>14,}")
    print(f"    {'All patches':30s}  {peak_inc.sum():>10,}  "
          f"{'':>9}  {total_inc.sum():>14,}")

print("\nDone.")
```

---

**Design decisions and key facts:**

**β derivation** — The infectious-period distribution is Gamma(5, 2), giving mean = shape × scale = 10 days. With R0 = 5, β = R0 / mean_inf = **0.5 /day**.

**`calc_capacity`** — Called in step 6 explicitly for visibility, and again internally by the `Model` constructor (which receives `birthrates`). It applies geometric-Brownian-motion growth over 3 650 ticks with CBR = 30, yielding ~34 % more agent slots than the initial population.

**Component initialisation order** matters because each component's `__init__` converts agents from the default SUSCEPTIBLE state:
1. `Susceptible` — assigns all agents to nodes, sets default state = S
2. `Exposed` — converts scenario.E agents S→E and assigns `etimer`
3. `Infectious` (InfectiousIR) — converts scenario.I agents S→I and assigns `itimer`
4. `Recovered` — converts scenario.R agents S→R (permanent immunity, no timer)
5. `Transmission` (TransmissionSE) — adds FOI arrays only, no state changes

**Network override** — `Model.__init__` builds its own default gravity network; assigning to `model.network` after construction replaces it with the custom k=0.01 / c=1.5 / row-normalised (≤15 %) matrix before `model.run()` is called.

**`BirthsByCBR(pyramid=None, track=False)`** — `track=False` skips date-of-birth recording, which avoids the need for an `AliasedDistribution` age-pyramid while still adding newborns as susceptible agents each tick.
