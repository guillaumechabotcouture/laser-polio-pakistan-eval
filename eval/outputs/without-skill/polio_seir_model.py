"""
Spatial SEIR Model for Polio Transmission Dynamics
===================================================
Gates Foundation Polio Eradication Modeling Project

Simulates wild poliovirus type 1 (WPV1) transmission across 10 Pakistani
districts using the LASER framework (laser-generic) with a standard SEIR
compartmental model and gravity-based spatial coupling.

Parameters:
  R0 ≈ 6          (wild poliovirus basic reproduction number)
  Latent period:   3 days mean (gamma-distributed)
  Infectious:      28 days mean (gamma-distributed)
  Immunity:        95% recovered at baseline per patch
  Seed cases:      5 infectious per patch
  Duration:        10 years (3,650 days)
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point

import laser.core.distributions as dists
from laser.core import PropertySet
from laser.generic import Model, SEIR

# ── Configuration ──────────────────────────────────────────────────────

NUM_PATCHES = 10
POP_PER_PATCH = 100_000
TOTAL_POP = NUM_PATCHES * POP_PER_PATCH

# Epidemiological parameters for WPV1
R0 = 6.0
LATENT_PERIOD_DAYS = 3        # Mean incubation (E → I)
INFECTIOUS_PERIOD_DAYS = 28   # Mean shedding duration (I → R)
BETA = R0 / INFECTIOUS_PERIOD_DAYS  # ≈ 0.2143 per day

# Simulation duration
YEARS = 10
NTICKS = 365 * YEARS  # 3,650 days

# Initial conditions per patch
INITIAL_INFECTIOUS = 5
FRACTION_IMMUNE = 0.95
INITIAL_RECOVERED = int(POP_PER_PATCH * FRACTION_IMMUNE)  # 95,000
INITIAL_SUSCEPTIBLE = POP_PER_PATCH - INITIAL_RECOVERED - INITIAL_INFECTIOUS  # 4,995
INITIAL_EXPOSED = 0

SEED = 20260220  # PRNG seed for reproducibility

# ── Scenario: 10 Pakistani districts (polio high-risk) ────────────────
# Approximate geographic coordinates for gravity-model spatial coupling.

district_data = [
    ("Karachi",          24.86, 67.01),
    ("Quetta",           30.18, 66.99),
    ("Peshawar",         34.01, 71.58),
    ("Lahore",           31.55, 74.35),
    ("Rawalpindi",       33.60, 73.05),
    ("Bannu",            32.99, 70.60),
    ("Tank",             32.22, 70.38),
    ("Dera Ismail Khan", 31.83, 70.90),
    ("Chaman",           30.92, 66.45),
    ("Zhob",             31.35, 69.45),
]

names = [d[0] for d in district_data]

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     list(range(NUM_PATCHES)),
        "population": [POP_PER_PATCH] * NUM_PATCHES,
        "S":          [INITIAL_SUSCEPTIBLE] * NUM_PATCHES,
        "E":          [INITIAL_EXPOSED] * NUM_PATCHES,
        "I":          [INITIAL_INFECTIOUS] * NUM_PATCHES,
        "R":          [INITIAL_RECOVERED] * NUM_PATCHES,
        "geometry":   [Point(lon, lat) for _, lat, lon in district_data],
    },
    index=names,
)

print("District scenario (initial conditions):")
print(scenario[["population", "S", "E", "I", "R"]].to_string())
print()

# ── Model parameters ──────────────────────────────────────────────────

params = PropertySet(
    {
        "nticks":    NTICKS,
        "beta":      BETA,
        "gravity_k": 500,   # Gravity model scaling constant
        "gravity_a": 1.0,   # Origin population exponent
        "gravity_b": 1.0,   # Destination population exponent
        "gravity_c": 2.0,   # Distance decay exponent
        "prng_seed": SEED,
    }
)

print(f"Model parameters:")
print(f"  R0                  = {R0}")
print(f"  beta                = {BETA:.4f} /day")
print(f"  Latent period       = {LATENT_PERIOD_DAYS} days (mean)")
print(f"  Infectious period   = {INFECTIOUS_PERIOD_DAYS} days (mean)")
print(f"  Simulation duration = {YEARS} years ({NTICKS} days)")
print(f"  Total population    = {TOTAL_POP:,}")
print(f"  Patches             = {NUM_PATCHES}")
print()

# ── Duration distributions (Numba JIT-compiled) ───────────────────────

# Latent period: gamma(shape=3, scale=1) → mean = 3 days, CV ≈ 0.58
exp_dur_dist = dists.gamma(shape=3.0, scale=1.0)

# Infectious period: gamma(shape=14, scale=2) → mean = 28 days, CV ≈ 0.27
inf_dur_dist = dists.gamma(shape=14.0, scale=2.0)

# ── Build model ───────────────────────────────────────────────────────

print("Building model...")
model = Model(scenario, params, name="pakistan-polio-seir")

# Register SEIR components (order matters for initialization)
susceptible  = SEIR.Susceptible(model)
exposed      = SEIR.Exposed(model, exp_dur_dist, inf_dur_dist)
infectious   = SEIR.Infectious(model, inf_dur_dist)
recovered    = SEIR.Recovered(model)
transmission = SEIR.Transmission(model, exp_dur_dist)

model.components = [susceptible, exposed, infectious, recovered, transmission]

print(f"  Nodes (patches):  {model.nodes.count}")
print(f"  Agents (people):  {model.people.count:,}")
print()

# ── Run simulation ────────────────────────────────────────────────────

print(f"Running {YEARS}-year simulation...")
model.run("Pakistan Polio SEIR")
print("Simulation complete.\n")

# ── Results ───────────────────────────────────────────────────────────

# Node-level time series: shape (nticks + 1, nnodes)
S = np.array(model.nodes.S)
E = np.array(model.nodes.E)
I = np.array(model.nodes.I)
R = np.array(model.nodes.R)

# Aggregate across all patches
S_total = S.sum(axis=1)
E_total = E.sum(axis=1)
I_total = I.sum(axis=1)
R_total = R.sum(axis=1)

# Incidence
incidence = np.array(model.nodes.newly_infected)
daily_incidence = incidence.sum(axis=1)
cumulative_cases = daily_incidence.cumsum()

# Peak timing
peak_day = daily_incidence.argmax()
peak_year = peak_day // 365 + 1
peak_day_of_year = peak_day % 365

print("=" * 65)
print("  SIMULATION RESULTS SUMMARY")
print("=" * 65)
print(f"  Total population:         {TOTAL_POP:>12,}")
print(f"  Final susceptible:        {S_total[-1]:>12,}")
print(f"  Final exposed:            {E_total[-1]:>12,}")
print(f"  Final infectious:         {I_total[-1]:>12,}")
print(f"  Final recovered/immune:   {R_total[-1]:>12,}")
print(f"  Cumulative infections:    {cumulative_cases[-1]:>12,}")
print(f"  Peak daily incidence:     {daily_incidence.max():>12,}")
print(f"  Peak at day {peak_day} (year {peak_year}, day {peak_day_of_year})")
print()

# Per-district final state
print("Per-district final state:")
print(f"  {'District':<20s} {'S':>8s} {'E':>6s} {'I':>6s} {'R':>10s}")
print(f"  {'-'*20} {'-'*8} {'-'*6} {'-'*6} {'-'*10}")
for idx, name in enumerate(names):
    print(
        f"  {name:<20s} "
        f"{S[-1, idx]:>8,} "
        f"{E[-1, idx]:>6,} "
        f"{I[-1, idx]:>6,} "
        f"{R[-1, idx]:>10,}"
    )
print()

# Annual incidence
print("Annual incidence:")
for year in range(YEARS):
    start = year * 365
    end = min((year + 1) * 365, len(daily_incidence))
    annual = daily_incidence[start:end].sum()
    print(f"  Year {year + 1:>2d}: {annual:>10,} cases")
print()

# Per-district cumulative incidence
district_cumulative = incidence.sum(axis=0)
print("Cumulative incidence by district:")
for idx, name in enumerate(names):
    print(f"  {name:<20s} {district_cumulative[idx]:>10,}")
