The filesystem is read-only, so I'll deliver the complete code here. Every design decision is annotated inline.

---

```python
#!/usr/bin/env python3
"""
Republic of Eastland — 4-Patch Spatial SEIR Model for Respiratory Illness

Uses the LASER framework (laser-generic) for spatially-explicit
agent-based simulation over 1 year (365 ticks, 1 tick = 1 day).

Specification:
    Patches    : 4 districts — Eastport (100k), Northton (200k),
                 Centerburg (150k), Westham (80k)
    Disease    : R0~5, beta = R0 / mean_infectious_period = 5/10 = 0.5
    Latent     : gamma(shape=4, scale=1)  — mean = 4 days
    Infectious : gamma(shape=5, scale=2)  — mean = 10 days
    Init       : 90% S, 0% E, 1% I, 9% R
    Components : Susceptible → Exposed → Infectious → Recovered → Transmission
    Network    : Gravity model (k=0.01, a=1, b=1, c=1.5), row-normalised ≤15%
    Duration   : 1 year (365 ticks)
    Vital dyn. : None (fixed population)

Usage:
    python3 eastland_respiratory_model.py
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.core.migration import gravity, row_normalizer, distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================================
# 1. Patch Configuration — 4 districts of the Republic of Eastland
# ============================================================================

patches = [
    {"nodeid": 0, "name": "Eastport",   "population": 100_000, "lat": 30.0, "lon": 50.0},
    {"nodeid": 1, "name": "Northton",   "population": 200_000, "lat": 30.5, "lon": 50.5},
    {"nodeid": 2, "name": "Centerburg", "population": 150_000, "lat": 31.0, "lon": 51.0},
    {"nodeid": 3, "name": "Westham",    "population":  80_000, "lat": 31.5, "lon": 51.5},
]

geometry = [Point(p["lon"], p["lat"]) for p in patches]
scenario = gpd.GeoDataFrame(patches, geometry=geometry, crs="EPSG:4326")

populations = np.array(scenario["population"])
nnodes = len(scenario)
names = list(scenario["name"])


# ============================================================================
# 2. Initial Conditions — 90% S, 0% E, 1% I, 9% R
#    E and I are rounded first; R absorbs the remainder so that
#    S + E + I + R == population exactly in every patch.
# ============================================================================

scenario["E"] = np.zeros(nnodes, dtype=np.uint32)
scenario["I"] = np.round(0.01 * populations).astype(np.uint32)
scenario["R"] = np.round(0.09 * populations).astype(np.uint32)
scenario["S"] = (populations
                 - scenario["E"]
                 - scenario["I"]
                 - scenario["R"]).astype(np.uint32)

assert (scenario.S + scenario.E + scenario.I + scenario.R
        == scenario.population).all(), \
    "Initial S+E+I+R must equal population in every patch"
assert (scenario.I > 0).any(), "At least one patch must have initial infections"

print("Initial conditions:")
for _, row in scenario.iterrows():
    print(f"  {row['name']:<12s}: pop={row['population']:>7,}  "
          f"S={row['S']:>6,} ({row['S'] / row['population'] * 100:.0f}%)  "
          f"E={row['E']:>4,}  "
          f"I={row['I']:>4,} ({row['I'] / row['population'] * 100:.0f}%)  "
          f"R={row['R']:>5,} ({row['R'] / row['population'] * 100:.0f}%)")


# ============================================================================
# 3. Parameters
#    beta = R0 / mean_infectious_period = 5 / (5×2) = 0.5
# ============================================================================

NTICKS = 365           # 1 year (1 tick = 1 day)
R0_VALUE = 5.0
MEAN_INF_PERIOD = 5 * 2    # gamma(shape=5, scale=2) → mean = 10 days
BETA = R0_VALUE / MEAN_INF_PERIOD   # 0.5

params = PropertySet({
    "prng_seed": 42,
    "nticks": NTICKS,
    "beta": BETA,
    "capacity_safety_factor": 1.1,   # small buffer; no births in this model
})

print(f"\nParameters: R0={R0_VALUE}, beta={BETA:.3f}")
print(f"  Latent period    : gamma(shape=4, scale=1) → mean = 4 days")
print(f"  Infectious period: gamma(shape=5, scale=2) → mean = 10 days")
print(f"  Duration         : {NTICKS} ticks (1 year)")


# ============================================================================
# 4. Duration Distributions
# ============================================================================

expdurdist = dists.gamma(shape=4, scale=1.0)   # Latent period    — mean = 4 days
infdurdist = dists.gamma(shape=5, scale=2.0)   # Infectious period — mean = 10 days


# ============================================================================
# 5. Build LASER Model (no vital dynamics — fixed population for 1 year)
# ============================================================================

model = Model(scenario, params, birthrates=None)


# ============================================================================
# 6. Gravity Migration Network
#    M_{i,j} = k · p_i^a · p_j^b / d_{ij}^c  (k=0.01, a=1, b=1, c=1.5)
#    Row-normalised so no patch exports more than 15% of its FOI.
# ============================================================================

lats = np.array([p["lat"] for p in patches])
lons = np.array([p["lon"] for p in patches])

dist_matrix = np.zeros((nnodes, nnodes))
for i in range(nnodes):
    for j in range(nnodes):
        if i != j:
            dist_matrix[i, j] = distance(lats[i], lons[i], lats[j], lons[j])

network = gravity(populations.astype(np.float64), dist_matrix, k=0.01, a=1, b=1, c=1.5)
network = row_normalizer(network, max_fraction=0.15)
model.network = network

# Validate network
assert model.network.sum() > 0, \
    "Gravity network is all zeros — check gravity params"
row_sums = model.network.sum(axis=1)
assert row_sums.max() < 0.3, \
    f"Max row sum {row_sums.max():.4f} exceeds safety limit"

print("\nGravity network (k=0.01, a=1, b=1, c=1.5, max_export=15%):")
for i in range(nnodes):
    print(f"  {names[i]:<12s}: row_sum={row_sums[i]:.5f}")


# ============================================================================
# 7. Assemble Components in Correct Order
#
#    Susceptible and Recovered propagate compartment counts tick-to-tick,
#    preserving S + E + I + R = N at all times.
#
#    Ordering rationale:
#      1. Susceptible  — copy S[t] → S[t+1] before Transmission decrements it
#      2. Exposed      — decrement etimer; transition E→I when timer expires
#      3. Infectious   — decrement itimer; transition I→R when timer expires
#      4. Recovered    — copy R[t] → R[t+1] after Infectious has incremented it
#      5. Transmission — compute spatially-coupled FOI; draw Bernoulli trials;
#                        decrement S[t+1], increment E[t+1], assign etimer
# ============================================================================

model.components = [
    SEIR.Susceptible(model),                      # 1. Propagate S counts
    SEIR.Exposed(model, expdurdist, infdurdist),  # 2. E→I transitions
    SEIR.Infectious(model, infdurdist),           # 3. I→R transitions
    SEIR.Recovered(model),                        # 4. Propagate R counts
    SEIR.Transmission(model, expdurdist),         # 5. S→E (spatially-coupled FOI)
]


# ============================================================================
# 8. Run Simulation
# ============================================================================

total_pop = int(scenario.population.sum())
print(f"\nRunning 1-year Eastland respiratory SEIR simulation ...")
print(f"  Patches          : {nnodes} ({', '.join(names)})")
print(f"  Total population : {total_pop:,}")
print(f"  Agent capacity   : {model.people.capacity:,}")

model.run("Eastland Respiratory SEIR")
print("Simulation complete.")


# ============================================================================
# 9. Post-Run Verification
# ============================================================================

print("\n--- Verification ---")

# Compartment non-negativity
any_negative = False
for arr, cname in [(model.nodes.S, "S"), (model.nodes.E, "E"),
                   (model.nodes.I, "I"), (model.nodes.R, "R")]:
    if np.any(arr[:NTICKS] < 0):
        print(f"  WARNING: Negative {cname} values detected!")
        any_negative = True
if not any_negative:
    print("  All compartments non-negative: OK")

# Population invariance at sampled ticks
invariant_ok = True
for tick in [0, NTICKS // 4, NTICKS // 2, 3 * NTICKS // 4, NTICKS - 1]:
    total_t = (int(model.nodes.S[tick].sum()) + int(model.nodes.E[tick].sum())
               + int(model.nodes.I[tick].sum()) + int(model.nodes.R[tick].sum()))
    if total_t != total_pop:
        print(f"  WARNING: S+E+I+R={total_t} != {total_pop} at tick {tick}")
        invariant_ok = False
if invariant_ok:
    print(f"  Population invariant (S+E+I+R={total_pop:,}): OK at all sampled ticks")

# Epidemic dynamics
total_infections = int(model.nodes.newly_infected[:NTICKS].sum())
print(f"  Total infections (1 year): {total_infections:,} "
      f"({total_infections / total_pop * 100:.1f}% attack rate)")

patches_with_cases = int(np.sum(model.nodes.newly_infected[:NTICKS].sum(axis=0) > 0))
print(f"  Patches with infections: {patches_with_cases}/{nnodes}")

print(f"  Network sum: {model.network.sum():.4f}, max row: {row_sums.max():.4f}: OK")

# Compartment table
print(f"\nCompartment counts (sum across all patches):")
print(f"  {'Tick':<8} {'S':>10} {'E':>8} {'I':>8} {'R':>10} {'Total':>10}")
print(f"  {'-' * 57}")
for tick in [0, NTICKS // 4, NTICKS // 2, 3 * NTICKS // 4, NTICKS - 1]:
    S = int(model.nodes.S[tick].sum())
    E = int(model.nodes.E[tick].sum())
    I = int(model.nodes.I[tick].sum())
    R = int(model.nodes.R[tick].sum())
    print(f"  {tick:<8d} {S:>10,} {E:>8,} {I:>8,} {R:>10,} {S + E + I + R:>10,}")

# Final state by patch with attack rates
print(f"\nFinal state by patch (tick {NTICKS - 1}):")
print(f"  {'Patch':<12} {'Pop':>8} {'S':>8} {'E':>5} {'I':>5} {'R':>8} {'Attack%':>8}")
print(f"  {'-' * 58}")
for i in range(nnodes):
    S_f = int(model.nodes.S[NTICKS - 1, i])
    E_f = int(model.nodes.E[NTICKS - 1, i])
    I_f = int(model.nodes.I[NTICKS - 1, i])
    R_f = int(model.nodes.R[NTICKS - 1, i])
    init_S = int(scenario.S.iloc[i])
    init_R = int(scenario.R.iloc[i])
    new_R = R_f - init_R
    attack_pct = new_R / init_S * 100 if init_S > 0 else 0.0
    print(f"  {names[i]:<12} {int(populations[i]):>8,} "
          f"{S_f:>8,} {E_f:>5,} {I_f:>5,} {R_f:>8,} {attack_pct:>7.1f}%")


# ============================================================================
# 10. Diagnostic Plots
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
colors = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12"]
days_axis = np.arange(NTICKS)

# (A) Daily incidence by patch — 7-day centred rolling average
ax = axes[0, 0]
kernel = np.ones(7) / 7
for i in range(nnodes):
    daily = model.nodes.newly_infected[:NTICKS, i].astype(np.float64)
    smoothed = np.convolve(daily, kernel, mode="same")
    ax.plot(days_axis, smoothed, color=colors[i], lw=1.8,
            label=f"{names[i]} ({int(populations[i] / 1000)}k)", alpha=0.9)
ax.set_xlabel("Day")
ax.set_ylabel("New infections (7-day avg)")
ax.set_title("(A) Daily Incidence by Patch")
ax.legend(fontsize=9)
ax.set_xlim(0, NTICKS - 1)

# (B) Susceptible fraction S/N over time
ax = axes[0, 1]
for i in range(nnodes):
    S = model.nodes.S[:NTICKS, i].astype(np.float64)
    N = np.maximum(S + model.nodes.E[:NTICKS, i]
                   + model.nodes.I[:NTICKS, i]
                   + model.nodes.R[:NTICKS, i], 1.0)
    ax.plot(days_axis, S / N, color=colors[i], lw=1.8, label=names[i], alpha=0.9)
ax.axhline(1.0 / R0_VALUE, color="black", ls=":", lw=1.5,
           alpha=0.6, label=f"1/R0 = {1.0 / R0_VALUE:.2f}")
ax.set_xlabel("Day")
ax.set_ylabel("S / N")
ax.set_title("(B) Susceptible Fraction")
ax.legend(fontsize=9)
ax.set_xlim(0, NTICKS - 1)
ax.set_ylim(0, 1)

# (C) Stacked SEIR fractions — Northton (largest patch, index 1)
ax = axes[1, 0]
i_large = 1
pop_i = populations[i_large]
ax.stackplot(
    days_axis,
    model.nodes.S[:NTICKS, i_large] / pop_i,
    model.nodes.E[:NTICKS, i_large] / pop_i,
    model.nodes.I[:NTICKS, i_large] / pop_i,
    model.nodes.R[:NTICKS, i_large] / pop_i,
    labels=["S (Susceptible)", "E (Exposed)", "I (Infectious)", "R (Recovered)"],
    colors=["#3498db", "#f39c12", "#e74c3c", "#2ecc71"],
    alpha=0.85,
)
ax.set_xlabel("Day")
ax.set_ylabel("Fraction of population")
ax.set_title(f"(C) SEIR Compartments — {names[i_large]} ({int(pop_i / 1000)}k)")
ax.legend(loc="upper right", fontsize=9)
ax.set_xlim(0, NTICKS - 1)
ax.set_ylim(0, 1)

# (D) Effective reproduction number R_eff = R0 × S/N
ax = axes[1, 1]
for i in range(nnodes):
    S = model.nodes.S[:NTICKS, i].astype(np.float64)
    N = np.maximum(S + model.nodes.E[:NTICKS, i]
                   + model.nodes.I[:NTICKS, i]
                   + model.nodes.R[:NTICKS, i], 1.0)
    ax.plot(days_axis, R0_VALUE * S / N, color=colors[i], lw=1.8,
            label=names[i], alpha=0.9)
ax.axhline(1.0, color="black", ls="-", lw=1.5, alpha=0.4,
           label="R_eff = 1 (epidemic threshold)")
ax.set_xlabel("Day")
ax.set_ylabel("R_eff = R0 × S/N")
ax.set_title("(D) Effective Reproduction Number")
ax.legend(fontsize=9)
ax.set_xlim(0, NTICKS - 1)

plt.suptitle(
    "Republic of Eastland — Respiratory SEIR  (R0=5, 1 year)\n"
    "latent=gamma(4,1)·mean=4d,  infectious=gamma(5,2)·mean=10d,  "
    "spatial gravity network",
    fontsize=11, fontweight="bold",
)
plt.tight_layout()

outdir = Path(__file__).parent / "eval" / "outputs"
outdir.mkdir(parents=True, exist_ok=True)
outpath = outdir / "eastland_respiratory_seir.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nDiagnostic plot saved to {outpath}")
plt.close(fig)
```

---

**Key design decisions explained:**

**Beta derivation** — `beta = R0 / mean_infectious_period = 5 / 10 = 0.5`. The mean of a Gamma(shape=5, scale=2) is shape×scale = 10 days.

**Component ordering** — The 5 components must execute in this exact sequence each tick:
1. `Susceptible` — copies `S[t] → S[t+1]` so Transmission has a value to decrement
2. `Exposed` — decrements `etimer`; when it hits 0, moves agent to Infectious and assigns `itimer`
3. `Infectious` — decrements `itimer`; when it hits 0, moves agent to Recovered
4. `Recovered` — copies `R[t] → R[t+1]` after Infectious has already incremented it
5. `Transmission` — computes spatially-coupled FOI via `model.network`, runs Bernoulli trials on susceptibles, decrements `S[t+1]` and increments `E[t+1]`

**Gravity network** — Built manually with `gravity()` + `row_normalizer()` so the network can be validated and inspected before the run. The network is checked (`sum > 0`, `max_row < 0.3`) before `model.run()`.

**No vital dynamics** — For a 1-year simulation with no births/deaths requested, `birthrates=None` is passed to `Model`. This keeps the population fixed.

**Distributions** — `dists.gamma()` returns a Numba-compatible callable used inside the LASER transmission loop. `expdurdist` is passed to both `Exposed` (assigns `etimer` on exposure) and `Transmission` (same distribution, used when S→E occurs in the transmission step).
