The complete code has been written to `/tmp/eastland_seir.py`. Here is the full file:

---

```python
"""Spatial SEIR model for the Republic of Eastland using the LASER framework.

Spatial setup — 4 patches:
  Patch 0  Eastburg   100,000 people
  Patch 1  Westfield  200,000 people
  Patch 2  Northvale  150,000 people
  Patch 3  Southport   80,000 people

Disease parameters:
  R0 ~ 5
  Latent period:     Gamma(shape=4, scale=1)  → mean = 4 days
  Infectious period: Gamma(shape=5, scale=2)  → mean = 10 days
  beta = R0 / mean_infectious_period = 5 / 10 = 0.5 d⁻¹

Initial conditions (uniform across all patches):
  Susceptible (S):  90%
  Exposed (E):       0%
  Infectious (I):    1%
  Recovered (R):     9%

Simulation: 365 daily time steps (1 year)

Component order for SEIR model:
  [Susceptible, Exposed, Infectious, Recovered, Transmission]

  1. Susceptible  — initialises S counts and assigns agents to nodes
  2. Exposed      — initialises E counts; drives E→I transitions each tick
  3. Infectious   — initialises I counts; drives I→R transitions each tick
  4. Recovered    — initialises R counts (permanent immunity, no waning)
  5. Transmission — computes force-of-infection; drives S→E transitions each tick
"""

import numpy as np
import matplotlib.pyplot as plt
import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.utils import grid
from laser.generic import SEIR, Model


# ── Disease Parameters ────────────────────────────────────────────────────────

R0               = 5.0

# Latent (exposed) period: Gamma(shape=4, scale=1) → mean = 4 × 1 = 4 days
LATENT_SHAPE     = 4
LATENT_SCALE     = 1

# Infectious period: Gamma(shape=5, scale=2) → mean = 5 × 2 = 10 days
INFECTIOUS_SHAPE = 5
INFECTIOUS_SCALE = 2
INFECTIOUS_MEAN  = INFECTIOUS_SHAPE * INFECTIOUS_SCALE   # 10 days

# Transmission rate derived from R0 and mean infectious period
BETA = R0 / INFECTIOUS_MEAN      # 5 / 10 = 0.5 per day

NTICKS = 365                     # 1 year of daily time steps


# ── Spatial Setup: Republic of Eastland ───────────────────────────────────────

PATCH_NAMES = ["Eastburg", "Westfield", "Northvale", "Southport"]
POPULATIONS  = [100_000,   200_000,     150_000,     80_000]
N_PATCHES    = len(POPULATIONS)


# ── Initial Conditions ────────────────────────────────────────────────────────

FRAC_S = 0.90   # fraction susceptible
FRAC_E = 0.00   # fraction exposed
FRAC_I = 0.01   # fraction infectious
FRAC_R = 0.09   # fraction recovered

assert abs(FRAC_S + FRAC_E + FRAC_I + FRAC_R - 1.0) < 1e-9


# ── Build Scenario GeoDataFrame ───────────────────────────────────────────────
# laser.core.utils.grid() creates an M×N spatial grid of rectangular patches.
# population_fn(row, col) is called for each node to set its population.

scenario = grid(
    M=1,
    N=N_PATCHES,
    node_size_degs=0.5,
    population_fn=lambda row, col: POPULATIONS[col],
    origin_x=25.0,   # longitude °E (fictional Eastland coordinates)
    origin_y=48.0,   # latitude  °N (fictional)
)

pops = np.array(POPULATIONS, dtype=np.int64)

# Compartment counts; derive S to enforce exact population conservation
init_I = np.round(FRAC_I * pops).astype(int)
init_R = np.round(FRAC_R * pops).astype(int)
init_E = np.round(FRAC_E * pops).astype(int)
init_S = pops - init_E - init_I - init_R

scenario["S"] = init_S
scenario["E"] = init_E
scenario["I"] = init_I
scenario["R"] = init_R


# ── Model Parameters ──────────────────────────────────────────────────────────

params = PropertySet({"nticks": NTICKS, "beta": BETA})


# ── Duration Distributions ────────────────────────────────────────────────────
# dists.gamma() returns a Numba-compatible callable f(tick, node_id) -> float

expdist = dists.gamma(shape=LATENT_SHAPE,     scale=LATENT_SCALE)
infdist = dists.gamma(shape=INFECTIOUS_SHAPE, scale=INFECTIOUS_SCALE)


# ── Assemble Model ────────────────────────────────────────────────────────────

model = Model(scenario, params, name="Eastland-SEIR")

# Components in the correct SEIR order:
#   Susceptible → Exposed → Infectious → Recovered → Transmission
#
# SEIR.Exposed(model, expdurdist, infdurdist)
#   expdurdist: samples incubation timers for agents newly moved S→E
#   infdurdist: samples infectious timers when the E→I transition fires
#
# SEIR.Infectious(model, infdurdist)
#   infdurdist: samples infectious timers for initially infectious agents
#
# SEIR.Transmission(model, expdurdist)
#   expdurdist: samples incubation timers for agents newly moved S→E

s  = SEIR.Susceptible(model)
e  = SEIR.Exposed(model, expdist, infdist)
i  = SEIR.Infectious(model, infdist)
r  = SEIR.Recovered(model)
tx = SEIR.Transmission(model, expdist)

model.components = [s, e, i, r, tx]


# ── Run Simulation ────────────────────────────────────────────────────────────

print("Running 1-year SEIR simulation ...")
model.run("Eastland SEIR")


# ── Results ───────────────────────────────────────────────────────────────────
# model.nodes.S shape: (nticks+1, n_nodes); index -1 = final time step (day 365)

print()
for j, name in enumerate(PATCH_NAMES):
    S_f = int(model.nodes.S[-1, j])
    E_f = int(model.nodes.E[-1, j])
    I_f = int(model.nodes.I[-1, j])
    R_f = int(model.nodes.R[-1, j])
    attack = 100.0 * (init_S[j] - S_f) / pops[j]
    print(f"  {name}: S={S_f:,}  E={E_f:,}  I={I_f:,}  R={R_f:,}  attack={attack:.1f}%")


# ── Plot: SEIR Dynamics per Patch ─────────────────────────────────────────────

ticks = np.arange(NTICKS + 1)   # 0, 1, ..., 365

fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
fig.suptitle(
    "Republic of Eastland — Spatial SEIR Model\n"
    f"R0={R0}  ·  Latent: Gamma(4,1)  ·  Infectious: Gamma(5,2)  ·  1 year",
    fontsize=12, fontweight="bold",
)

COLORS = {"S": "#2196F3", "E": "#FF9800", "I": "#F44336", "R": "#4CAF50"}

for j, (ax, name) in enumerate(zip(axes.flat, PATCH_NAMES)):
    pop_j = int(pops[j])

    def pct(arr, p=pop_j):
        return arr / p * 100

    ax.plot(ticks, pct(model.nodes.S[:, j]), color=COLORS["S"], lw=2, label="S — Susceptible")
    ax.plot(ticks, pct(model.nodes.E[:, j]), color=COLORS["E"], lw=2, label="E — Exposed")
    ax.plot(ticks, pct(model.nodes.I[:, j]), color=COLORS["I"], lw=2, label="I — Infectious")
    ax.plot(ticks, pct(model.nodes.R[:, j]), color=COLORS["R"], lw=2, label="R — Recovered")

    ax.set_title(f"Patch {j}: {name}  (N = {pop_j:,})", fontsize=10)
    ax.set_ylabel("Population (%)")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9, loc="right")
    ax.grid(True, alpha=0.3, linestyle=":")

for ax in axes[1]:
    ax.set_xlabel("Day")

plt.tight_layout()
plt.savefig("eastland_seir.png", dpi=150, bbox_inches="tight")
print("\nPlot saved -> eastland_seir.png")
plt.show()
```

---

**Key design decisions explained:**

**`beta` calculation:** `beta = R0 / mean_infectious_period = 5.0 / 10 = 0.5 d⁻¹`. This is the standard relationship in an SIR/SEIR model where `R0 = beta × D_infectious`.

**Scenario GeoDataFrame:** Built with `laser.core.utils.grid(M=1, N=4)` — a single-row four-column grid. The `population_fn(row, col)` lambda maps each column index to the corresponding patch population. Fictional Eastern-European coordinates are used so the framework can compute a realistic distance/gravity network between patches.

**Initial conditions (`scenario["S/E/I/R"]`):** `S` is derived last (`S = pop − E − I − R`) so that integer rounding never violates population conservation.

**Component constructor arguments:**
- `SEIR.Exposed(model, expdist, infdist)` — needs *both* distributions: `expdist` to assign incubation timers to newly-seeded exposed agents, and `infdist` because the Exposed component itself fires the E→I transition and must immediately assign the infectious timer.
- `SEIR.Infectious(model, infdist)` — assigns infectious timers to the agents initialised as infectious at `t=0`.
- `SEIR.Transmission(model, expdist)` — uses `expdist` to assign incubation timers whenever a susceptible agent is exposed (S→E).

**Component order `[s, e, i, r, tx]`:** This is the canonical order established in the laser-generic test suite. Within each tick the flow propagation is: existing E→I and I→R transitions fire first, then new S→E exposures are computed from the current I count.
