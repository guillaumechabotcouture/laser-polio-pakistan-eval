Now I have all the information I need. `laser_core.LaserFrame` and `laser_core.PropertySet` work correctly. The `laser.generic` high-level components fail due to a `laser.core` namespace mismatch in this environment, so I'll implement the SEIR dynamics following the exact same LASER component pattern (`step(tick)` methods, `_initialize_flows`, time-first state arrays).

```python
"""
4-Patch SEIR Model using the LASER framework (laser-core).

Disease parameters
------------------
  R0                = 5.0
  Latent period     = 4 days   → sigma = 1/4  = 0.25 /day  (E → I rate)
  Infectious period = 10 days  → gamma = 1/10 = 0.10 /day  (I → R rate)
  beta              = R0 * gamma = 0.50 /day

Patches (no spatial coupling)
------------------------------
  Patch 0: 100 000 people
  Patch 1: 200 000 people
  Patch 2: 150 000 people
  Patch 3:  80 000 people

Initial conditions (applied to every patch)
--------------------------------------------
  S₀ = 90%,  E₀ = 0%,  I₀ = 1%,  R₀ = 9%

Duration: 365 days (1 year)

Framework notes
---------------
  * laser_core.PropertySet  — parameter container with attribute-style access
  * laser_core.LaserFrame   — columnar store; node state arrays are vector
    properties with shape (nticks+1, n_patches), i.e. time-first indexing,
    matching the laser.generic convention.
  * Components follow the laser.generic pattern: each has a step(tick) method
    that reads state[tick] and accumulates changes into state[tick+1].
  * _initialize_flows copies state[tick] → state[tick+1] at the start of every
    tick so each component only has to add or subtract its own flow.
"""

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point

# laser_core ships the core primitives used by laser-generic
from laser_core import LaserFrame, PropertySet


# =========================================================================
# 1.  PARAMETERS
# =========================================================================

R0                = 5.0
latent_period     = 4.0    # days  (mean time E → I)
infectious_period = 10.0   # days  (mean time I → R)

gamma = 1.0 / infectious_period   # 0.10 /day  recovery rate
sigma = 1.0 / latent_period       # 0.25 /day  progression rate
beta  = R0 * gamma                # 0.50 /day  transmission rate

params = PropertySet({
    "nticks":             365,          # 1 year
    "beta":               beta,
    "sigma":              sigma,
    "gamma":              gamma,
    "R0":                 R0,
    "latent_period":      latent_period,
    "infectious_period":  infectious_period,
    "prng_seed":          20260101,
})


# =========================================================================
# 2.  SCENARIO  —  4-patch GeoDataFrame (laser Model convention)
# =========================================================================

patch_populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int64)
n_patches         = len(patch_populations)

# Initial fractions: 90% S, 0% E, 1% I, 9% R
S0    = np.floor(patch_populations * 0.90).astype(np.int64)
E0    = np.zeros(n_patches, dtype=np.int64)
I0    = np.floor(patch_populations * 0.01).astype(np.int64)
R0arr = patch_populations - S0 - E0 - I0        # absorbs rounding → ~9%

# Arbitrary point geometries (required by laser Model convention)
coords = [(-87.6, 41.8), (-88.0, 40.7), (-86.2, 39.8), (-89.4, 43.1)]

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(n_patches),
        "population": patch_populations,
        "S":          S0,
        "E":          E0,
        "I":          I0,
        "R":          R0arr,
    },
    geometry=[Point(lon, lat) for lon, lat in coords],
    crs="EPSG:4326",
)

print("=" * 60)
print("4-Patch SEIR Model — LASER Framework (laser-core)")
print("=" * 60)
print(f"\nParameters: R0={params.R0}  beta={params.beta:.4f}"
      f"  sigma={params.sigma:.4f}  gamma={params.gamma:.4f}")
print(f"Duration  : {params.nticks} days\n")
print("Initial scenario:")
print(scenario[["nodeid", "population", "S", "E", "I", "R"]].to_string(index=False))
print()


# =========================================================================
# 3.  NODE STATE STORE  —  LaserFrame with time-series vector properties
#
#     nodes.<state>[tick, patch]  →  shape (nticks+1, n_patches)
#     This matches the laser.generic convention for compartmental models.
# =========================================================================

nticks = params.nticks

nodes = LaserFrame(n_patches)

# Compartment time-series
nodes.add_vector_property("S", nticks + 1, dtype=np.float64)
nodes.add_vector_property("E", nticks + 1, dtype=np.float64)
nodes.add_vector_property("I", nticks + 1, dtype=np.float64)
nodes.add_vector_property("R", nticks + 1, dtype=np.float64)

# Incidence tracking (flows between compartments)
nodes.add_vector_property("newly_exposed",    nticks + 1, dtype=np.float64)
nodes.add_vector_property("newly_infectious", nticks + 1, dtype=np.float64)
nodes.add_vector_property("newly_recovered",  nticks + 1, dtype=np.float64)

# Seed tick 0 from scenario
nodes.S[0] = S0.astype(np.float64)
nodes.E[0] = E0.astype(np.float64)
nodes.I[0] = I0.astype(np.float64)
nodes.R[0] = R0arr.astype(np.float64)


# =========================================================================
# 4.  SEIR COMPONENTS  (laser.generic component pattern)
#
#     Each component owns one transition.  It reads state[tick] (unchanged
#     original values) and accumulates its delta into state[tick+1] (which
#     _initialize_flows already pre-filled with state[tick]).
# =========================================================================

class TransmissionSE:
    """
    S → E  :  new exposures driven by force of infection.

    FOI = beta * I[t] / N   (no spatial coupling between patches).
    new_E = FOI * S[t]

    Updates S[t+1] -= new_E,  E[t+1] += new_E.
    """

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes  = nodes
        self.params = params

    def step(self, tick: int) -> None:
        n   = self.nodes
        N   = n.S[tick] + n.E[tick] + n.I[tick] + n.R[tick]
        foi = self.params.beta * n.I[tick] / N   # force of infection (per patch)
        new_E = foi * n.S[tick]

        n.newly_exposed[tick + 1]  = new_E
        n.S[tick + 1]             -= new_E
        n.E[tick + 1]             += new_E


class ExposedProgression:
    """
    E → I  :  end of latent period, agents become infectious.

    new_I = sigma * E[t]

    Updates E[t+1] -= new_I,  I[t+1] += new_I.
    """

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes  = nodes
        self.params = params

    def step(self, tick: int) -> None:
        n     = self.nodes
        new_I = self.params.sigma * n.E[tick]

        n.newly_infectious[tick + 1] = new_I
        n.E[tick + 1]               -= new_I
        n.I[tick + 1]               += new_I


class InfectiousRecovery:
    """
    I → R  :  recovery after infectious period.

    new_R = gamma * I[t]

    Updates I[t+1] -= new_R,  R[t+1] += new_R.
    """

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes  = nodes
        self.params = params

    def step(self, tick: int) -> None:
        n     = self.nodes
        new_R = self.params.gamma * n.I[tick]

        n.newly_recovered[tick + 1] = new_R
        n.I[tick + 1]              -= new_R
        n.R[tick + 1]              += new_R


# =========================================================================
# 5.  MODEL ORCHESTRATOR
#
#     Mirrors laser.generic.model.Model.run():
#       _initialize_flows  — copy state[t] → state[t+1]  (LASER convention)
#       component.step(t)  — each component adds/subtracts its flow
# =========================================================================

class SEIRModel:
    """
    Patch-level SEIR orchestrator following the laser.generic.Model design.

    Attributes
    ----------
    scenario   : GeoDataFrame  — patch metadata (as in laser Model)
    nodes      : LaserFrame    — time-series state store
    params     : PropertySet   — simulation parameters
    components : list          — ordered list of SEIR components
    """

    def __init__(
        self,
        scenario:   gpd.GeoDataFrame,
        nodes:      LaserFrame,
        params:     PropertySet,
        components: list,
    ) -> None:
        self.scenario   = scenario
        self.nodes      = nodes
        self.params     = params
        self.components = components

    def _initialize_flows(self, tick: int) -> None:
        """
        Forward-fill all state arrays from tick to tick+1.

        This is the laser.generic convention: every component then only needs
        to apply its own delta to state[tick+1].
        """
        n = self.nodes
        n.S[tick + 1] = n.S[tick]
        n.E[tick + 1] = n.E[tick]
        n.I[tick + 1] = n.I[tick]
        n.R[tick + 1] = n.R[tick]

    def run(self) -> None:
        for tick in range(self.params.nticks):
            self._initialize_flows(tick)
            for component in self.components:
                component.step(tick)
        print(f"Simulation complete: {self.params.nticks} days, "
              f"{self.nodes.count} patches.\n")


# =========================================================================
# 6.  BUILD AND RUN
# =========================================================================

model = SEIRModel(
    scenario=scenario,
    nodes=nodes,
    params=params,
    components=[
        TransmissionSE(nodes, params),
        ExposedProgression(nodes, params),
        InfectiousRecovery(nodes, params),
    ],
)

model.run()


# =========================================================================
# 7.  RESULTS SUMMARY
# =========================================================================

print("=== Final state  (day 365) ===")
header = f"{'Patch':>5}  {'N':>8}  {'S':>10}  {'E':>6}  {'I':>6}  {'R':>10}  {'Attack%':>8}"
print(header)
print("-" * len(header))

for p in range(n_patches):
    N   = int(patch_populations[p])
    Sf  = nodes.S[-1, p]
    Ef  = nodes.E[-1, p]
    If  = nodes.I[-1, p]
    Rf  = nodes.R[-1, p]
    # Attack rate = new recoveries as fraction of initially susceptible
    new_R       = Rf - float(R0arr[p])
    attack_rate = new_R / float(S0[p]) * 100.0
    print(f"{p:>5}  {N:>8,}  {Sf:>10.0f}  {Ef:>6.1f}  {If:>6.1f}  "
          f"{Rf:>10.0f}  {attack_rate:>7.1f}%")

print("\n=== Peak infectious ===")
for p in range(n_patches):
    peak_day = int(nodes.I[:, p].argmax())
    peak_I   = nodes.I[peak_day, p]
    print(f"  Patch {p} (N={patch_populations[p]:,}): "
          f"peak I = {peak_I:.0f}  ({peak_I / patch_populations[p] * 100:.1f}%)  "
          f"on day {peak_day}")


# =========================================================================
# 8.  VISUALISATION
# =========================================================================

days   = np.arange(nticks + 1)
colors = {"S": "steelblue", "E": "darkorange", "I": "firebrick", "R": "forestgreen"}

fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
axes = axes.flatten()

for p, ax in enumerate(axes):
    N = patch_populations[p]
    ax.plot(days, nodes.S[:, p] / N * 100, lw=2, label="S", color=colors["S"])
    ax.plot(days, nodes.E[:, p] / N * 100, lw=2, label="E", color=colors["E"])
    ax.plot(days, nodes.I[:, p] / N * 100, lw=2, label="I", color=colors["I"])
    ax.plot(days, nodes.R[:, p] / N * 100, lw=2, label="R", color=colors["R"])
    ax.set_title(f"Patch {p}  (N = {N:,})", fontsize=11)
    ax.set_ylabel("Population (%)")
    ax.set_xlabel("Day")
    ax.set_xlim(0, nticks)
    ax.set_ylim(0, 100)
    ax.legend(loc="center right", fontsize=9)
    ax.grid(True, alpha=0.3)

fig.suptitle(
    f"4-Patch SEIR  |  R₀ = {R0}  ·  latent = {latent_period} d  "
    f"·  infectious = {infectious_period} d  |  No spatial coupling",
    fontsize=12,
)
plt.tight_layout()
plt.savefig("seir_4patch.png", dpi=150, bbox_inches="tight")
print("\nPlot saved → seir_4patch.png")
plt.show()
```

**How it works:**

| Layer | What it uses |
|---|---|
| Parameters | `laser_core.PropertySet` — dict with attribute access |
| Node state | `laser_core.LaserFrame` — vector properties shaped `(nticks+1, n_patches)` |
| Scenario | `geopandas.GeoDataFrame` — standard laser Model input convention |
| Components | Custom classes with `step(tick)` matching the `laser.generic` pattern |
| Orchestrator | `SEIRModel._initialize_flows` + `component.step(tick)` loop, mirroring `laser.generic.model.Model.run` |

**SEIR dynamics (per tick, per patch, no coupling):**

```
FOI   = β · I[t] / N
ΔS→E  = FOI · S[t]          (TransmissionSE)
ΔE→I  = σ · E[t]            (ExposedProgression)
ΔI→R  = γ · I[t]            (InfectiousRecovery)
```

With β = 0.50, σ = 0.25, γ = 0.10 → R₀ = β/γ = 5. Each patch runs an independent epidemic; expect ~≥95% attack rate over 1 year given R₀ = 5.
