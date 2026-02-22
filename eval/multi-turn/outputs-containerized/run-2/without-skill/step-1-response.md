Now I have all the details needed. Here is the complete implementation:

```python
"""
Basic 4-patch SEIR model for a respiratory disease using the LASER framework
(laser-generic package).

Disease parameters:
  - R0 = 5
  - Latent period = 4 days (mean)
  - Infectious period = 10 days (mean)
  - beta = R0 / infectious_period = 0.5 /day

Patch populations:  100k, 200k, 150k, 80k
Initial conditions: 90% S, 0% E, 1% I, 9% R (uniform across all patches)
Duration: 365 days (1 year)

No spatial coupling (no migration), no demographics, no seasonal forcing.
"""

import numpy as np
import matplotlib.pyplot as plt
import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.utils import grid, initialize_population
from laser.generic import SEIR, Model


# ─────────────────────────────────────────────────────────────
# 1. Disease and simulation parameters
# ─────────────────────────────────────────────────────────────
R0                = 5.0
LATENT_PERIOD     = 4    # days (mean exposed duration)
INFECTIOUS_PERIOD = 10   # days (mean infectious duration)
NTICKS            = 365  # 1 year

# Transmission rate: R0 = beta * infectious_period  →  beta = R0 / infectious_period
BETA = R0 / INFECTIOUS_PERIOD   # 0.5 /day

# ─────────────────────────────────────────────────────────────
# 2. Patch populations
# ─────────────────────────────────────────────────────────────
POPULATIONS = [100_000, 200_000, 150_000, 80_000]   # one entry per patch
N_PATCHES   = len(POPULATIONS)

def _pop_fn(pops):
    """Factory that returns a (row, col) → population function for grid()."""
    def fn(row, col):
        return pops[col]          # 1-row grid: col index = patch index
    return fn

# Create a 1×4 spatial grid; geometry is Polygon cells in geographic degrees.
# The gravity network is computed automatically from patch centroids and
# populations but is never activated (no migration component is added).
scenario = grid(M=1, N=4, population_fn=_pop_fn(POPULATIONS))

# Set initial compartment fractions uniformly across all patches.
# S is assigned last (as the remainder) to preserve integer totals exactly.
#   S=90%, E=0%, I=1%, R=9%
scenario = initialize_population(
    scenario,
    [0.90, 0.00, 0.01, 0.09],   # [S_frac, E_frac, I_frac, R_frac]
    states=["S", "E", "I", "R"],
)

# ─────────────────────────────────────────────────────────────
# 3. Model parameters
# ─────────────────────────────────────────────────────────────
params = PropertySet({
    "nticks":    NTICKS,
    "beta":      BETA,
    "prng_seed": 42,
})

# ─────────────────────────────────────────────────────────────
# 4. Initialise the Model
# ─────────────────────────────────────────────────────────────
model = Model(scenario, params, name="seir_4patch_respiratory")

# ─────────────────────────────────────────────────────────────
# 5. Duration distributions
# ─────────────────────────────────────────────────────────────
# Gamma(shape=k, scale=θ)  →  mean = k·θ
#   shape=4, scale=1  →  mean = 4 days  (latent period)
exp_dur_dist = dists.gamma(shape=float(LATENT_PERIOD), scale=1.0)

# Normal(loc=μ, scale=σ)   →  mean = μ
#   loc=10, scale=2  →  mean = 10 days  (infectious period)
inf_dur_dist = dists.normal(loc=float(INFECTIOUS_PERIOD), scale=2.0)

# ─────────────────────────────────────────────────────────────
# 6. Assemble SEIR components
# ─────────────────────────────────────────────────────────────
# Susceptible  – initialises agent nodeid/state and node-level S counts
# Exposed      – manages etimer (E→I transition) using exp_dur_dist;
#                also seeds itimer for when agents enter I
# Infectious   – manages itimer (I→R transition) using inf_dur_dist
# Recovered    – passive tracker; agents remain R indefinitely (no waning)
# Transmission – force-of-infection per patch: S→E at rate β·I/N
#                (exp_dur_dist used to assign etimer to newly exposed agents)

s_comp  = SEIR.Susceptible(model)
e_comp  = SEIR.Exposed(model, exp_dur_dist, inf_dur_dist)
i_comp  = SEIR.Infectious(model, inf_dur_dist)
r_comp  = SEIR.Recovered(model)
tx_comp = SEIR.Transmission(model, exp_dur_dist)

# No migration component → each patch's force of infection is self-contained;
# patches evolve as independent SEIR epidemics.
model.components = [s_comp, e_comp, i_comp, r_comp, tx_comp]

# ─────────────────────────────────────────────────────────────
# 7. Pre-run summary
# ─────────────────────────────────────────────────────────────
PATCH_NAMES = [f"Patch {i} (N={p//1000}k)" for i, p in enumerate(POPULATIONS)]

print("=" * 62)
print("  SEIR 4-Patch Respiratory Disease Model (LASER framework)")
print("=" * 62)
print(f"  R0 = {R0}  |  beta = {BETA:.3f} /day")
print(f"  Latent period    = {LATENT_PERIOD} days  "
      f"[Gamma(shape={LATENT_PERIOD}, scale=1)]")
print(f"  Infectious period = {INFECTIOUS_PERIOD} days  "
      f"[Normal(loc={INFECTIOUS_PERIOD}, scale=2)]")
print(f"  Duration          = {NTICKS} days (1 year)")
print()
print("  Initial patch states:")
hdr0 = f"  {'Patch':<22} {'Pop':>10} {'S':>10} {'E':>8} {'I':>8} {'R':>10}"
print(hdr0)
print("  " + "-" * (len(hdr0) - 2))
for idx in range(N_PATCHES):
    row = scenario.iloc[idx]
    print(f"  {PATCH_NAMES[idx]:<22} {int(row.population):>10,} "
          f"{int(row.S):>10,} {int(row.E):>8,} "
          f"{int(row.I):>8,} {int(row.R):>10,}")
print()

# ─────────────────────────────────────────────────────────────
# 8. Run the simulation
# ─────────────────────────────────────────────────────────────
model.run("SEIR 4-Patch Respiratory Disease")

# ─────────────────────────────────────────────────────────────
# 9. Extract results
#    Arrays are shaped (nticks+1, n_patches) = (366, 4):
#      axis-0 = time  (day 0 … day 365)
#      axis-1 = patch (patch 0 … patch 3)
# ─────────────────────────────────────────────────────────────
S_arr = model.nodes.S
E_arr = model.nodes.E
I_arr = model.nodes.I
R_arr = model.nodes.R

S_total = S_arr.sum(axis=1)
E_total = E_arr.sum(axis=1)
I_total = I_arr.sum(axis=1)
R_total = R_arr.sum(axis=1)

days = np.arange(NTICKS + 1)   # 0, 1, …, 365

# ─────────────────────────────────────────────────────────────
# 10. Post-run summary
# ─────────────────────────────────────────────────────────────
init_S = [int(scenario["S"].iloc[i]) for i in range(N_PATCHES)]

print("\n=== Final state (Day 365) ===")
hdr = (f"{'Patch':<22} {'Pop':>10} {'Final S':>10} {'Final E':>8} "
       f"{'Final I':>8} {'Final R':>10} {'Attack%':>9}")
sep = "-" * len(hdr)
print(hdr)
print(sep)
for idx in range(N_PATCHES):
    pop    = POPULATIONS[idx]
    fs     = int(S_arr[-1, idx])
    fe     = int(E_arr[-1, idx])
    fi     = int(I_arr[-1, idx])
    fr     = int(R_arr[-1, idx])
    attack = (init_S[idx] - fs) / init_S[idx] * 100
    print(f"{PATCH_NAMES[idx]:<22} {pop:>10,} {fs:>10,} {fe:>8,} "
          f"{fi:>8,} {fr:>10,} {attack:>8.1f}%")

total_pop    = sum(POPULATIONS)
total_final  = int(S_total[-1] + E_total[-1] + I_total[-1] + R_total[-1])
total_attack = (sum(init_S) - int(S_total[-1])) / sum(init_S) * 100
print(sep)
print(f"{'All patches':<22} {total_pop:>10,} {int(S_total[-1]):>10,} "
      f"{int(E_total[-1]):>8,} {int(I_total[-1]):>8,} "
      f"{int(R_total[-1]):>10,} {total_attack:>8.1f}%")

peak_day   = int(np.argmax(I_total))
peak_count = int(I_total[peak_day])
print(f"\nPopulation conservation: initial={total_pop:,}, final={total_final:,}, "
      f"diff={abs(total_pop - total_final)}")
print(f"Peak infectious (all patches combined): {peak_count:,} on day {peak_day}")

# ─────────────────────────────────────────────────────────────
# 11. Visualisation: 4 per-patch panels + 1 aggregate panel
# ─────────────────────────────────────────────────────────────
COLORS = {
    "S": "#1f77b4",   # blue
    "E": "#9467bd",   # purple
    "I": "#d62728",   # red
    "R": "#2ca02c",   # green
}

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle(
    f"SEIR 4-Patch Respiratory Disease Model  "
    f"(R0={R0}, latent={LATENT_PERIOD}d, infectious={INFECTIOUS_PERIOD}d, β={BETA:.2f})",
    fontsize=13, fontweight="bold",
)

# Per-patch subplots: layout is [[patch0, patch1, empty], [patch2, patch3, agg]]
LAYOUT = [(0, 0), (0, 1), (1, 0), (1, 1)]
for idx, (r, c) in enumerate(LAYOUT):
    ax  = axes[r, c]
    pop = POPULATIONS[idx]
    ax.plot(days, S_arr[:, idx] / 1e3, label="S", color=COLORS["S"], lw=2)
    ax.plot(days, E_arr[:, idx] / 1e3, label="E", color=COLORS["E"], lw=2)
    ax.plot(days, I_arr[:, idx] / 1e3, label="I", color=COLORS["I"], lw=2)
    ax.plot(days, R_arr[:, idx] / 1e3, label="R", color=COLORS["R"], lw=2)
    ax.set_title(PATCH_NAMES[idx], fontsize=11)
    ax.set_xlabel("Day")
    ax.set_ylabel("Population (thousands)")
    ax.legend(loc="center right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, NTICKS)

# Aggregate panel (bottom-right)
ax_agg = axes[1, 2]
ax_agg.plot(days, S_total / 1e3, label="S", color=COLORS["S"], lw=2)
ax_agg.plot(days, E_total / 1e3, label="E", color=COLORS["E"], lw=2)
ax_agg.plot(days, I_total / 1e3, label="I", color=COLORS["I"], lw=2)
ax_agg.plot(days, R_total / 1e3, label="R", color=COLORS["R"], lw=2)
ax_agg.axvline(peak_day, color="gray", ls="--", lw=1,
               label=f"I peak day {peak_day}")
ax_agg.set_title("All Patches Combined", fontsize=11)
ax_agg.set_xlabel("Day")
ax_agg.set_ylabel("Population (thousands)")
ax_agg.legend(loc="center right", fontsize=9)
ax_agg.grid(True, alpha=0.3)
ax_agg.set_xlim(0, NTICKS)

# Hide the unused top-right cell
axes[0, 2].axis("off")

plt.tight_layout()
plt.savefig("seir_4patch_results.png", dpi=150, bbox_inches="tight")
print("\nFigure saved: seir_4patch_results.png")
plt.show()
```

---

**Key design decisions:**

| Choice | Rationale |
|---|---|
| `beta = R0 / infectious_period = 0.5` | Standard relation for discrete-time per-contact rate; matches LASER test convention |
| `exp_dur_dist = dists.gamma(shape=4, scale=1)` | Mean = 4 days (latent); Gamma is epidemiologically appropriate for incubation |
| `inf_dur_dist = dists.normal(loc=10, scale=2)` | Mean = 10 days (infectious); matches LASER test style |
| `initialize_population([0.90, 0.00, 0.01, 0.09])` | S gets remainder, ensuring S+E+I+R = pop exactly for each patch |
| No migration component in `model.components` | The gravity network is computed inside `Model.__init__` but never consulted; patches evolve independently |
| `grid(M=1, N=4, ...)` | Creates GeoDataFrame with `nodeid`, `population`, Polygon `geometry` columns — the format `Model` expects |
