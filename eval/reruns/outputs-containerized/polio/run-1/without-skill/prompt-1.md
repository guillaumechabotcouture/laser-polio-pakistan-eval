Here is the complete simulation code, also printed below:

---

```python
"""
Spatial SEIR Model for Polio Transmission Dynamics across Pakistan Districts
============================================================================
Gates Foundation / Institute for Disease Modeling – LASER Framework

Configuration
─────────────
  Patches          : 10 districts, each 100,000 population
  R0               : 6  (poliovirus)
  Latent period    : 3  days mean  [Gamma(shape=3, scale=1)]
  Infectious period: 28 days mean  [Normal(loc=28, scale=5)]
  beta             : R0 / infectious_period ≈ 0.2143 per day
  Initial state    : 95 % recovered/immune, 5 infectious per patch
  Simulation       : 10 years (3 650 daily time-steps)
"""

import numpy as np
import matplotlib.pyplot as plt

import laser.core.distributions as dists
from laser.core import PropertySet
from laser.core.utils import grid
from laser.generic import SEIR, Model

# ── Epidemiological parameters ────────────────────────────────────────────────
N_PATCHES     = 10
POP_PER_PATCH = 100_000
R0            = 6.0
LATENT_MEAN   = 3.0   # days (mean latent / pre-infectious period)
INFECT_MEAN   = 28.0  # days (mean infectious period)
BETA          = R0 / INFECT_MEAN  # ≈ 0.2143 per day

# ── Simulation parameters ─────────────────────────────────────────────────────
YEARS  = 10
NTICKS = YEARS * 365  # 3 650 daily steps

# ── Initial conditions per patch ──────────────────────────────────────────────
INIT_INFECTIOUS  = 5
INIT_RECOVERED   = int(0.95 * POP_PER_PATCH)  # 95 000 immune
INIT_EXPOSED     = 0
INIT_SUSCEPTIBLE = POP_PER_PATCH - INIT_RECOVERED - INIT_INFECTIOUS - INIT_EXPOSED
# = 100 000 - 95 000 - 5 = 4 995

# ── Effective reproduction number at start ─────────────────────────────────────
S_FRAC    = INIT_SUSCEPTIBLE / POP_PER_PATCH   # fraction susceptible
REFF_INIT = R0 * S_FRAC                         # initial effective R

print("=" * 65)
print("Polio SEIR Spatial Model – 10 Pakistan Districts")
print("=" * 65)
print(f"  Population per patch   : {POP_PER_PATCH:,}")
print(f"  R0                     : {R0}")
print(f"  beta  (= R0 / D_I)     : {BETA:.4f} per day")
print(f"  Latent period mean     : {LATENT_MEAN} days  [Gamma(shape=3, scale=1)]")
print(f"  Infectious period mean : {INFECT_MEAN} days  [Normal(loc=28, scale=5)]")
print(f"  Initial S per patch    : {INIT_SUSCEPTIBLE:,}  ({100 * S_FRAC:.2f} %)")
print(f"  Initial E per patch    : {INIT_EXPOSED}")
print(f"  Initial I per patch    : {INIT_INFECTIOUS}")
print(f"  Initial R per patch    : {INIT_RECOVERED:,}  ({100 * INIT_RECOVERED / POP_PER_PATCH:.0f} %)")
print(f"  Reff (initial)         : {REFF_INIT:.3f}  "
      f"({'< 1 → no sustained outbreak' if REFF_INIT < 1 else '>= 1 → outbreak possible'})")
print(f"  Herd immunity threshold: {100 * (1 - 1 / R0):.1f} % immune required")
print(f"  Simulation             : {YEARS} years  ({NTICKS:,} days)")
print()

# ── Scenario: 10 patches in a 1 × 10 grid, NW Pakistan (KPK region) ──────────
# Origin: 33 °N, 71 °E; node_size_degs ≈ 0.08983 ° ≈ 10 km per cell
scenario = grid(
    M=1,
    N=N_PATCHES,
    node_size_degs=0.08983,
    population_fn=lambda row, col: POP_PER_PATCH,
    origin_x=71.0,   # longitude (°E)
    origin_y=33.0,   # latitude  (°N)
)

# Assign initial compartment counts for every patch
scenario["S"] = INIT_SUSCEPTIBLE
scenario["E"] = INIT_EXPOSED
scenario["I"] = INIT_INFECTIOUS
scenario["R"] = INIT_RECOVERED

print("Scenario – one row per district:")
print(scenario[["nodeid", "population", "S", "E", "I", "R"]].to_string(index=False))
print()

# ── Model parameters ──────────────────────────────────────────────────────────
params = PropertySet({
    "nticks"    : NTICKS,
    "beta"      : BETA,
    "prng_seed" : 20260101,
    # Gravity model for spatial coupling (standard defaults)
    "gravity_k" : 500,
    "gravity_a" : 1,
    "gravity_b" : 1,
    "gravity_c" : 2,
})

# ── Initialise LASER Model ────────────────────────────────────────────────────
print("Initialising LASER Model …")
model = Model(scenario, params)

# ── Duration distributions (callables: (tick, node_id) → float) ───────────────
# Latent (exposed) period – Gamma(shape=3, scale=1)  → mean = 3 days
expdist = dists.gamma(shape=LATENT_MEAN, scale=1.0)

# Infectious period – Normal(loc=28, scale=5)  → mean = 28 days
infdist = dists.normal(loc=INFECT_MEAN, scale=5.0)

# ── SEIR components ───────────────────────────────────────────────────────────
#   Susceptible  – records S counts per node/tick
#   Exposed      – E→I transition via etimer (expdist) then itimer (infdist)
#   Infectious   – I→R transition via itimer (infdist)
#   Recovered    – records R counts per node/tick
#   Transmission – force of infection; drives S→E, assigns etimer (expdist)
s  = SEIR.Susceptible(model)
e  = SEIR.Exposed(model, expdist, infdist)
i  = SEIR.Infectious(model, infdist)
r  = SEIR.Recovered(model)
tx = SEIR.Transmission(model, expdist)

model.components = [s, e, i, r, tx]

# ── Run ───────────────────────────────────────────────────────────────────────
print("Running 10-year simulation …")
model.run("Pakistan Polio SEIR – 10 districts, 10 years")

# ── Extract results ───────────────────────────────────────────────────────────
# Shape: (nticks + 1, n_patches);  axis 0 = tick,  axis 1 = patch
S_arr = model.nodes.S
E_arr = model.nodes.E
I_arr = model.nodes.I
R_arr = model.nodes.R

# Sum across all patches
S_tot = S_arr.sum(axis=1)
E_tot = E_arr.sum(axis=1)
I_tot = I_arr.sum(axis=1)
R_tot = R_arr.sum(axis=1)
N_tot = S_tot + E_tot + I_tot + R_tot

ticks = np.arange(len(S_tot))
years = ticks / 365.0

# ── Summary statistics ────────────────────────────────────────────────────────
print()
print("=" * 65)
print("Simulation Results Summary")
print("=" * 65)
print(f"  Initial total population   : {N_tot[0]:>10,.0f}")
print(f"  Final total population     : {N_tot[-1]:>10,.0f}")
print(f"  Population drift (ΔN)      : {N_tot[-1] - N_tot[0]:>+10.1f}  (conservation check)")
print()
print(f"  Initial S : {S_tot[0]:>8,.0f}  ({100 * S_tot[0] / N_tot[0]:.2f} %)")
print(f"  Final   S : {S_tot[-1]:>8,.0f}  ({100 * S_tot[-1] / N_tot[-1]:.2f} %)")
print()
peak_I_val = int(I_tot.max())
peak_I_day = int(I_tot.argmax())
print(f"  Peak infectious (all patches) : {peak_I_val:,}  at day {peak_I_day} "
      f"(year {peak_I_day / 365:.1f})")
print(f"  Final infectious              : {int(I_tot[-1]):,}")
print()
delta_R_total = int(R_tot[-1] - R_tot[0])
print(f"  Cumulative new infections (ΔR): {delta_R_total:,}")
print(f"  Final recovered / immune      : {int(R_tot[-1]):,}  "
      f"({100 * R_tot[-1] / N_tot[-1]:.2f} %)")
print()
print(f"  Effective R at start: Reff = R0 × (S/N) = {R0} × {S_FRAC:.4f} = {REFF_INIT:.3f}")
print(f"  {'No sustained outbreak expected (Reff < 1)' if REFF_INIT < 1 else 'Outbreak possible (Reff >= 1)'}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    "Polio SEIR Transmission Model – 10 Pakistan Districts\n"
    f"R₀ = {R0},  Latent = {LATENT_MEAN} d,  Infectious = {INFECT_MEAN} d,  "
    f"Initial immunity = 95 %",
    fontsize=13, fontweight="bold",
)

# Panel 1 – All-patch SEIR compartments
ax = axes[0, 0]
ax.plot(years, S_tot / 1e3, label="Susceptible (S)", color="steelblue",   lw=1.8)
ax.plot(years, E_tot / 1e3, label="Exposed (E)",     color="darkorange",  lw=1.8)
ax.plot(years, I_tot / 1e3, label="Infectious (I)",  color="crimson",     lw=1.8)
ax.plot(years, R_tot / 1e3, label="Immune / Rec. (R)", color="forestgreen", lw=1.8)
ax.set_xlabel("Time (years)")
ax.set_ylabel("Individuals (thousands)")
ax.set_title("SEIR Compartments – All 10 Districts")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2 – Infectious cases per district
ax = axes[0, 1]
cmap_colors = plt.cm.tab10(np.linspace(0, 1, N_PATCHES))
for p in range(N_PATCHES):
    ax.plot(years, I_arr[:, p], color=cmap_colors[p], lw=1.0, alpha=0.85,
            label=f"District {p + 1}")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Infectious individuals")
ax.set_title("Infectious Cases by District")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.3)

# Panel 3 – Susceptible fraction over time
ax = axes[1, 0]
pct_S = 100.0 * S_tot / N_tot
ax.plot(years, pct_S, color="steelblue", lw=1.8, label="S fraction")
ax.axhline(100.0 * S_FRAC, color="steelblue", ls="--", alpha=0.5,
           label=f"Initial S = {100 * S_FRAC:.1f} %")
herd_pct = 100.0 * (1.0 - 1.0 / R0)
ax.axhline(herd_pct, color="dimgray", ls=":", lw=1.2,
           label=f"Herd immunity threshold = {herd_pct:.1f} %")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Susceptible (%)")
ax.set_title("Susceptible Population Fraction")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=0)

# Panel 4 – Cumulative new infections
ax = axes[1, 1]
delta_R = R_tot - R_tot[0]
ax.plot(years, delta_R, color="forestgreen", lw=1.8)
ax.set_xlabel("Time (years)")
ax.set_ylabel("Cumulative new infections")
ax.set_title("Cumulative New Infections (ΔR)")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = "polio_seir_simulation.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nPlot saved → {out_path}")
```

---

Here is what the code does, section by section:

**Parameter derivation**

| Quantity | Value | Derivation |
|---|---|---|
| `BETA` | 0.2143 / day | `R0 / INFECT_MEAN = 6 / 28` |
| `INIT_SUSCEPTIBLE` | 4 995 / patch | `100 000 − 95 000 − 5` |
| `Reff` (initial) | ≈ 0.30 | `R0 × (S/N) = 6 × 0.04995` |

With 95 % initial immunity the effective reproduction number is well below 1, so no sustained outbreak occurs — this is the herd-immunity scenario that vaccination campaigns aim to achieve.

**LASER components wired up**

```
Transmission (S→E, assigns etimer ~ Gamma(3,1))
Exposed      (E→I, decrements etimer; assigns itimer ~ Normal(28,5))
Infectious   (I→R, decrements itimer)
Susceptible / Recovered  (bookkeeping arrays)
```

**Spatial coupling** — `Model` automatically builds a gravity network from the 10-node GeoDataFrame and `gravity_k/a/b/c` parameters, so cross-district transmission is handled without extra code.

**Results** — after `model.run()` the per-district time-series live in `model.nodes.S/E/I/R`, shaped `(nticks+1, 10)`, ready for further analysis or visualisation.
