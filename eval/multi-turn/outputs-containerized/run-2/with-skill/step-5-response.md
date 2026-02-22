Here is the complete script with vital dynamics added:

```python
"""
4-Patch Gravity SEIR Model — Seasonal Forcing + Vital Dynamics
==============================================================
Disease: respiratory (R0=5, latent=4d, infectious=10d)
Patches: pop 100k / 200k / 150k / 80k, arranged in a line 75 km apart
Gravity: k=0.01, a=1, b=1, c=1.5; row-normalised, max 15% export per patch
Seasonal: winter peak days 0–90 → 1.3x; summer trough days 150–240 → 0.7x
          linear transitions between; implemented via LASER ValuesMap
Vital:   CBR=30/1000/yr, CDR=10/1000/yr → ~2% net annual growth
         BirthsByCBR + MortalityByCDR; capacity pre-allocated via calc_capacity
Run:     3650 days (10 years) — gravity spatial coupling + vital dynamics
Package: laser-generic v1.0.0
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.utils import calc_capacity
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Simulation settings ────────────────────────────────────────────────────
NTICKS = 10 * 365   # 3650 days (10 years)
POPS   = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
NNODES = len(POPS)
SEED   = 42

# ── Disease parameters ─────────────────────────────────────────────────────
R0       = 5.0
LATENT_D = 4
INFECT_D = 10
BETA     = R0 / INFECT_D  # = 0.5 per day

# ── Vital dynamics rates (per 1000/year) ───────────────────────────────────
CBR = 30   # crude birth rate  — per 1000/year
CDR = 10   # crude death rate  — per 1000/year
# Net growth ≈ (CBR - CDR) / 1000 = 2%/yr → ~22% over 10 years

assert 1 <= CBR <= 60, f"CBR={CBR} not in expected per-1000/year range"
assert 1 <= CDR <= 60, f"CDR={CDR} not in expected per-1000/year range"

# ── Initial compartment fractions ──────────────────────────────────────────
F_S, F_E, F_I, F_R = 0.90, 0.00, 0.01, 0.09
assert abs(F_S + F_E + F_I + F_R - 1.0) < 1e-9, "Fractions must sum to 1"

# ── Compute initial counts ─────────────────────────────────────────────────
I_init = np.maximum(1, np.round(F_I * POPS).astype(np.int32))
R_init = np.round(F_R * POPS).astype(np.int32)
E_init = np.zeros(NNODES, dtype=np.int32)
S_init = POPS - E_init - I_init - R_init

assert np.all(S_init + E_init + I_init + R_init == POPS)
assert np.all(I_init > 0)
assert np.all(S_init >= 0)

# ── Scenario GeoDataFrame ──────────────────────────────────────────────────
PATCH_LON = np.array([0.000, 0.675, 1.351, 2.027])
scenario = gpd.GeoDataFrame({
    "nodeid":     np.arange(NNODES),
    "name":       [f"patch_{i}" for i in range(NNODES)],
    "population": POPS,
    "S": S_init, "E": E_init, "I": I_init, "R": R_init,
    "geometry":   [Point(lon, 0.0) for lon in PATCH_LON],
}, crs="EPSG:4326")

# ── Vital dynamics rate arrays — shape (NTICKS, NNODES), per-1000/year ─────
birthrates     = np.full((NTICKS, NNODES), float(CBR), dtype=np.float64)
mortalityrates = np.full((NTICKS, NNODES), float(CDR), dtype=np.float64)

# Sanity-check units (the #1 silent failure is passing daily per-capita here)
assert np.all(birthrates   >= 1) and np.all(birthrates   <= 60), \
    f"Birthrates out of per-1000/year range: {birthrates.min():.4f}–{birthrates.max():.4f}"
assert np.all(mortalityrates >= 1) and np.all(mortalityrates <= 60), \
    f"Mortality rates out of per-1000/year range: {mortalityrates.min():.4f}–{mortalityrates.max():.4f}"

# ── Age pyramid for newborns ───────────────────────────────────────────────
# Single-bin distribution: every newborn is assigned age 0 days.
# (MortalityByCDR does not use age, so exact pyramid shape is not critical here.)
pyramid = AliasedDistribution(np.array([1.0]))

# ── Explicit capacity pre-allocation for 10 years ─────────────────────────
# calc_capacity projects population growth over NTICKS to size the LaserFrame.
# safety_factor=2.0 gives headroom for stochastic variation.
# Model.__init__() also calls this internally when birthrates are provided —
# the explicit call here is to surface the numbers before the model is built.
capacity = calc_capacity(birthrates, POPS, safety_factor=2.0)
print("Pre-allocated capacity per patch (safety_factor=2.0):")
for i in range(NNODES):
    pct = 100.0 * (capacity[i] / POPS[i] - 1.0)
    print(f"  Patch {i} (N₀={POPS[i]:>7,}):  capacity = {capacity[i]:>7,}  "
          f"({pct:+.0f}% headroom)")
print(f"  Total pre-allocated: {capacity.sum():,}")
print()

# ── Model parameters ───────────────────────────────────────────────────────
params = PropertySet({
    "prng_seed":             SEED,
    "nticks":                NTICKS,
    "beta":                  BETA,
    "capacity_safety_factor": 2.0,  # passed to Model's internal calc_capacity call
})

# ── Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.gamma(shape=16, scale=0.25)   # E→I: mean = 4 d
infdurdist = dists.gamma(shape=40, scale=0.25)   # I→R: mean = 10 d

# ── Gravity coupling network ───────────────────────────────────────────────
PATCH_X_KM = np.array([0.0, 75.0, 150.0, 225.0])
GRAV_K, GRAV_A, GRAV_B, GRAV_C = 0.01, 1.0, 1.0, 1.5
MAX_EXPORT = 0.15

DIST_KM = np.abs(PATCH_X_KM[:, None] - PATCH_X_KM[None, :])
np.fill_diagonal(DIST_KM, np.inf)

W = (GRAV_K
     * np.power(POPS[:, None].astype(np.float64), GRAV_A)
     * np.power(POPS[None, :].astype(np.float64), GRAV_B)
     / np.power(DIST_KM, GRAV_C))
np.fill_diagonal(W, 0.0)

network = W / POPS[:, None].astype(np.float64)

row_sums = network.sum(axis=1)
scale = np.minimum(1.0, MAX_EXPORT / row_sums)
network = network * scale[:, None]

assert np.all(network >= 0),                        "Negative network weight"
assert np.all(network.sum(axis=1) <= MAX_EXPORT + 1e-9), \
    "Export fraction exceeds 15% cap"
assert network.sum() > 0,                           "Network is all zeros"

print("Gravity network row sums (export fractions):")
for i in range(NNODES):
    row = network[i]
    print(f"  Patch {i} → {[f'{row[j]:.4f}' for j in range(NNODES)]}  "
          f"(total export = {row.sum():.4f})")
print()

# ── Seasonal forcing profile — tiled across 10 years ──────────────────────
# Piecewise-linear profile (one year):
#   Days   0– 90: 1.3x   (winter peak)
#   Days  90–150: 1.3→0.7 (spring transition)
#   Days 150–240: 0.7x   (summer trough)
#   Days 240–365: 0.7→1.3 (autumn/winter transition)
_d = np.arange(365, dtype=np.float64)
season_365 = np.where(
    _d < 90,   1.3,
    np.where(
        _d < 150,  1.3 + (0.7 - 1.3) * (_d - 90.0)  / (150.0 - 90.0),
        np.where(
            _d < 240, 0.7,
            0.7 + (1.3 - 0.7) * (_d - 240.0) / (365.0 - 240.0)
        )
    )
)

assert abs(season_365.mean() - 1.0) < 0.01, \
    f"Seasonal profile mean={season_365.mean():.4f}, must be ~1.0"
assert abs(season_365[0]   - 1.3) < 1e-9, "Day 0 should be 1.3x"
assert abs(season_365[89]  - 1.3) < 1e-9, "Day 89 should be 1.3x"
assert abs(season_365[150] - 0.7) < 1e-9, "Day 150 should be 0.7x"
assert abs(season_365[239] - 0.7) < 1e-9, "Day 239 should be 0.7x"
assert np.all(season_365 > 0), "Negative seasonality values"

# Tile to cover all 10 simulation years
assert NTICKS % 365 == 0, "NTICKS must be an integer multiple of 365"
season_tiled = np.tile(season_365, NTICKS // 365)   # shape: (3650,)
seasonality  = ValuesMap.from_timeseries(season_tiled, NNODES)

print(f"Seasonal profile: mean={season_365.mean():.4f}, "
      f"peak={season_365.max():.2f}x (days 0–90 each year), "
      f"trough={season_365.min():.2f}x (days 150–240 each year)")
print()

# ── Construct model ────────────────────────────────────────────────────────
# Passing birthrates here triggers Model.__init__() to call calc_capacity
# internally using capacity_safety_factor from params.
model = Model(scenario, params, birthrates=birthrates)
model.network = network

# ── Assemble components ────────────────────────────────────────────────────
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model,    birthrates=birthrates,         pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortalityrates),
]

# ── Run ────────────────────────────────────────────────────────────────────
print(f"Running gravity 4-patch SEIR model — seasonal forcing + vital dynamics")
print(f"  β={BETA:.3f} d⁻¹,  R₀={R0},  latent={LATENT_D}d,  infectious={INFECT_D}d")
print(f"  Gravity: k={GRAV_K}, a={GRAV_A}, b={GRAV_B}, c={GRAV_C}")
print(f"  Spacing: 75 km;  max export fraction: {MAX_EXPORT:.0%}")
print(f"  Patches (pop): {POPS.tolist()}")
print(f"  Init: {F_S:.0%} S / {F_E:.0%} E / {F_I:.0%} I / {F_R:.0%} R")
print(f"  Seasonal: winter 1.3x (d0–90/yr), summer 0.7x (d150–240/yr)")
print(f"  Vital: CBR={CBR}/1000/yr, CDR={CDR}/1000/yr "
      f"(net ~{(CBR - CDR) / 10:.0f}%/yr)")
model.run("Gravity 4-Patch SEIR — Seasonal + Vital Dynamics")
print("Simulation complete.\n")

# ── Extract compartment time series ───────────────────────────────────────
years = np.arange(NTICKS) / 365.0
S = model.nodes.S[:NTICKS, :]
E = model.nodes.E[:NTICKS, :]
I = model.nodes.I[:NTICKS, :]
R = model.nodes.R[:NTICKS, :]
N = S + E + I + R

newly_infected = model.nodes.newly_infected[:NTICKS, :].sum(axis=1)

# ── Population growth check ────────────────────────────────────────────────
pop_start = int(N[0].sum())
pop_end   = int(N[-1].sum())
observed_annual_growth = (pop_end / pop_start) ** (1.0 / 10.0) - 1.0
expected_annual_growth = (CBR - CDR) / 1000.0   # 0.020

print("Population at start and end:")
print(f"  Day    0:  {pop_start:>10,}")
print(f"  Day {NTICKS}:  {pop_end:>10,}")
print(f"  Change:    {pop_end - pop_start:>+10,}  "
      f"({100.0 * (pop_end / pop_start - 1.0):+.1f}% over 10 yr)")
print(f"  Observed annual growth: {100.0 * observed_annual_growth:.2f}%  "
      f"(expected ~{100.0 * expected_annual_growth:.2f}% from "
      f"CBR–CDR = {CBR}–{CDR} = {CBR - CDR}/1000/yr)")
print()

# ── Validation ─────────────────────────────────────────────────────────────
assert np.all(S >= 0), "Negative S"
assert np.all(E >= 0), "Negative E"
assert np.all(I >= 0), "Negative I"
assert np.all(R >= 0), "Negative R"
assert np.all(N >  0), "Zero or negative N"

# Population should grow (CBR > CDR)
assert pop_end > pop_start, \
    f"Population shrank: start={pop_start:,}, end={pop_end:,}"

# Epidemic should grow beyond seed in every patch
for p in range(NNODES):
    assert I[:, p].max() > I_init[p], \
        f"No epidemic growth in patch {p} (pop={POPS[p]:,})"

# Observed growth should be within 0.5 pp of expected
assert abs(observed_annual_growth - expected_annual_growth) < 0.005, (
    f"Annual growth {100*observed_annual_growth:.2f}% deviates more than 0.5 pp "
    f"from expected {100*expected_annual_growth:.2f}%"
)

print("All validation checks passed.\n")

# ── Summary statistics ─────────────────────────────────────────────────────
print("Summary:")
print(f"  Peak total infectious : {I.sum(axis=1).max():,}")
print()
for p in range(NNODES):
    peak_ip  = I[:, p].max()
    day_peak = int(I[:, p].argmax())
    n_end    = int(N[-1, p])
    print(f"  Patch {p} (N₀={POPS[p]//1000}k → N₁₀={n_end//1000}k): "
          f"peak I = {peak_ip:>6,}  day {day_peak:4d}")

# ── Figure 1: SEIR compartments per patch (10 years) ──────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
fig.suptitle(
    f"4-Patch SEIR — Gravity + Seasonal Forcing + Vital Dynamics (10 years)\n"
    f"R₀ = {R0},  β = {BETA:.2f} d⁻¹,  "
    f"Latent = {LATENT_D} d,  Infectious = {INFECT_D} d  |  "
    f"CBR = {CBR}, CDR = {CDR} /1000/yr",
    fontsize=11,
)

STYLES = {
    "S": dict(color="steelblue",  lw=1.5, ls="-",  label="S"),
    "E": dict(color="darkorange", lw=1.0, ls="--", label="E"),
    "I": dict(color="firebrick",  lw=1.5, ls="-",  label="I"),
    "R": dict(color="seagreen",   lw=1.0, ls=":",  label="R"),
}

for p, ax in enumerate(axes.flat):
    for key, arr in (("S", S), ("E", E), ("I", I), ("R", R)):
        ax.plot(years, arr[:, p], **STYLES[key])
    ax.set_title(
        f"Patch {p}  (N₀ = {POPS[p]:,}  →  N₁₀ = {int(N[-1, p]):,})", fontsize=10
    )
    ax.set_ylabel("Count")
    ax.legend(fontsize=8, loc="center right")
    ax.grid(alpha=0.3)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):,}")
    )

for ax in axes[1]:
    ax.set_xlabel("Year")

plt.tight_layout()
plt.savefig("seir_4patch_vital_dynamics.png", dpi=120)
print("\nPlot saved → seir_4patch_vital_dynamics.png")
plt.close()

# ── Figure 2: Aggregate daily incidence (10 years) ────────────────────────
fig2, ax2 = plt.subplots(figsize=(12, 4))
ax2.fill_between(years, newly_infected, alpha=0.25, color="firebrick")
ax2.plot(years, newly_infected, color="firebrick", lw=1.5, label="New infections")
ax2.set_title(
    f"Daily New Infections — All 4 Patches Combined  "
    f"(R₀ = {R0}, Gravity + Seasonal + Vital Dynamics, 10 yr)",
    fontsize=11,
)
ax2.set_xlabel("Year")
ax2.set_ylabel("New infections per day")
ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax2.legend(fontsize=9)
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("seir_4patch_gravity_incidence.png", dpi=120)
print("Plot saved → seir_4patch_gravity_incidence.png")
plt.close()

# ── Figure 3: Gravity network heatmap ─────────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(5, 4))
im = ax3.imshow(network, cmap="YlOrRd", vmin=0, vmax=MAX_EXPORT)
plt.colorbar(im, ax=ax3, label="Outflow fraction")
ax3.set_xticks(range(NNODES))
ax3.set_yticks(range(NNODES))
ax3.set_xticklabels([f"P{i}\n({POPS[i]//1000}k)" for i in range(NNODES)])
ax3.set_yticklabels([f"P{i} ({POPS[i]//1000}k)" for i in range(NNODES)])
ax3.set_title(
    f"Gravity Network (row-normalised, cap={MAX_EXPORT:.0%})\n"
    f"k={GRAV_K}, a={GRAV_A}, b={GRAV_B}, c={GRAV_C}",
    fontsize=11,
)
ax3.set_xlabel("Destination patch")
ax3.set_ylabel("Source patch")
for i in range(NNODES):
    for j in range(NNODES):
        ax3.text(j, i, f"{network[i, j]:.3f}", ha="center", va="center",
                 fontsize=9,
                 color="black" if network[i, j] < MAX_EXPORT * 0.7 else "white")
plt.tight_layout()
plt.savefig("seir_4patch_gravity_network.png", dpi=120)
print("Plot saved → seir_4patch_gravity_network.png")
plt.close()

# ── Figure 4: Population trajectory (10 years) ────────────────────────────
fig4, ax4 = plt.subplots(figsize=(10, 4))
patch_colors = ["steelblue", "darkorange", "seagreen", "firebrick"]
for p in range(NNODES):
    ax4.plot(years, N[:, p] / 1000, lw=1.8, color=patch_colors[p],
             label=f"Patch {p} (N₀={POPS[p]//1000}k → N₁₀={int(N[-1, p])//1000}k)")
ax4.plot(years, N.sum(axis=1) / 1000, lw=2.5, ls="--", color="black",
         label=f"Total (N₀={pop_start//1000}k → N₁₀={pop_end//1000}k)")
ax4.set_title(
    f"Population Trajectory — Vital Dynamics  "
    f"(CBR={CBR}, CDR={CDR} /1000/yr)\n"
    f"Start: {pop_start:,}  →  End: {pop_end:,}  "
    f"({100*(pop_end/pop_start - 1):+.1f}% over 10 yr,  "
    f"~{100*observed_annual_growth:.2f}%/yr observed  vs  "
    f"{100*expected_annual_growth:.2f}%/yr expected)",
    fontsize=10,
)
ax4.set_xlabel("Year")
ax4.set_ylabel("Population (thousands)")
ax4.legend(fontsize=9)
ax4.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("seir_4patch_vital_dynamics_pop.png", dpi=120)
print("Plot saved → seir_4patch_vital_dynamics_pop.png")
plt.close()

# ── Figure 5: Seasonal forcing profile (one representative year) ───────────
fig5, ax5 = plt.subplots(figsize=(10, 3))
_days_1y = np.arange(365)
ax5.axvspan(0,   90,  alpha=0.12, color="steelblue", label="Winter peak (1.3x)")
ax5.axvspan(90,  150, alpha=0.06, color="gray",      label="Transition")
ax5.axvspan(150, 240, alpha=0.12, color="tomato",    label="Summer trough (0.7x)")
ax5.axvspan(240, 365, alpha=0.06, color="gray")
ax5.plot(_days_1y, season_365, color="black", lw=2.0, label="β multiplier")
ax5.axhline(1.0, color="gray", lw=1.0, ls="--", label="Baseline (1.0x)")
ax5.set_title(
    "Seasonal Forcing Profile (ValuesMap, repeated annually)  —  "
    "Winter 1.3x | Summer 0.7x | Linear transitions",
    fontsize=12,
)
ax5.set_xlabel("Day of year")
ax5.set_ylabel("β multiplier")
ax5.set_ylim(0.5, 1.5)
ax5.set_xlim(0, 364)
ax5.legend(fontsize=9, loc="upper right")
ax5.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("seir_4patch_seasonal_profile.png", dpi=120)
print("Plot saved → seir_4patch_seasonal_profile.png")
plt.close()
```

**Key changes from the 1-year script:**

| What changed | Detail |
|---|---|
| `NTICKS` | `365` → `10 * 365` (3650) |
| New imports | `AliasedDistribution`, `calc_capacity`, `BirthsByCBR`, `MortalityByCDR` |
| Rate arrays | `birthrates` and `mortalityrates` as `(3650, 4)` float64, value 30 / 10 per-1000/yr |
| `pyramid` | `AliasedDistribution(np.array([1.0]))` — single bin, all newborns age 0 |
| `calc_capacity` | Called explicitly before `Model()` to print capacity; also called internally by `Model.__init__()` via `capacity_safety_factor=2.0` in params |
| `Model()` | `birthrates=birthrates` passed to trigger internal capacity pre-allocation |
| Components | `BirthsByCBR` and `MortalityByCDR` appended after `Transmission` |
| Seasonal profile | `np.tile(season_365, 10)` → 3650-element array before `ValuesMap.from_timeseries` |
| Removed assertion | `N[:, p] == POPS[p]` — population is no longer constant |
| Population print | Start/end totals + observed vs expected ~2%/yr growth |
| Output figures | x-axis in years; new Figure 4 shows population trajectory |

The expected console output for the population check will look like:
```
Population at start and end:
  Day    0:     530,000
  Day 3650:     ~645,000
  Change:       +~115,000  (+~21.7% over 10 yr)
  Observed annual growth: ~2.00%  (expected ~2.00% from CBR–CDR = 30–10 = 20/1000/yr)
```
