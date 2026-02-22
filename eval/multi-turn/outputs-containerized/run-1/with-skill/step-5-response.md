Here is the complete script with vital dynamics added:

```python
"""
4-patch SEIR model with gravity-model spatial coupling, seasonal forcing,
and vital dynamics: births (CBR=30/1000/yr) + deaths (CDR=10/1000/yr).
Extended to 10 years. (laser-generic v1.0.0)

Disease parameters:
  R0 = 5, latent period = 4 days, infectious period = 10 days

Vital dynamics:
  CBR = 30 per 1000/year (LASER unit convention)
  CDR = 10 per 1000/year
  Net growth = (30-10)/1000 = 2%/year → 10-yr factor ≈ 1.02^10 ≈ 1.219

  Newborns enter as susceptibles. Deaths drawn uniformly from all compartments.
  BirthsByCBR and MortalityByCDR both expect per-1000/year rates — passing
  daily per-capita values is a silent failure (no births, static population).

Capacity pre-allocation:
  calc_capacity(birthrates, populations, safety_factor=2.0) is called
  explicitly to show the headroom, and implicitly by Model.__init__() via
  the `birthrates` argument with capacity_safety_factor=2.0 in params.
  CBR=30 over 10 years → cumulative births ≈ 30% of initial N; ×2 safety = 60%.

Initial conditions:
  Patch 0 (source): 90% S | 1% I | 9% R
  Patches 1–3:     100% S

Spatial coupling: gravity model k=0.01, a=1, b=1, c=1.5, row-capped at 15%.
Seasonal forcing: winter peak 1.3×, summer trough 0.7×, cosine transitions.
Simulation: 10 years (3650 ticks).
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.core.demographics import AliasedDistribution
from laser.core.utils import calc_capacity
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
import matplotlib.pyplot as plt

# ── 1. Simulation parameters ───────────────────────────────────────────────────
R0              = 5.0
latent_days     = 4.0    # mean exposed (E) period
infectious_days = 10.0   # mean infectious (I) period

beta   = R0 / infectious_days   # 0.5 per day
nticks = 10 * 365               # 10 years

# Vital dynamics — LASER expects per-1000/year; see assertion below
CBR = 30.0   # crude birth rate
CDR = 10.0   # crude death rate

print(f"beta = {beta:.4f}  (R0={R0}, D_I={infectious_days}d)")
print(f"CBR={CBR:.0f}/1000/yr, CDR={CDR:.0f}/1000/yr  "
      f"→  net growth = {(CBR - CDR) / 10:.1f}%/yr")

# Guard against accidental daily-rate units (the #1 silent failure)
assert 1 <= CBR <= 60, f"CBR={CBR} looks wrong — must be per-1000/yr (range 10–50)"
assert 1 <= CDR <= 60, f"CDR={CDR} looks wrong — must be per-1000/yr (range 5–20)"

# ── 2. Patch populations and initial conditions ────────────────────────────────
populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
nnodes      = len(populations)

seed_frac = 0.01   # 1% initially infectious in source patch

I_init = np.zeros(nnodes, dtype=np.int32)
E_init = np.zeros(nnodes, dtype=np.int32)
R_init = np.zeros(nnodes, dtype=np.int32)

# Patch 0: 1% I (active seed) + 9% R (prior immunity)
I_init[0] = int(round(seed_frac * populations[0]))
R_init[0] = int(round(0.09    * populations[0]))

S_init = (populations - I_init - E_init - R_init).astype(np.int32)

assert np.all(S_init + E_init + I_init + R_init == populations), \
    "S + E + I + R must equal N"
assert np.all(S_init >= 0), "Negative S"
assert I_init[0] > 0, "Source patch must have infectious agents"

print(f"\nInitial total population: {populations.sum():,}")
print(f"  Patch 0 — S={S_init[0]:,} ({S_init[0]/populations[0]:.0%}), "
      f"I={I_init[0]:,} ({I_init[0]/populations[0]:.1%}), "
      f"R={R_init[0]:,} ({R_init[0]/populations[0]:.0%})")
print(f"  Patches 1–3 — 100% S")

# ── 3. Gravity-model spatial network ──────────────────────────────────────────
positions_km = np.array([0.0, 75.0, 150.0, 225.0])

dist_km = np.abs(positions_km[:, None] - positions_km[None, :])
np.fill_diagonal(dist_km, np.inf)

gravity_k = 0.01
gravity_a = 1.0
gravity_b = 1.0
gravity_c = 1.5

F_raw = (
    gravity_k
    * populations[:, None] ** gravity_a
    * populations[None, :] ** gravity_b
    / dist_km ** gravity_c
)
np.fill_diagonal(F_raw, 0.0)

export_frac_matrix = F_raw / populations[:, None]

max_export_frac = 0.15
export_total    = export_frac_matrix.sum(axis=1)
cap_scale       = np.where(
    export_total > max_export_frac,
    max_export_frac / export_total,
    1.0,
)
network = export_frac_matrix * cap_scale[:, None]

assert np.all(network >= 0), "Negative network weight"
assert network.sum(axis=1).max() <= max_export_frac + 1e-10, "Export cap violated"
assert network.sum() > 0, "Network is all zeros"

row_sums = network.sum(axis=1)
print(f"\nGravity network row sums (FOI export fraction): {row_sums.round(4)}")

# ── 4. Seasonal forcing profile (1 year, tiled to nticks) ─────────────────────
season_hi = 1.3
season_lo = 0.7

season_365 = np.ones(365, dtype=np.float64)
season_365[0:91]    = season_hi    # days   0–90  winter peak
season_365[150:241] = season_lo    # days 150–240 summer trough

n_down = 150 - 91                  # 59 steps: ramp winter→summer
season_365[91:150] = (
    0.5 * (season_hi + season_lo)
    + 0.5 * (season_hi - season_lo)
    * np.cos(np.pi * np.arange(n_down) / (n_down - 1))
)

n_up = 365 - 241                   # 124 steps: ramp summer→winter
season_365[241:365] = (
    0.5 * (season_hi + season_lo)
    - 0.5 * (season_hi - season_lo)
    * np.cos(np.pi * np.arange(n_up) / (n_up - 1))
)

assert abs(season_365.mean() - 1.0) < 0.02, \
    f"Seasonal profile not normalized: mean={season_365.mean():.4f}"

# Tile the 365-day profile across all 10 years
season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality  = ValuesMap.from_timeseries(season_tiled, nnodes)
print(f"\nSeasonal profile: mean={season_365.mean():.4f}, "
      f"min={season_365.min():.2f}, max={season_365.max():.2f}  "
      f"→ tiled ValuesMap shape {seasonality.values.shape}")

# ── 5. Vital dynamics rate arrays ──────────────────────────────────────────────
# Shape (nticks, nnodes), constant rates; BirthsByCBR/MortalityByCDR expect
# per-1000/year — the framework divides internally.
birthrates     = np.full((nticks, nnodes), CBR, dtype=np.float64)
mortalityrates = np.full((nticks, nnodes), CDR, dtype=np.float64)

# ── 6. Explicit calc_capacity call — shows pre-allocated slots ────────────────
# CBR=30 over 10 years → each patch births ≈ 30% of initial N cumulatively.
# safety_factor=2 doubles that headroom to guard against stochastic spikes.
capacity_per_node = calc_capacity(birthrates, populations, safety_factor=2.0)

print(f"\ncalc_capacity (CBR={CBR:.0f}/1000/yr × {nticks//365} yr, safety_factor=2.0):")
for i in range(nnodes):
    print(f"  Patch {i}: initial={populations[i]:>9,}  "
          f"capacity={capacity_per_node[i]:>9,}  "
          f"({capacity_per_node[i] / populations[i]:.2f}×)")
print(f"  Total:   initial={populations.sum():>9,}  "
      f"capacity={capacity_per_node.sum():>9,}  "
      f"({capacity_per_node.sum() / populations.sum():.2f}×)")

# ── 7. Age pyramid for BirthsByCBR ────────────────────────────────────────────
# Stable exponential age distribution: P(alive at age a days) ∝ exp(-CDR/1000/365 × a).
# This approximates the equilibrium demography under constant mortality.
cdr_daily  = CDR / 1000.0 / 365.0
ages_days  = np.arange(100 * 365, dtype=np.float64)
stable_age = np.exp(-cdr_daily * ages_days)
pyramid    = AliasedDistribution(stable_age)

# ── 8. Scenario GeoDataFrame ───────────────────────────────────────────────────
lons = np.full(nnodes, 70.0)
lats = 30.0 + positions_km / 111.0

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(nnodes, dtype=np.int32),
        "name":       [f"Patch_{i}" for i in range(nnodes)],
        "population": populations,
        "S":          S_init,
        "E":          E_init,
        "I":          I_init,
        "R":          R_init,
        "geometry":   [Point(lon, lat) for lon, lat in zip(lons, lats)],
    },
    crs="EPSG:4326",
)

# ── 9. Model parameters ────────────────────────────────────────────────────────
params = PropertySet(
    {
        "prng_seed":              42,
        "nticks":                 nticks,
        "beta":                   beta,
        # capacity_safety_factor is read by Model.__init__() when it calls
        # calc_capacity internally (same value as the explicit call above).
        "capacity_safety_factor": 2.0,
    }
)

# ── 10. Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.gamma(shape=16, scale=0.25)   # latent    mean = 4 d
infdurdist = dists.gamma(shape=40, scale=0.25)   # infectious mean = 10 d

# ── 11. Build model ─────────────────────────────────────────────────────────────
# Passing birthrates triggers Model.__init__() → calc_capacity for the correct
# pre-allocated LaserFrame capacity (guided by capacity_safety_factor in params).
model = Model(scenario, params, birthrates=birthrates)

# Overwrite default network with the manually-constructed gravity network
model.network = network
assert model.network.sum() > 0, "Network must be non-zero"

# ── 12. Assemble components ─────────────────────────────────────────────────────
# Ordering: SEIR transitions → Transmission → Births → Deaths.
# BirthsByCBR adds new susceptibles each tick; MortalityByCDR removes agents.
# Note: MortalityByCDR parameter is `mortalityrates=`, not `deathrates=`.
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model, birthrates=birthrates, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortalityrates),
]

# ── 13. Run ────────────────────────────────────────────────────────────────────
print(f"\nRunning {nticks}-tick ({nticks // 365}-year) simulation ...")
model.run("4-Patch SEIR + Gravity + Seasonality + Vital Dynamics (10yr)")
print("Done.\n")

# ── 14. Extract output arrays ──────────────────────────────────────────────────
S = model.nodes.S[:nticks, :]    # shape (nticks, nnodes)
E = model.nodes.E[:nticks, :]
I = model.nodes.I[:nticks, :]
R = model.nodes.R[:nticks, :]
N = S + E + I + R

# ── 15. Population growth check ────────────────────────────────────────────────
pop_start = int(N[0, :].sum())
pop_end   = int(N[-1, :].sum())
years     = nticks / 365.0

theoretical_end = pop_start * ((1.0 + (CBR - CDR) / 1000.0) ** years)
observed_annual = (pop_end / pop_start) ** (1.0 / years) - 1.0
target_annual   = (CBR - CDR) / 1000.0   # 0.02 (2%)

print("─" * 56)
print(f"  Population at start  (day     0):  {pop_start:>12,}")
print(f"  Population at end    (day {nticks:>5}):  {pop_end:>12,}")
print(f"  Theoretical end pop (2%/yr × 10yr): {theoretical_end:>12,.0f}")
print(f"  Ratio simulated/theoretical:        {pop_end / theoretical_end:>12.4f}")
print(f"")
print(f"  Observed annual growth rate:  {observed_annual * 100:>7.3f}%")
print(f"  Target  annual growth rate:   {target_annual * 100:>7.3f}%")
print(f"  Difference:                   {(observed_annual - target_annual) * 100:>+7.3f} pp")
print("─" * 56)

# ── 16. Verification ───────────────────────────────────────────────────────────
assert np.all(S >= 0) and np.all(E >= 0) and np.all(I >= 0) and np.all(R >= 0), \
    "Negative compartment count"

# Growth rate within ±0.5 percentage points of 2%/yr
assert abs(observed_annual - target_annual) < 0.005, (
    f"Annual growth {observed_annual * 100:.3f}% deviates from "
    f"target {target_annual * 100:.3f}% by more than 0.5 pp"
)

cumulative_infections = model.nodes.newly_infected[:nticks, :].sum(axis=0)
assert np.all(cumulative_infections > 0), (
    f"No infections in patch(es): {np.where(cumulative_infections == 0)[0]}"
)

print("\nVerification passed.")

# Epidemic summary
print(f"\n{'Patch':<10} {'N start':>12} {'N end':>12} "
      f"{'Total Infected':>16} {'Attack Rate':>12}")
print("-" * 64)
for i in range(nnodes):
    ar = cumulative_infections[i] / populations[i]
    print(
        f"Patch {i:<4}  {populations[i]:>12,}  {N[-1, i]:>11,}  "
        f"{cumulative_infections[i]:>14,.0f}  {ar:>11.1%}"
    )
print("-" * 64)
print(f"{'Total':<10} {pop_start:>12,}  {pop_end:>11,}  "
      f"{cumulative_infections.sum():>14,.0f}  "
      f"{cumulative_infections.sum() / pop_start:>11.1%}")

# ── 17. Diagnostic plots ───────────────────────────────────────────────────────
days         = np.arange(nticks)
patch_labels = [f"Patch {i}  (N₀={populations[i] // 1000}k)" for i in range(nnodes)]
colors       = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

fig, axes = plt.subplots(2, 3, figsize=(18, 8))
fig.suptitle(
    f"4-Patch SEIR + Gravity + Seasonality + Vital Dynamics  |  "
    f"R₀={R0}, β={beta}  |  "
    f"CBR={CBR:.0f}/CDR={CDR:.0f} per 1000/yr (net +2%/yr)  |  "
    f"{nticks // 365}-year simulation  |  seed: patch 0 only",
    fontsize=9,
)

# Panel A: total population trajectory vs theoretical compound growth
ax = axes[0, 0]
total_N   = N.sum(axis=1)
theory_N  = pop_start * ((1.0 + (CBR - CDR) / 1000.0) ** (days / 365.0))
ax.plot(days, total_N,  color="tab:blue", linewidth=1.5, label="Simulated")
ax.plot(days, theory_N, color="tomato", linestyle="--", linewidth=1.5,
        label=f"Theoretical (2%/yr)")
ax.set_title("Total Population: Simulated vs Theoretical")
ax.set_xlabel("Day")
ax.set_ylabel("Population")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel B: per-patch population trajectories
ax = axes[0, 1]
for i in range(nnodes):
    ax.plot(days, N[:, i], label=patch_labels[i], color=colors[i])
ax.set_title("Population per Patch")
ax.set_xlabel("Day")
ax.set_ylabel("Population")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel C: susceptible fraction S/N by patch
ax = axes[0, 2]
for i in range(nnodes):
    ax.plot(days, S[:, i] / N[:, i], label=patch_labels[i], color=colors[i])
ax.set_title("Susceptible Fraction  S/N  by Patch")
ax.set_xlabel("Day")
ax.set_ylabel("S / N")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel D: aggregate SEIR fractions (all patches combined)
ax = axes[1, 0]
N_total = N.sum(axis=1)
ax.plot(days, S.sum(axis=1) / N_total, label="S", color="tab:blue")
ax.plot(days, E.sum(axis=1) / N_total, label="E", color="tab:orange")
ax.plot(days, I.sum(axis=1) / N_total, label="I", color="tab:red")
ax.plot(days, R.sum(axis=1) / N_total, label="R", color="tab:green")
ax.set_title("Aggregate SEIR Dynamics  (All Patches)")
ax.set_xlabel("Day")
ax.set_ylabel("Fraction of Population")
ax.legend()
ax.grid(alpha=0.3)

# Panel E: daily new infections by patch
ax = axes[1, 1]
for i in range(nnodes):
    ax.plot(days, model.nodes.newly_infected[:nticks, i],
            label=patch_labels[i], color=colors[i], alpha=0.8)
ax.set_title("Daily New Infections by Patch")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel F: total daily incidence (all patches, 7-day smoothed overlay)
ax = axes[1, 2]
total_inc = model.nodes.newly_infected[:nticks, :].sum(axis=1)
ax.plot(days, total_inc, color="black", alpha=0.4, linewidth=0.8, label="Daily")
smoothed  = np.convolve(total_inc, np.ones(7) / 7, mode="same")
ax.plot(days, smoothed, color="black", linewidth=1.5, label="7-day avg")
ax.set_title("Total Daily New Infections  (All Patches)")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("seir_4patch_vital_dynamics.png", dpi=100, bbox_inches="tight")
plt.show()
print("\nPlot saved: seir_4patch_vital_dynamics.png")
```

Key changes from the 1-year script:

**Vital dynamics additions:**
- `birthrates` / `mortalityrates`: `(nticks, nnodes)` float64 arrays, constant at 30 and 10 per 1000/year. The `MortalityByCDR` parameter is `mortalityrates=`, not `deathrates=`.
- `AliasedDistribution(stable_age)`: stable exponential age pyramid (proportional to `exp(-CDR/1000/365 × age)`).
- `BirthsByCBR` and `MortalityByCDR` appended after `Transmission` in the component list.

**Capacity pre-allocation:**
- `calc_capacity(birthrates, populations, safety_factor=2.0)` is called explicitly to print the per-node slot headroom (≈1.6× initial N), then `capacity_safety_factor: 2.0` in `params` directs the same calculation inside `Model.__init__()`.
- `birthrates` is also passed to `Model(scenario, params, birthrates=birthrates)` so the constructor performs the internal pre-allocation.

**Expected output:** The ~2% annual growth assertion checks `|observed_rate − 0.020| < 0.005 pp`. With net growth = (30−10)/1000 = 0.020/yr and 530k agents, stochastic variation over 10 years is well within that tolerance.
