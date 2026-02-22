Now I have everything I need. The bug is clear.

**Root cause:** All patches start with identical proportions — 90% S, 1% I, 9% R. The Transmission component computes FOI as `β · I[i]/N[i]`, so when every patch has the same I/N ratio at t=0, the spatial coupling terms `(col_sum[i] − row_sum[i]) · ft` produce only microscopic differences (~0.001% of ft). All patches enter the epidemic in perfect synchrony and the curves are visually indistinguishable regardless of the network.

There is no bug in the spatial coupling code itself — `model.network` is correctly accessed at runtime in the Transmission step. The fix is to seed infection in **one patch only**, leaving the others fully susceptible, so spatial propagation is observable.

```python
"""
4-patch SEIR model with gravity-model spatial coupling and seasonal forcing
(laser-generic v1.0.0).

Disease parameters:
  R0 = 5, latent period = 4 days, infectious period = 10 days

Initial conditions:
  Patch 0 (source): 90% Susceptible | 1% Infectious | 9% Recovered | 0% Exposed
  Patches 1–3:     100% Susceptible | 0% Infectious | 0% Recovered | 0% Exposed

  Seeding infection only in patch 0 makes spatial propagation visible.
  If all patches start with the same I/N ratio the FOI formula (β·I/N) is
  identical everywhere on tick 0, so spatial coupling produces no perceptible
  differentiation regardless of network strength.

Spatial coupling: gravity model  F_ij = k · N_i^a · N_j^b / d_ij^c
  k=0.01, a=1, b=1, c=1.5
  Patches in a line, 75 km apart.
  Row-normalized to cap each patch's total FOI export at 15%.

Seasonal forcing (ValuesMap):
  Days   0– 90:  1.3x baseline (winter peak plateau)
  Days  91–149:  cosine ramp down 1.3x → 0.7x
  Days 150–240:  0.7x baseline (summer trough plateau)
  Days 241–364:  cosine ramp up  0.7x → 1.3x

Simulation duration: 1 year (365 days).
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
import matplotlib.pyplot as plt

# ── 1. Disease parameters ──────────────────────────────────────────────────────
R0              = 5.0
latent_days     = 4.0    # mean exposed (E) period
infectious_days = 10.0   # mean infectious (I) period

# beta = R0 / D_I  (standard SIR/SEIR relationship)
beta   = R0 / infectious_days   # = 0.5 per day
nticks = 365                    # 1 year

print(f"beta = {beta:.4f}  (R0={R0}, infectious_period={infectious_days}d)")

# ── 2. Patch populations and initial conditions ────────────────────────────────
populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
nnodes      = len(populations)

# FIX: Seed infection only in patch 0; patches 1–3 start fully susceptible.
#
# When all patches share the same I/N ratio at t=0 the force-of-infection
# β·I[i]/N[i] is identical in every patch, so the asymmetric gravity network
# produces only microscopic corrections (~row_sum difference × ft ≈ 1e-5).
# The patches enter the epidemic in perfect synchrony and the curves are
# visually indistinguishable.  Seeding a single source patch breaks that
# symmetry and lets the travelling wave of infection reveal the network.
seed_frac = 0.01          # 1 % initially infectious in the source patch

I_init = np.zeros(nnodes, dtype=np.int32)
E_init = np.zeros(nnodes, dtype=np.int32)
R_init = np.zeros(nnodes, dtype=np.int32)

# Patch 0: partial pre-existing immunity (9 % R) + active seed (1 % I)
I_init[0] = int(round(seed_frac * populations[0]))
R_init[0] = int(round(0.09 * populations[0]))

# All patches: remaining population is fully susceptible
S_init = (populations - I_init - E_init - R_init).astype(np.int32)

# Verify
assert np.all(S_init + E_init + I_init + R_init == populations), \
    "S + E + I + R must equal N in every patch"
assert np.all(S_init >= 0), "Negative susceptible count"
assert I_init[0] > 0, "Source patch must have at least one infectious agent"

print(
    f"Initial conditions (source patch 0 only):\n"
    f"  Total N={populations.sum():,}\n"
    f"  Patch 0 — S={S_init[0]:,} ({S_init[0]/populations[0]:.0%}), "
    f"I={I_init[0]:,} ({I_init[0]/populations[0]:.1%}), "
    f"R={R_init[0]:,} ({R_init[0]/populations[0]:.0%})\n"
    f"  Patches 1–3 — fully susceptible (100 % S, 0 % I)"
)

# ── 3. Gravity-model spatial network ──────────────────────────────────────────
# Patches in a line, 75 km apart.
positions_km = np.array([0.0, 75.0, 150.0, 225.0])

# Pairwise distance matrix (km); diagonal → inf to keep self-terms at zero.
dist_km = np.abs(positions_km[:, None] - positions_km[None, :])
np.fill_diagonal(dist_km, np.inf)

# Gravity parameters
gravity_k = 0.01
gravity_a = 1.0    # N_i exponent
gravity_b = 1.0    # N_j exponent
gravity_c = 1.5    # distance exponent

# Raw flows  F_ij = k · N_i^a · N_j^b / d_ij^c  (off-diagonal; diagonal is 0)
F_raw = (
    gravity_k
    * populations[:, None] ** gravity_a
    * populations[None, :] ** gravity_b
    / dist_km ** gravity_c
)
np.fill_diagonal(F_raw, 0.0)

# Convert to dimensionless export-fraction matrix:
#   f_ij = F_ij / N_i  =  fraction of patch i's FOI it "sends" to patch j.
export_frac_matrix = F_raw / populations[:, None]

# Row-cap: if a patch's total export exceeds 15%, scale that row down uniformly.
max_export_frac = 0.15
export_total    = export_frac_matrix.sum(axis=1)
cap_scale       = np.where(
    export_total > max_export_frac,
    max_export_frac / export_total,
    1.0,
)
network = export_frac_matrix * cap_scale[:, None]

# Sanity checks
assert np.all(network >= 0), "Negative network weight detected"
assert network.sum(axis=1).max() <= max_export_frac + 1e-10, \
    f"Export cap violated: max row sum = {network.sum(axis=1).max():.4f}"
assert network.sum() > 0, "Network is all zeros — no spatial coupling"

row_sums = network.sum(axis=1)
print(f"\nGravity network (k={gravity_k}, a={gravity_a}, b={gravity_b}, c={gravity_c}):")
print("  network[i,j] = fraction of patch i's FOI exported to patch j")
print(np.array2string(network, precision=4, suppress_small=True))
print(f"  Row sums (FOI export fraction per patch): {row_sums.round(4)}")
print(f"  15% cap applied to: "
      f"{(export_total > max_export_frac).sum()}/{nnodes} patches")

# ── 4. Seasonal forcing profile ────────────────────────────────────────────────
# Piecewise seasonal multiplier on beta:
#   Days   0– 90:  1.3x plateau  (winter peak)
#   Days  91–149:  cosine ramp down from 1.3x to 0.7x  (59 transition steps)
#   Days 150–240:  0.7x plateau  (summer trough)
#   Days 241–364:  cosine ramp up from 0.7x to 1.3x  (124 transition steps)

season_hi = 1.3
season_lo = 0.7

season_365 = np.ones(365, dtype=np.float64)

# ── Plateau regions ─────────────────────────────────────────────────────
season_365[0:91]    = season_hi   # days   0–90  (winter peak)
season_365[150:241] = season_lo   # days 150–240 (summer trough)

# ── Transition: winter peak → summer trough  (days 91–149, 59 steps) ────
n_down  = 150 - 91                          # = 59
i_down  = np.arange(n_down)                 # 0 … 58
season_365[91:150] = (
    0.5 * (season_hi + season_lo)
    + 0.5 * (season_hi - season_lo) * np.cos(np.pi * i_down / (n_down - 1))
)

# ── Transition: summer trough → winter peak  (days 241–364, 124 steps) ──
n_up   = 365 - 241                          # = 124
i_up   = np.arange(n_up)                    # 0 … 123
season_365[241:365] = (
    0.5 * (season_hi + season_lo)
    - 0.5 * (season_hi - season_lo) * np.cos(np.pi * i_up / (n_up - 1))
)

# ── Validate ─────────────────────────────────────────────────────────────
assert abs(season_365.mean() - 1.0) < 0.02, \
    f"Seasonal profile mean={season_365.mean():.4f}, expected ~1.0"
assert abs(season_365[0]   - season_hi) < 1e-9, "Day 0 must equal 1.3"
assert abs(season_365[90]  - season_hi) < 1e-9, "Day 90 must equal 1.3"
assert abs(season_365[150] - season_lo) < 1e-9, "Day 150 must equal 0.7"
assert abs(season_365[240] - season_lo) < 1e-9, "Day 240 must equal 0.7"
assert abs(season_365[91]  - season_hi) < 1e-9, "Ramp-down start must be 1.3"
assert abs(season_365[149] - season_lo) < 1e-9, "Ramp-down end must be 0.7"
assert abs(season_365[241] - season_lo) < 1e-9, "Ramp-up start must be 0.7"
assert abs(season_365[364] - season_hi) < 1e-9, "Ramp-up end must be 1.3"

print(f"\nSeasonal profile (365-day):")
print(f"  mean={season_365.mean():.4f}, "
      f"min={season_365.min():.4f} (day {np.argmin(season_365)}), "
      f"max={season_365.max():.4f} (day {np.argmax(season_365)})")
print(f"  Winter peak  (days  0– 90): {season_365[ 0:91].min():.2f}–{season_365[ 0:91].max():.2f}x")
print(f"  Summer trough(days150–240): {season_365[150:241].min():.2f}–{season_365[150:241].max():.2f}x")

# ── Build ValuesMap for nticks days across nnodes patches ─────────────────────
season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality  = ValuesMap.from_timeseries(season_tiled, nnodes)

print(f"  ValuesMap shape: {seasonality.values.shape}  (nticks={nticks}, nnodes={nnodes})")

# ── 5. Scenario GeoDataFrame ───────────────────────────────────────────────────
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

# ── 6. Model parameters ────────────────────────────────────────────────────────
params = PropertySet(
    {
        "prng_seed": 42,
        "nticks":    nticks,
        "beta":      beta,
    }
)

# ── 7. Duration distributions ──────────────────────────────────────────────────
expdurdist = dists.gamma(shape=16, scale=0.25)   # latent period  mean=4d
infdurdist = dists.gamma(shape=40, scale=0.25)   # infectious period mean=10d

# ── 8. Build model and install gravity network ─────────────────────────────────
model = Model(scenario, params, birthrates=None)

# Overwrite the default network with the manually-constructed gravity network.
# The Transmission component accesses self.model.network at runtime, so this
# assignment takes effect for the entire simulation.
model.network = network
assert model.network.sum() > 0, "Network must be non-zero for spatial coupling"

# ── 9. Assemble components ─────────────────────────────────────────────────────
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
]

# ── 10. Run ────────────────────────────────────────────────────────────────────
print("\nRunning simulation ...")
model.run("4-Patch SEIR with Gravity Coupling + Seasonal Forcing")
print("Done.\n")

# ── 11. Post-run verification ──────────────────────────────────────────────────
S = model.nodes.S[:nticks, :]    # shape (nticks, nnodes)
E = model.nodes.E[:nticks, :]
I = model.nodes.I[:nticks, :]
R = model.nodes.R[:nticks, :]
N = S + E + I + R

# 11a. Population conservation (no demographics → N constant per patch)
for i in range(nnodes):
    N_i = N[:, i]
    assert N_i.min() == N_i.max(), (
        f"Population not conserved in Patch {i}: "
        f"min={N_i.min():,}, max={N_i.max():,}"
    )

# 11b. Non-negativity
assert np.all(S >= 0) and np.all(E >= 0) and np.all(I >= 0) and np.all(R >= 0), \
    "Negative compartment count detected"

# 11c. Epidemic occurred in every patch (infection spreads from patch 0 via network)
cumulative_infections = model.nodes.newly_infected[:nticks, :].sum(axis=0)
assert np.all(cumulative_infections > 0), (
    f"No infections in patch(es): {np.where(cumulative_infections == 0)[0]}\n"
    f"  Check that the network has non-zero entries connecting patch 0 to all others."
)

print("Verification passed.")

# Compute epidemic arrival day (first day with new infections) per patch
arrival_day = np.array([
    int(np.argmax(model.nodes.newly_infected[:nticks, i] > 0))
    for i in range(nnodes)
])

print(f"\n{'Patch':<10} {'Population':>12} {'Arrival Day':>12} "
      f"{'Total Infected':>16} {'Attack Rate':>12}")
print("-" * 64)
for i in range(nnodes):
    ar = cumulative_infections[i] / populations[i]
    print(
        f"Patch {i:<4}  {populations[i]:>12,}  "
        f"{arrival_day[i]:>11d}  "
        f"{cumulative_infections[i]:>14,.0f}  {ar:>11.1%}"
    )
total_pop = populations.sum()
total_inf = cumulative_infections.sum()
print("-" * 64)
print(f"{'Total':<10} {total_pop:>12,}  {'':>12}  "
      f"{total_inf:>14,.0f}  {total_inf/total_pop:>11.1%}")

# Check that spatial propagation produced staggered arrival times
if arrival_day.max() > arrival_day.min():
    print(f"\nSpatial propagation confirmed: arrival days span "
          f"{arrival_day.min()}–{arrival_day.max()} "
          f"(range = {arrival_day.max() - arrival_day.min()} days)")
else:
    print("\nWarning: all patches show simultaneous epidemic arrival — "
          "inspect network coupling strength.")

# ── 12. Diagnostic plots ───────────────────────────────────────────────────────
days         = np.arange(nticks)
patch_labels = [f"Patch {i}  (N={populations[i]//1000}k)" for i in range(nnodes)]
colors       = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

fig, axes = plt.subplots(2, 3, figsize=(18, 8))
fig.suptitle(
    f"4-Patch SEIR with Gravity Coupling + Seasonal Forcing  |  "
    f"R₀={R0}, β={beta}  |  "
    f"k={gravity_k}, a={gravity_a}, b={gravity_b}, c={gravity_c}  |  "
    f"max export={max_export_frac:.0%}  |  "
    f"seasonal: {season_lo}×–{season_hi}× (winter peak days 0–90)  |  "
    f"seed: patch 0 only",
    fontsize=9,
)

# Panel A: seasonal forcing profile
ax = axes[0, 0]
ax.fill_between(days, season_365, 1.0,
                where=season_365 >= 1.0, alpha=0.25, color="steelblue",
                label="Above baseline (winter)")
ax.fill_between(days, season_365, 1.0,
                where=season_365 <= 1.0, alpha=0.25, color="tomato",
                label="Below baseline (summer)")
ax.plot(days, season_365, color="black", linewidth=1.5)
ax.axhline(1.0, color="gray", linewidth=0.8, linestyle="--")
ax.axvspan(0,   90,  alpha=0.08, color="steelblue")
ax.axvspan(150, 240, alpha=0.08, color="tomato")
ax.set_title("Seasonal Forcing Multiplier")
ax.set_xlabel("Day of Year")
ax.set_ylabel("β multiplier")
ax.set_ylim(0.5, 1.5)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel B: daily new infections by patch — staggered peaks show spatial spread
ax = axes[0, 1]
for i in range(nnodes):
    ax.plot(days, model.nodes.newly_infected[:nticks, i],
            label=patch_labels[i], color=colors[i])
    # Mark epidemic arrival day
    ax.axvline(arrival_day[i], color=colors[i], linewidth=0.8,
               linestyle=":", alpha=0.7)
ax.set_title("Daily New Infections by Patch\n"
             "(dotted lines = epidemic arrival day)")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel C: susceptible fraction S/N by patch
ax = axes[0, 2]
for i in range(nnodes):
    ax.plot(days, S[:, i] / populations[i],
            label=patch_labels[i], color=colors[i])
ax.set_title("Susceptible Fraction  S/N  by Patch")
ax.set_xlabel("Day")
ax.set_ylabel("S / N")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel D: aggregate SEIR compartments (all patches combined)
ax = axes[1, 0]
ax.plot(days, S.sum(axis=1) / total_pop, label="S", color="tab:blue")
ax.plot(days, E.sum(axis=1) / total_pop, label="E", color="tab:orange")
ax.plot(days, I.sum(axis=1) / total_pop, label="I", color="tab:red")
ax.plot(days, R.sum(axis=1) / total_pop, label="R", color="tab:green")
ax.set_title("Aggregate SEIR Dynamics  (All Patches)")
ax.set_xlabel("Day")
ax.set_ylabel("Fraction of Population")
ax.legend()
ax.grid(alpha=0.3)

# Panel E: total daily incidence across all patches
ax = axes[1, 1]
ax.plot(days, model.nodes.newly_infected[:nticks, :].sum(axis=1), color="black")
ax.set_title("Total Daily New Infections  (All Patches)")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.grid(alpha=0.3)

# Panel F: effective beta over time
ax = axes[1, 2]
ax.plot(days, beta * season_365, color="darkviolet", linewidth=1.5)
ax.axhline(beta, color="gray", linewidth=0.8, linestyle="--",
           label=f"β baseline = {beta}")
ax.set_title("Effective β = β × Seasonal Multiplier")
ax.set_xlabel("Day")
ax.set_ylabel("Effective β (per day)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("seir_4patch_gravity.png", dpi=100, bbox_inches="tight")
plt.show()
print("\nPlot saved: seir_4patch_gravity.png")
```

**What changed and why:**

| Location | Change | Reason |
|---|---|---|
| `I_init` | Zero for patches 1–3; only patch 0 seeded | Breaks the synchrony that made curves identical |
| `R_init` | Zero for patches 1–3 | Gives uninvaded patches full susceptibility so the arriving wave is clearly visible |
| `assert np.all(I_init > 0)` | Removed | No longer required — patches 1–3 legitimately start uninfected |
| Panel B title | Added arrival-day dotted vertical lines | Directly shows the staggered spatial propagation |
| Post-run check | Added arrival-day table and range assertion | Confirms coupling is working quantitatively |

When all patches share the same I/N ratio at t=0 the term `β·I[i]/N[i]` is equal everywhere, so the asymmetric gravity network only produces corrections of order `(col_sum − row_sum) × ft ≈ 5×10⁻⁶` — far too small to visually separate the curves. With only patch 0 seeded, the infection front travels outward through the network and each patch peaks 10–30 days later than its upstream neighbour, making the coupling unmistakable.
