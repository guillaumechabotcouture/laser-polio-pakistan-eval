The workspace filesystem is read-only, so I'll present the complete code here. Save it as `polio_seir_gravity_seasonal.py`.

---

```python
#!/usr/bin/env python3
"""
Polio SEIR Spatial Model — Gravity Migration + Monsoon Seasonality
===================================================================
Extends a basic LASER polio SEIR model with:

  1. Gravity-model migration network
       M[i,j] = k · pop_i^a · pop_j^b / dist[i,j]^c
       a=1, b=0.5, c=1.5

  2. Monsoon seasonal forcing
       peak  Jul–Oct:  1.3× baseline  (cosine centre ~Sept 1, day 243)
       trough Dec–Mar: 0.7× baseline  (cosine centre ~Mar 2,  day 61)
       Profile normalised so mean == 1.0 → no bias on R_eff.

  3. Row-normalised network: no patch exports > 15% of its FOI.

Geography : 10 patches in a straight line, 50 km apart.
Disease   : Wild poliovirus type 1 (R0≈6, inf ~28 d, exp ~3 d).
Simulation: 20 years (10 yr burn-in + 10 yr analysis).

Run:
    python polio_seir_gravity_seasonal.py
Output:
    polio_gravity_seasonal.png
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists

from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.generic.importation import Infect_Random_Agents
from laser.core.migration import gravity, row_normalizer

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ─── 0. CONFIGURATION ────────────────────────────────────────────────────────

N_PATCHES    = 10
SPACING_KM   = 50.0
NTICKS       = 20 * 365       # total simulation (days)
BURNIN       = 10 * 365       # burn-in before analysis (days)

# Polio epidemiological parameters (wild poliovirus type 1)
R0          = 6.0
INF_MEAN    = 28.0            # infectious period mean (days)
INF_SIGMA   = 3.0             # infectious period std dev
EXP_SHAPE   = 3               # exposed: gamma(3, 1) → mean ≈ 3 days
EXP_SCALE   = 1.0
BETA        = R0 / INF_MEAN   # ≈ 0.214 /day

# Gravity network
GRAVITY_K   = 0.01            # mean export fraction per patch
GRAVITY_A   = 1.0             # source population exponent
GRAVITY_B   = 0.5             # destination population exponent  (user-specified)
GRAVITY_C   = 1.5             # distance decay exponent           (user-specified)
MAX_ROW_SUM = 0.15            # hard cap on row sums              (user-specified)

# Seasonal forcing (monsoon)
# s(d) = 1 + A·cos(2π(d − PEAK_DAY)/365)  — mean 1.0, range [0.7, 1.3]
PEAK_DAY    = 243             # September 1 (centre of Jul–Oct monsoon)
AMPLITUDE   = 0.30            # (1.3 − 0.7) / 2

# Vital dynamics
CBR         = 30.0            # crude birth rate per-1,000/year
CDR         = 7.0             # crude death rate per-1,000/year

# Initial conditions
INIT_IMMUNE_FRAC = 0.50
INIT_SEED        = 3          # infectious agents per patch

# Importation (prevents stochastic extinction in small patches)
IMPORT_PERIOD = 60            # ticks between events
IMPORT_COUNT  = 1             # agents infected per event


# ─── 1. GEOGRAPHIC SCENARIO ──────────────────────────────────────────────────

populations = np.array([
    200_000,   # Patch 00 — small
    500_000,   # Patch 01 — medium
    150_000,   # Patch 02 — small
    800_000,   # Patch 03 — large hub
    300_000,   # Patch 04 — medium
    100_000,   # Patch 05 — smallest
    600_000,   # Patch 06 — large
    250_000,   # Patch 07 — medium
    400_000,   # Patch 08 — medium
    350_000,   # Patch 09 — medium
], dtype=np.int64)
assert len(populations) == N_PATCHES

# Linear arrangement along latitude starting near south Pakistan (30°N, 70°E)
BASE_LAT, BASE_LON = 30.0, 70.0
KM_PER_DEG = 111.0
lats = BASE_LAT + np.arange(N_PATCHES) * (SPACING_KM / KM_PER_DEG)
lons = np.full(N_PATCHES, BASE_LON)

scenario = gpd.GeoDataFrame(
    [
        {
            "nodeid":     i,
            "name":       f"District_{i:02d}",
            "population": int(populations[i]),
            "lat":        float(lats[i]),
            "lon":        float(lons[i]),
            "geometry":   Point(lons[i], lats[i]),
        }
        for i in range(N_PATCHES)
    ],
    crs="EPSG:4326",
)

# Initial compartments — strict: S + E + I + R = population
scenario["R"] = np.minimum(
    np.round(INIT_IMMUNE_FRAC * scenario["population"]).astype(int),
    scenario["population"] - INIT_SEED,
)
scenario["I"] = INIT_SEED
scenario["E"] = 0
scenario["S"] = scenario["population"] - scenario["R"] - scenario["I"] - scenario["E"]

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Compartment sum mismatch"
assert (scenario.S >= 0).all(), "Negative S in initial conditions"
assert scenario.I.sum() > 0,   "No initial infectious agents"

print("─" * 64)
print(f"{'Polio SEIR  ·  Gravity + Monsoon  ·  Linear Corridor':^64}")
print("─" * 64)
print(f"  Patches      : {N_PATCHES} × {int(SPACING_KM)} km")
print(f"  Total pop    : {scenario.population.sum():,}")
print(f"  Init S/I/R   : {scenario.S.sum():,} / {scenario.I.sum()} / {scenario.R.sum():,}")


# ─── 2. SEASONAL FORCING PROFILE ─────────────────────────────────────────────
#
# s(d) = 1 + 0.30·cos(2π(d − 243)/365)
#
# Extrema (before normalisation):
#   d = 243  →  1.30  (Sept 1  — centre of Jul–Oct monsoon)
#   d =  61  →  0.70  (Mar 2   — centre of Dec–Mar dry trough)
#
# The cosine integrates to exactly zero over [0,365), so mean = 1.0.
# We still renormalise explicitly to eliminate floating-point drift.

days_of_year = np.arange(365)
season_365   = 1.0 + AMPLITUDE * np.cos(2.0 * np.pi * (days_of_year - PEAK_DAY) / 365.0)
season_365  /= season_365.mean()   # enforce mean == 1.0

# Validation
assert abs(season_365.mean() - 1.0) < 1e-6, \
    f"Season mean {season_365.mean():.8f} ≠ 1.0"
assert season_365.max() <= 1.31 and season_365.min() >= 0.69, \
    "Seasonal extrema outside [0.7, 1.3] — check AMPLITUDE"

_ref       = pd.Timestamp("2001-01-01")
_peak_day  = int(season_365.argmax())
_trough    = int(season_365.argmin())
print(f"\n  Seasonal forcing  (cosine, amplitude ±{AMPLITUDE}):")
print(f"    Peak   {season_365.max():.3f}×  day {_peak_day:3d}"
      f"  ({(_ref + pd.Timedelta(days=_peak_day)):%d %b})")
print(f"    Trough {season_365.min():.3f}×  day {_trough:3d}"
      f"  ({(_ref + pd.Timedelta(days=_trough)):%d %b})")

# Tile over full simulation and wrap as a ValuesMap (same series for every node)
nnodes       = len(scenario)
season_tiled = np.tile(season_365, NTICKS // 365 + 1)[:NTICKS]
seasonality  = ValuesMap.from_timeseries(season_tiled, nnodes)


# ─── 3. SIMULATION PARAMETERS ────────────────────────────────────────────────

# CRITICAL: CBR/CDR must be per-1,000/year.  Daily per-capita rates cause
# calc_capacity to see near-zero growth → no births → model silently dies out.
assert 1.0 <= CBR <= 60.0, f"CBR={CBR} must be per-1,000/year"
assert 1.0 <= CDR <= 60.0, f"CDR={CDR} must be per-1,000/year"

parameters = PropertySet({
    "prng_seed":  42,
    "nticks":     NTICKS,
    "exp_shape":  EXP_SHAPE,
    "exp_scale":  EXP_SCALE,
    "inf_mean":   INF_MEAN,
    "inf_sigma":  INF_SIGMA,
    "beta":       BETA,
    "cbr":        CBR,
    # Infect_Random_Agents reads these four keys from model.params:
    "importation_period": IMPORT_PERIOD,
    "importation_count":  IMPORT_COUNT,
    "importation_start":  0,
    "importation_end":    NTICKS,
    # Headroom for 20 years of CBR > CDR growth
    "capacity_safety_factor": 2.5,
})

print(f"\n  Epi: R0={R0}  beta={BETA:.4f}/d  "
      f"inf={INF_MEAN}±{INF_SIGMA}d  exp~{EXP_SHAPE * EXP_SCALE:.0f}d")


# ─── 4. DISTRIBUTIONS ────────────────────────────────────────────────────────

expdurdist = dists.gamma(shape=EXP_SHAPE, scale=EXP_SCALE)
infdurdist = dists.normal(loc=INF_MEAN,   scale=INF_SIGMA)


# ─── 5. VITAL DYNAMICS ARRAYS ────────────────────────────────────────────────

# Shape (NTICKS+1, nnodes); values in per-1,000/year.
# The extra row ensures calc_capacity (called in Model.__init__) has full
# coverage for the entire simulation horizon.
birthrate_values = np.full((NTICKS + 1, nnodes), CBR, dtype=np.float32)
deathrate_values = np.full((NTICKS + 1, nnodes), CDR, dtype=np.float32)

# Stable age pyramid for initialising agent ages at t = 0.
# P(age) ∝ exp(−μ · age_years) where μ = CDR/1000/year.
_daily_mort  = CDR / (1_000.0 * 365.0)
_age_weights = np.exp(-_daily_mort * 365.0 * np.arange(80)).astype(np.float64)
pyramid = AliasedDistribution(_age_weights)


# ─── 6. BUILD MODEL ──────────────────────────────────────────────────────────

print(f"\n  Building model…")
model = Model(scenario, parameters, birthrates=birthrate_values)

# NOTE: gravity_k/a/b/c are intentionally absent from params so that
# Model.__init__() does NOT auto-compute the network.  We set model.network
# manually below to apply the custom 15% row-sum cap.


# ─── 7. GRAVITY MIGRATION NETWORK ────────────────────────────────────────────
#
# M[i,j] = k · pop_i^a · pop_j^b / dist[i,j]^c
#   a = 1.0  (source population exponent)
#   b = 0.5  (destination population exponent — user-specified)
#   c = 1.5  (distance decay exponent         — user-specified)
#
# Steps:
#   1. Exact inter-patch distances (km) — linear layout makes this trivial.
#   2. Raw gravity matrix with k=1 (unnormalised).
#   3. Scale so that GRAVITY_K == mean fraction exported per patch.
#   4. Cap every row at MAX_ROW_SUM = 0.15 with row_normalizer().
#
# row_normalizer(M, cap) rescales rows whose sum > cap to exactly cap;
# rows already below cap are left unchanged.

# Step 1 — exact distances
patch_idx = np.arange(N_PATCHES)
dist_km   = np.abs(patch_idx[:, None] - patch_idx[None, :]).astype(float) * SPACING_KM
np.fill_diagonal(dist_km, np.inf)   # infinity → no self-migration term

# Step 2 — raw gravity (k=1)
raw_net = gravity(
    np.array(scenario.population, dtype=float),
    dist_km,
    1,           # k  (we rescale below)
    GRAVITY_A,   # a = 1.0
    GRAVITY_B,   # b = 0.5
    GRAVITY_C,   # c = 1.5
)

# Step 3 — scale to desired mean export fraction
mean_export = raw_net.sum(axis=1).mean()
assert mean_export > 0, "Raw gravity matrix is all-zero; check dist_km"
network_scaled = raw_net * (GRAVITY_K / mean_export)

# Step 4 — enforce 15% cap
model.network = row_normalizer(network_scaled, MAX_ROW_SUM)

# Validation
row_sums = model.network.sum(axis=1)
assert model.network.sum() > 0,              "Network is all-zero — patches isolated"
assert row_sums.max() <= MAX_ROW_SUM + 1e-9, f"Row sum cap exceeded: {row_sums.max():.6f}"
assert (model.network >= 0).all(),           "Negative network weights"
assert np.allclose(model.network.diagonal(), 0.0), "Non-zero diagonal (self-migration)"

print(f"\n  Gravity network (a={GRAVITY_A}, b={GRAVITY_B}, c={GRAVITY_C}, cap={MAX_ROW_SUM:.0%}):")
print(f"    Row sums — min:{row_sums.min():.4f}  "
      f"mean:{row_sums.mean():.4f}  max:{row_sums.max():.4f}")
with np.printoptions(precision=2, suppress=True, linewidth=90):
    print("    M[i→j] in % FOI:")
    print((model.network * 100).__repr__())


# ─── 8. ASSEMBLE COMPONENTS ──────────────────────────────────────────────────
#
# Execution order per tick:
#   Susceptible    — propagate S[t+1]=S[t]; anchors S+E+I+R=N invariant
#   Exposed        — decrement etimer; E→I + sample itimer
#   Infectious     — decrement itimer; I→R transitions
#   Recovered      — propagate R[t+1]=R[t]; permanent immunity
#   BirthsByCBR    — add newborns; fires on_birth callbacks
#   MortalityByCDR — remove deaths; decrement compartment counts
#   Importation    — periodic seeding (prevents extinction in small patches)
#   Transmission   — compute seasonal FOI with spatial network coupling
#                    (must run LAST so I[t] is fully updated)

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model,    expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    BirthsByCBR(model,    birthrates=birthrate_values, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=deathrate_values),
    Infect_Random_Agents(model),   # reads importation_* keys from model.params
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
]

print(f"\n  Components ({len(model.components)}):")
for c in model.components:
    print(f"    {type(c).__name__}")


# ─── 9. RUN ──────────────────────────────────────────────────────────────────

print(f"\n  Running {NTICKS // 365}-year simulation…")
model.run("Polio SEIR — Gravity + Monsoon Seasonality")
print("  Complete.\n")


# ─── 10. POST-RUN VERIFICATION ───────────────────────────────────────────────

newly = model.nodes.newly_infected[:NTICKS]     # (NTICKS, N_PATCHES)
S_arr = model.nodes.S[:NTICKS]
E_arr = model.nodes.E[:NTICKS]
I_arr = model.nodes.I[:NTICKS]
R_arr = model.nodes.R[:NTICKS]
N_arr = S_arr + E_arr + I_arr + R_arr

print("Verification:")

# 1. Population growth (CBR=30 > CDR=7 → population should increase)
pop_chg = (N_arr[-1].sum() - N_arr[0].sum()) / N_arr[0].sum() * 100.0
print(f"  Population:   {N_arr[0].sum():,} → {N_arr[-1].sum():,}  ({pop_chg:+.1f}%)")
assert 0.0 < pop_chg < 200.0, \
    f"Population change {pop_chg:.1f}% unrealistic — check CBR/CDR units"

# 2. No negative compartment counts
for label, arr in [("S", S_arr), ("E", E_arr), ("I", I_arr), ("R", R_arr)]:
    assert (arr >= 0).all(), f"Negative {label} detected — depletion bug"
print("  Compartments: all non-negative ✓")

# 3. Disease persists in post-burn-in window
post_newly  = newly[BURNIN:]
total_cases = int(post_newly.sum())
assert total_cases > 0, "Disease extinct in post-burn-in — check beta / importation"
print(f"  Post-burn-in cases: {total_cases:,}")

# 4. Spatial coupling active
active_patches = int((post_newly.sum(axis=0) > 0).sum())
print(f"  Active patches: {active_patches}/{N_PATCHES}")
assert active_patches > 1, "Only 1 active patch — network may be all-zero"

# 5. Seasonal signal: monsoon peak (Jul–Oct) > dry trough (Dec–Mar)
ticks_post  = np.arange(BURNIN, NTICKS)
doy         = ticks_post % 365
daily_inc   = post_newly.sum(axis=1)

peak_mask   = (doy >= 182) & (doy <= 304)    # Jul 1 – Oct 31
trough_mask = (doy <= 89)  | (doy >= 335)    # Dec 1 – Mar 31

peak_mean   = float(daily_inc[peak_mask].mean())
trough_mean = float(daily_inc[trough_mask].mean())
print(f"  Seasonal ratio: {peak_mean/max(trough_mean,1e-9):.2f}×  "
      f"(peak {peak_mean:.1f}/d, trough {trough_mean:.1f}/d)")
assert peak_mean > trough_mean, \
    "Seasonal forcing not visible — peak incidence ≤ trough incidence"

print("  All checks passed ✓\n")


# ─── 11. DIAGNOSTIC PLOTS ────────────────────────────────────────────────────

years      = np.arange(NTICKS) / 365.0
post_years = np.arange(BURNIN, NTICKS) / 365.0

fig = plt.figure(figsize=(14, 10))
fig.suptitle(
    f"Polio SEIR · {N_PATCHES} patches @ {int(SPACING_KM)} km\n"
    f"Gravity (a={GRAVITY_A}, b={GRAVITY_B}, c={GRAVITY_C}, cap {MAX_ROW_SUM:.0%}) · "
    f"Monsoon forcing (1.3× peak, 0.7× trough)",
    fontsize=10.5,
)
gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.48, wspace=0.35)

# ── Panel A: Seasonal forcing profile ────────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])
ax_a.plot(days_of_year, season_365, color="steelblue", lw=1.8)
ax_a.axhline(1.0, color="gray", lw=0.8, ls="--", label="Baseline (1.0×)")
ax_a.axhspan(182, 304, alpha=0.15, color="coral",   label="Jul–Oct  peak (1.3×)")
ax_a.axhspan(  0,  90, alpha=0.12, color="skyblue", label="Dec–Mar trough (0.7×)")
ax_a.axhspan(335, 365, alpha=0.12, color="skyblue")
ax_a.axvline(PEAK_DAY, color="coral",   lw=0.9, ls=":", alpha=0.7)
ax_a.axvline(_trough,  color="skyblue", lw=0.9, ls=":", alpha=0.7)
ax_a.set_xlim(0, 364); ax_a.set_ylim(0.55, 1.50)
ax_a.set_xlabel("Day of year"); ax_a.set_ylabel("Transmission multiplier")
ax_a.set_title("Seasonal Forcing Profile")
ax_a.legend(fontsize=7)

# ── Panel B: Gravity network heatmap ─────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
im_b = ax_b.imshow(model.network * 100.0, cmap="YlOrRd", aspect="auto", vmin=0)
plt.colorbar(im_b, ax=ax_b, label="FOI exported (%)")
ax_b.set_title(f"Gravity Network  (cap = {MAX_ROW_SUM:.0%})")
ax_b.set_xlabel("Destination patch"); ax_b.set_ylabel("Source patch")
ax_b.set_xticks(range(N_PATCHES)); ax_b.set_yticks(range(N_PATCHES))

# ── Panel C: National epidemic curve (post-burn-in) ──────────────────────────
ax_c = fig.add_subplot(gs[1, :])
raw_inc = post_newly.sum(axis=1).astype(float)
smooth  = np.convolve(raw_inc, np.ones(14) / 14.0, mode="same")
ax_c.fill_between(post_years, raw_inc, alpha=0.20, color="firebrick")
ax_c.plot(post_years, smooth, color="firebrick", lw=1.2, label="14-day rolling avg")
for yr in range(BURNIN // 365, NTICKS // 365):
    ax_c.axvspan(yr + 182/365, yr + 304/365, alpha=0.07, color="coral", lw=0)
ax_c.axhline(raw_inc.mean(), color="gray", ls="--", lw=0.8,
             label=f"Mean = {raw_inc.mean():.1f}/d")
ax_c.set_xlabel("Year"); ax_c.set_ylabel("New infections / day")
ax_c.set_title("National Epidemic Curve (post-burn-in; coral bands = monsoon)")
ax_c.legend(fontsize=8)

# ── Panel D: Spatial incidence heatmap ───────────────────────────────────────
ax_d = fig.add_subplot(gs[2, 0])
n_wk   = (NTICKS - BURNIN) // 7
weekly = post_newly[:n_wk * 7].reshape(n_wk, 7, N_PATCHES).sum(axis=1).T
im_d   = ax_d.imshow(
    weekly, aspect="auto", cmap="hot",
    extent=[BURNIN/365, BURNIN/365 + n_wk/52.18, N_PATCHES - 0.5, -0.5],
)
plt.colorbar(im_d, ax=ax_d, label="Weekly cases")
ax_d.set_xlabel("Year"); ax_d.set_ylabel("Patch  (50 km apart)")
ax_d.set_title("Spatial Incidence Heatmap")

# ── Panel E: Susceptible fraction + network row sums ─────────────────────────
ax_e = fig.add_subplot(gs[2, 1])
frac_S = S_arr.sum(axis=1) / N_arr.sum(axis=1)
ax_e.plot(years, frac_S, color="darkgreen", lw=1.2, label="S/N (national)")
ax_e.axvline(BURNIN/365, color="gray", ls="--", lw=0.8, label="End burn-in")
ax_e.set_xlabel("Year  /  Patch index")
ax_e.set_ylabel("Susceptible fraction", color="darkgreen")
ax_e.tick_params(axis="y", labelcolor="darkgreen")

ax_e2 = ax_e.twinx()
ax_e2.bar(range(N_PATCHES), row_sums * 100, color="steelblue",
          alpha=0.45, width=0.6, label="Row sum (%)")
ax_e2.axhline(MAX_ROW_SUM * 100, color="steelblue", ls=":", lw=1.0,
              label=f"Cap = {MAX_ROW_SUM:.0%}")
ax_e2.set_ylabel("Network row sum (%)", color="steelblue")
ax_e2.tick_params(axis="y", labelcolor="steelblue")
ax_e.set_title("Susceptible Fraction & Network Row Sums")

h1, l1 = ax_e.get_legend_handles_labels()
h2, l2 = ax_e2.get_legend_handles_labels()
ax_e.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")

plt.savefig("polio_gravity_seasonal.png", dpi=130, bbox_inches="tight")
print("Diagnostic figure → polio_gravity_seasonal.png")
plt.show()
```

---

Here is what each section does and why it is built that way:

**Gravity network (Section 7)**

```
M[i,j] = k · pop_i^a · pop_j^b / dist_ij^c    a=1, b=0.5, c=1.5
```

`gravity()` is called with `k=1` to get unnormalized weights, which are then scaled so that `GRAVITY_K` equals the *mean* row sum across all patches. `row_normalizer(M, 0.15)` then clips any row exceeding 0.15 down to exactly 0.15, leaving under-cap rows untouched. Patches that are close to a large hub (Patch 03 pop 800k) receive a disproportionately large destination-pull due to `pop^0.5`, while the `dist^1.5` decay makes the network sparse beyond a few hops.

**Seasonal forcing (Section 2)**

```
s(d) = 1 + 0.30·cos(2π(d − 243)/365)
```

A cosine with amplitude 0.3 centred on day 243 (Sept 1) gives `s_max = 1.30` at the peak and `s_min = 0.70` six months later (day 61, early March). Both extrema fall inside the requested windows. The profile is renormalized to `mean = 1.0` so it does not bias R_eff — an un-normalized profile is the most common silent error with seasonal forcing.

**Key correctness guards**

| What | Where | Why |
|---|---|---|
| `assert 1 ≤ CBR/CDR ≤ 60` | Section 3 | Catches daily per-capita rate (silent birth failure) |
| `assert season.mean() == 1.0` | Section 2 | Prevents systematic R_eff bias |
| `assert row_sums.max() ≤ 0.15` | Section 7 | Confirms cap is enforced |
| Post-run seasonal ratio check | Section 10 | Confirms forcing is epidemiologically visible |
| `SEIR.Transmission` runs last | Section 8 | Uses fully-updated I[t] for FOI |
