"""
Gravity-Model 4-Patch SEIR — LASER Framework  (Seasonal Forcing + Vital Dynamics)
===================================================================================
Respiratory disease: R0=5, latent period=4d, infectious period=10d
Patches:  100k / 200k / 150k / 80k (east-west line, 75 km apart)
Initial:  Patch 0: 1% I, 99% S  |  Patches 1-3: 100% S (fully naive)
Duration: 10 years (3 650 days)
Gravity coupling: k=0.01, a=1, b=1, c=1.5 | Row-normalised ≤ 15% export
Seasonality: winter peak (days 0-90) 1.3×, summer trough (days 150-240) 0.7×
             tiled annually for the full 10-year run
Vital dynamics: CBR=30 per 1000/year, CDR=10 per 1000/year → ~2% net growth
                Expected after 10 years: 530 000 × 1.02^10 ≈ 646 000
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
import matplotlib.pyplot as plt

# ── 1. Patch populations and initial conditions ───────────────────────────────
populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
n_patches   = len(populations)

I_init    = np.zeros(n_patches, dtype=np.int32)
I_init[0] = max(1, int(0.01 * populations[0]))   # 1% of Patch 0 infected (1 000)
R_init    = np.zeros(n_patches, dtype=np.int32)
E_init    = np.zeros(n_patches, dtype=np.int32)
S_init    = populations - E_init - I_init - R_init

assert np.all(S_init + E_init + I_init + R_init == populations), \
    "Initial S+E+I+R must equal population in every patch"
assert I_init[0] > 0,           "Patch 0 must have at least one initial infectious agent"
assert np.all(I_init[1:] == 0), "Only Patch 0 should be seeded with infection"

# ── 2. Gravity-model migration network ───────────────────────────────────────
GRAVITY_K  = 0.01
GRAVITY_A  = 1.0
GRAVITY_B  = 1.0
GRAVITY_C  = 1.5
MAX_EXPORT = 0.15

pos_km    = np.array([0.0, 75.0, 150.0, 225.0])
dist_km   = np.abs(pos_km[:, None] - pos_km[None, :])
pop_f     = populations.astype(np.float64)
dist_safe = dist_km + np.eye(n_patches)   # diagonal = 1.0 to avoid /0

network = (
    GRAVITY_K
    * (pop_f[:, None] ** GRAVITY_A)
    * (pop_f[None, :] ** GRAVITY_B)
    / (dist_safe ** GRAVITY_C)
)
np.fill_diagonal(network, 0.0)
network /= pop_f[:, None]   # convert to per-capita daily rates

row_sums = network.sum(axis=1)
scale    = np.where(row_sums > MAX_EXPORT, MAX_EXPORT / row_sums, 1.0)
network *= scale[:, None]

assert np.all(network >= 0),                             "Negative migration rate"
assert np.all(network.sum(axis=1) <= MAX_EXPORT + 1e-9), "Row sum exceeds cap"
np.testing.assert_array_almost_equal(np.diag(network), 0, err_msg="Diagonal non-zero")

print("Gravity migration network (fraction of population per day):")
print("         " + "  ".join(f"   Patch {j}" for j in range(n_patches)))
for i in range(n_patches):
    vals = "  ".join(f"{network[i, j]:.5f}" for j in range(n_patches))
    print(f"Patch {i}:  {vals}")
print(f"Row sums: {network.sum(axis=1).round(4)}")

# ── 3. Scenario GeoDataFrame ──────────────────────────────────────────────────
LON_PER_KM = 1.0 / (111.32 * np.cos(np.radians(33.0)))

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(n_patches),
        "name":       [f"Patch {i}" for i in range(n_patches)],
        "population": populations,
        "S":          S_init,
        "E":          E_init,
        "I":          I_init,
        "R":          R_init,
        "geometry":   [Point(72.0 + i * 75.0 * LON_PER_KM, 33.0)
                       for i in range(n_patches)],
    },
    crs="EPSG:4326",
)

print("\nInitial conditions:")
for _, row in scenario.iterrows():
    N = row["population"]
    print(
        f"  {row['name']}: N={N:,}  "
        f"S={row['S']:,} ({row['S']/N*100:.1f}%)  "
        f"I={row['I']:,} ({row['I']/N*100:.1f}%)  "
        f"R={row['R']:,} ({row['R']/N*100:.1f}%)"
    )

# ── 4. Vital-dynamics rates ────────────────────────────────────────────────────
# BirthsByCBR and MortalityByCDR both expect per-1000/year (NOT daily per-capita).
# Passing daily per-capita (e.g. 0.00008) is the #1 silent failure: calc_capacity
# sees near-zero growth, pre-allocates ≈ initial pop, and LaserFrame.add() has
# no free slots, so no births occur at all.
CBR = 30.0   # crude birth rate  (per 1 000 / year) — typical Sub-Saharan Africa
CDR = 10.0   # crude death rate  (per 1 000 / year)
# Net growth = (30 - 10) / 1000 = 2% / year
# After 10 years: 530 000 × 1.02^10 ≈ 646 000

assert 1 <= CBR <= 60, f"CBR must be per-1000/year; got {CBR}"
assert 1 <= CDR <= 60, f"CDR must be per-1000/year; got {CDR}"
assert CBR > CDR,      "Expected net growth (CBR > CDR); check values"

nticks = 10 * 365   # 3 650 ticks

# Constant-rate (nticks × n_patches) arrays — uniform across patches and time
birthrate_array = np.full((nticks, n_patches), CBR)   # per-1000/year
deathrate_array = np.full((nticks, n_patches), CDR)   # per-1000/year

# ── 5. Model parameters ───────────────────────────────────────────────────────
# capacity_safety_factor: Model.__init__ passes this to calc_capacity to
# pre-allocate agent slots.  10-year run at ~2% growth → max factor ≈ 1.22,
# so 2.0× is comfortable; increase to 3–4 for faster growth or longer runs.
params = PropertySet({
    "prng_seed":              42,
    "nticks":                 nticks,
    "beta":                   0.5,        # R0 / infectious_period = 5 / 10
    "capacity_safety_factor": 2.0,
})

# ── 5b. Seasonal forcing — tiled for 10 years ─────────────────────────────────
WINTER_MULT = 1.3
SUMMER_MULT = 0.7

days_365 = np.arange(365, dtype=np.float64)
season_365 = np.where(
    days_365 < 90,
    WINTER_MULT,
    np.where(
        days_365 < 150,
        WINTER_MULT + (SUMMER_MULT - WINTER_MULT) * (days_365 - 90.0) / 60.0,
        np.where(
            days_365 < 240,
            SUMMER_MULT,
            SUMMER_MULT + (WINTER_MULT - SUMMER_MULT) * (days_365 - 240.0) / 125.0,
        ),
    ),
)
season_365 /= season_365.mean()   # normalise → mean = 1.0 exactly

assert abs(season_365.mean() - 1.0) < 1e-10, \
    f"Seasonal mean = {season_365.mean():.8f}, expected 1.0"
assert np.all(season_365 >= 0.0), "Seasonal multiplier contains negative values"

# Tile the 365-day profile across 10 years
season     = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season, n_patches)

print(f"\nSeasonal profile tiled for {nticks // 365} years "
      f"(normalised mean = {season.mean():.6f})")
print(f"  Winter (days 0-90):    {season_365[:90].mean():.3f}× beta  "
      f"(eff. β = {params.beta * season_365[:90].mean():.4f})")
print(f"  Summer (days 150-240): {season_365[150:240].mean():.3f}× beta  "
      f"(eff. β = {params.beta * season_365[150:240].mean():.4f})")

# ── 6. Age pyramid for newborns ────────────────────────────────────────────────
# AliasedDistribution with a single bin at index 0 so every newborn is assigned
# age 0 (i.e. dob = current tick).  The parenthetical in the API —
# "Age pyramid for sampling newborn ages (usually age 0)" — confirms this use.
pyramid = AliasedDistribution(np.array([1.0]))

# ── 7. Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.constant_int(4)    # latent period:     exactly 4 days
infdurdist = dists.constant_int(10)   # infectious period: exactly 10 days

# ── 8. Construct model ────────────────────────────────────────────────────────
# Passing birthrates= to Model.__init__() causes it to call calc_capacity()
# internally, pre-allocating enough agent slots for the projected growth ×
# capacity_safety_factor.  If birthrates were wrong units (daily per-capita)
# calc_capacity would see ~0 growth → capacity ≈ initial pop → no births.
model = Model(scenario, params, birthrates=birthrate_array)
model.network = network   # assign gravity coupling matrix

# ── 9. Component list ─────────────────────────────────────────────────────────
# BirthsByCBR / MortalityByCDR are placed AFTER the epidemiological components
# so that transmission and state transitions happen first; deaths in the same
# tick do not prematurely remove just-exposed agents.
#
# IMPORTANT: MortalityByCDR uses keyword `mortalityrates=`, NOT `deathrates=`.
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model,    birthrates=birthrate_array,    pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=deathrate_array),
]

# ── 10. Record initial total population ──────────────────────────────────────
N_start = int(populations.sum())
print(f"\nTotal population at START (day 0): {N_start:,}")

# ── 11. Run ────────────────────────────────────────────────────────────────────
print(f"Running 4-patch SEIR + vital dynamics for {nticks} days (10 years)...")
model.run("GravitySEIR_4patch_VitalDynamics")
print("Done.\n")

# ── 12. Total population at end ───────────────────────────────────────────────
# LASER stores compartment arrays at indices 0 … nticks-1; final state = nticks-1.
last = nticks - 1
N_end_patches = (
    model.nodes.S[last, :] + model.nodes.E[last, :] +
    model.nodes.I[last, :] + model.nodes.R[last, :]
)
N_end = int(N_end_patches.sum())

# Naive compound-interest prediction: N(t) = N0 × (1 + r)^t,  r = (CBR-CDR)/1000
expected_factor = (1.0 + (CBR - CDR) / 1000.0) ** 10
N_expected      = int(round(N_start * expected_factor))

print("=== Population Summary ===")
print(f"  Total at START (day    0): {N_start:>10,}")
print(f"  Total at END   (day 3650): {N_end:>10,}")
print(f"  Expected (~2%/yr naive):   {N_expected:>10,}")
print(f"  Simulated growth factor:   {N_end / N_start:.4f}×  over 10 years")
print(f"  Equivalent annual rate:    {(N_end / N_start) ** 0.1 - 1:.3%}/year")
print()
print("  Per-patch breakdown:")
for idx in range(n_patches):
    n0 = int(populations[idx])
    ne = int(N_end_patches[idx])
    print(f"    Patch {idx}: {n0:>8,} → {ne:>8,}  (×{ne/n0:.3f})")

# Sanity: CBR > CDR → population must have grown
assert N_end > N_start, \
    f"Population should have grown (CBR={CBR} > CDR={CDR}) but {N_end} < {N_start}"

# Allow ±15% of the naive prediction to account for stochastic variation and the
# fact that epidemic deaths are slightly additive to background mortality
growth_deviation = abs(N_end / N_start - expected_factor) / expected_factor
assert growth_deviation < 0.15, (
    f"Simulated growth factor {N_end/N_start:.4f} deviates "
    f"{growth_deviation:.1%} from expected {expected_factor:.4f}"
)
print(f"\nPopulation growth sanity check passed: "
      f"{(N_end/N_start)**0.1 - 1:.2%}/yr matches expected ~2%/yr")

# ── 13. Epidemic summary ──────────────────────────────────────────────────────
print("\n=== Epidemic Summary (full 10 years) ===")
for idx in range(n_patches):
    I_trace   = model.nodes.I[:nticks, idx]
    peak_I    = int(I_trace.max())
    peak_day  = int(I_trace.argmax())
    total_inf = int(model.nodes.newly_infected[:nticks, idx].sum())
    print(
        f"  Patch {idx} (init N={int(populations[idx]):,}): "
        f"peak I = {peak_I:,} on day {peak_day} (year {peak_day // 365 + 1}) | "
        f"total infections (10 yr) = {total_inf:,}"
    )

for name, arr in [("S", model.nodes.S), ("E", model.nodes.E),
                  ("I", model.nodes.I), ("R", model.nodes.R)]:
    assert np.all(arr[:nticks, :] >= 0), f"Negative values in {name} compartment"
print("\nCompartment non-negativity check passed")

# Verify traveling-wave property is preserved (Patch 0 peaks first)
peak_days = [int(model.nodes.I[:nticks, idx].argmax()) for idx in range(n_patches)]
print(f"First-wave peak days by patch: {peak_days}")
assert peak_days[0] <= peak_days[-1], \
    "Expected Patch 0 (source) to peak no later than Patch 3 (distal)"
print("Spatial coupling check passed: infection wave propagates west → east")

# ── 14. Figures ───────────────────────────────────────────────────────────────
# Layout: 3 rows × 2 cols
#   (1,1) migration heatmap          |  (1,2) total population trajectory
#   (2,1) Patch 0 SEIR (10 years)    |  (2,2) Patch 1 SEIR
#   (3,1) Patch 2 SEIR               |  (3,2) Patch 3 SEIR

# Per-tick population per patch  (S+E+I+R — grows with vital dynamics)
N_ticks = (model.nodes.S[:nticks, :] + model.nodes.E[:nticks, :] +
           model.nodes.I[:nticks, :] + model.nodes.R[:nticks, :])   # (nticks, 4)
N_total = N_ticks.sum(axis=1)

PATCH_COLORS = ["steelblue", "darkorange", "forestgreen", "crimson"]
t            = np.arange(nticks)
year_ticks   = np.arange(0, nticks + 1, 365)
year_labels  = [str(y) for y in range(11)]   # '0' … '10'

fig = plt.figure(figsize=(15, 12))
fig.suptitle(
    "4-Patch SEIR — Gravity Coupling + Seasonal Forcing + Vital Dynamics (10 years)\n"
    r"R$_0$=5, Latent=4d, Inf=10d  |  k=0.01, a=b=1, c=1.5  |  "
    "CBR=30, CDR=10/1000/yr  |  Winter 1.3×, Summer 0.7×  |  Seed: Patch 0 only",
    fontsize=11,
)

# Panel 1 — migration network heatmap
ax_net = fig.add_subplot(3, 2, 1)
im = ax_net.imshow(network * 100, cmap="YlOrRd", aspect="auto", vmin=0)
plt.colorbar(im, ax=ax_net, label="% of pop / day")
tick_labels = [f"P{j}\n({populations[j]//1000}k)" for j in range(n_patches)]
ax_net.set_xticks(range(n_patches)); ax_net.set_yticks(range(n_patches))
ax_net.set_xticklabels(tick_labels, fontsize=8)
ax_net.set_yticklabels([f"P{i}" for i in range(n_patches)], fontsize=8)
ax_net.set_xlabel("Destination", fontsize=8); ax_net.set_ylabel("Origin", fontsize=8)
ax_net.set_title("Daily migration rates\n(row = origin patch)", fontsize=9)
for i in range(n_patches):
    for j in range(n_patches):
        txt_color = "white" if network[i, j] > 0.09 else "black"
        ax_net.text(j, i, f"{network[i, j]*100:.2f}%",
                    ha="center", va="center", fontsize=7, color=txt_color)

# Panel 2 — total population trajectory
ax_pop = fig.add_subplot(3, 2, 2)
N_exp_curve = N_start * (1.0 + (CBR - CDR) / 1000.0) ** (t / 365.0)
ax_pop.plot(t, N_total / 1_000, color="black",  lw=2.0, label="Simulated total")
ax_pop.plot(t, N_exp_curve / 1_000, color="grey", lw=1.5, ls="--",
            label=f"Expected compound ({(CBR-CDR)/10:.0f}%/yr)")
for idx in range(n_patches):
    ax_pop.plot(t, N_ticks[:, idx] / 1_000, color=PATCH_COLORS[idx],
                lw=1.0, alpha=0.65, label=f"Patch {idx}")
ax_pop.set_xticks(year_ticks); ax_pop.set_xticklabels(year_labels, fontsize=8)
ax_pop.set_xlim(0, nticks - 1)
ax_pop.set_xlabel("Year"); ax_pop.set_ylabel("Population (thousands)")
ax_pop.set_title(
    f"Population trajectory  (CBR={CBR:.0f}, CDR={CDR:.0f} per 1000/yr)\n"
    f"Start: {N_start:,}  →  End: {N_end:,}  "
    f"(×{N_end/N_start:.3f} in 10 yr  ≈  {(N_end/N_start)**0.1-1:.2%}/yr)",
    fontsize=9,
)
ax_pop.legend(fontsize=7, loc="upper left", ncol=2)
ax_pop.grid(alpha=0.3)

# Panels 3–6 — SEIR curves per patch  (% of current patch population)
# The denominator N grows each year due to vital dynamics; expressing compartments
# as % of current N keeps the y-axis on a 0-100 scale and makes endemic dynamics
# visible even as absolute counts grow.
for idx in range(n_patches):
    ax = fig.add_subplot(3, 2, idx + 3)

    # Shade repeating seasonal windows (winter=blue, summer=red)
    for yr in range(10):
        offset = yr * 365
        ax.axvspan(offset,       offset + 90,  alpha=0.05, color="cornflowerblue")
        ax.axvspan(offset + 150, offset + 240, alpha=0.05, color="tomato")

    N_patch = np.maximum(N_ticks[:, idx], 1)   # avoid /0
    pct     = 100.0 / N_patch

    ax.plot(t, model.nodes.S[:nticks, idx] * pct, label="S", color="steelblue",   lw=1.5)
    ax.plot(t, model.nodes.E[:nticks, idx] * pct, label="E", color="darkorange",  lw=1.5)
    ax.plot(t, model.nodes.I[:nticks, idx] * pct, label="I", color="crimson",     lw=1.5)
    ax.plot(t, model.nodes.R[:nticks, idx] * pct, label="R", color="forestgreen", lw=1.5)

    peak_day = int(model.nodes.I[:nticks, idx].argmax())
    ax.axvline(peak_day, color="crimson", ls="--", lw=1, alpha=0.5)

    seed_label = "  [seeded]" if idx == 0 else ""
    ax.set_title(
        f"Patch {idx}  (init N={int(populations[idx]):,}){seed_label}\n"
        f"Peak I day {peak_day} (yr {peak_day//365+1})  |  "
        f"End N={int(N_end_patches[idx]):,}",
        fontsize=9,
    )
    ax.set_xticks(year_ticks); ax.set_xticklabels(year_labels, fontsize=8)
    ax.set_xlabel("Year"); ax.set_ylabel("% of current population")
    ax.set_xlim(0, nticks - 1); ax.set_ylim(0, 100)
    ax.legend(loc="center right", fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("seir_4patch_gravity_seasonal_vital.png", dpi=150, bbox_inches="tight")
plt.show()
print("Plot saved → seir_4patch_gravity_seasonal_vital.png")