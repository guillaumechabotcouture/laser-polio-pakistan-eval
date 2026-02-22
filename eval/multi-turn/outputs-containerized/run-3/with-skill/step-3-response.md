I have everything I need. The mean of the piecewise-linear profile is exactly 1.0 by construction (both transition ramps average (1.3+0.7)/2 = 1.0), so `beta` retains its meaning as the average transmission rate.

Here is the complete modified script:

```python
"""
Gravity-Model 4-Patch SEIR — LASER Framework  (Seasonal Forcing)
=================================================================
Respiratory disease: R0=5, latent period=4d, infectious period=10d
Patches:  100k / 200k / 150k / 80k (east-west line, 75 km apart)
Initial:  90% S, 0% E, 1% I, 9% R
Duration: 365 days (1 year)
Gravity coupling: k=0.01, a=1, b=1, c=1.5 | Row-normalised ≤ 15% export
Seasonality: winter peak (days 0-90) 1.3×, summer trough (days 150-240) 0.7×
             linear transitions between windows; mean normalised to 1.0
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap          # ← seasonal forcing
import matplotlib.pyplot as plt

# ── 1. Patch populations and initial conditions ───────────────────────────────
populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
n_patches = len(populations)

I_init = np.round(0.01 * populations).astype(np.int32)
R_init = np.round(0.09 * populations).astype(np.int32)
E_init = np.zeros(n_patches, dtype=np.int32)
S_init = populations - E_init - I_init - R_init

assert np.all(S_init + E_init + I_init + R_init == populations), \
    "Initial S+E+I+R must equal population in every patch"
assert np.all(I_init > 0), "Every patch needs at least one initial infectious agent"

# ── 2. Gravity-model migration network ───────────────────────────────────────
# Patches arranged in an east-west line, 75 km apart:
#   Patch 0 ─── 75 km ─── Patch 1 ─── 75 km ─── Patch 2 ─── 75 km ─── Patch 3
#
# Gravity flow:  W[i,j] = k * pop_i^a * pop_j^b / dist[i,j]^c
# Per-capita:    network[i,j] = W[i,j] / pop_i
# Row-normalise: if sum_j(network[i,j]) > MAX_EXPORT, scale row down

GRAVITY_K  = 0.01
GRAVITY_A  = 1.0
GRAVITY_B  = 1.0
GRAVITY_C  = 1.5
MAX_EXPORT = 0.15   # max fraction of patch i population that leaves per tick

# Pairwise distances (km) for a 1-D line layout
pos_km  = np.array([0.0, 75.0, 150.0, 225.0])
dist_km = np.abs(pos_km[:, None] - pos_km[None, :])   # shape (4, 4)

# Gravity weights via broadcasting.
# Add identity to distance matrix so the zero diagonal doesn't cause /0;
# the diagonal is zeroed out immediately after.
pop_f     = populations.astype(np.float64)
dist_safe = dist_km + np.eye(n_patches)                # diagonal = 1.0 (dummy)

network = (
    GRAVITY_K
    * (pop_f[:, None] ** GRAVITY_A)
    * (pop_f[None, :] ** GRAVITY_B)
    / (dist_safe ** GRAVITY_C)
)
np.fill_diagonal(network, 0.0)   # no self-migration

# Convert raw weights to per-capita daily rates  (divide each row i by pop_i)
network /= pop_f[:, None]

# Row-normalise: clamp total daily outflow to MAX_EXPORT per patch
row_sums = network.sum(axis=1)
scale    = np.where(row_sums > MAX_EXPORT, MAX_EXPORT / row_sums, 1.0)
network *= scale[:, None]

# Sanity checks on the network
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
# East-west line: 75 km ≈ 0.803° longitude at 33°N  (for geometry display only;
# distances are computed analytically above in km)
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
        f"S={row['S']:,} ({row['S']/N*100:.0f}%)  "
        f"I={row['I']:,} ({row['I']/N*100:.1f}%)  "
        f"R={row['R']:,} ({row['R']/N*100:.0f}%)"
    )

# ── 4. Model parameters ───────────────────────────────────────────────────────
# beta = R0 / infectious_period = 5 / 10 = 0.5 per day
nticks = 365

params = PropertySet({
    "prng_seed": 42,
    "nticks":    nticks,
    "beta":      0.5,
})

# ── 4b. Seasonal forcing profile ──────────────────────────────────────────────
# Piecewise-linear annual cycle (all-patches identical):
#   Days   0– 90: winter peak  → 1.3× beta
#   Days  90–150: linear ramp  → 1.3 → 0.7 (60-day transition)
#   Days 150–240: summer trough → 0.7× beta
#   Days 240–365: linear ramp  → 0.7 → 1.3 (125-day transition back to winter)
#
# Mean is exactly 1.0 by construction:
#   both ramps average (1.3+0.7)/2 = 1.0, so the area-weighted mean = 1.0
# This preserves beta as the true annual-average transmission rate.

WINTER_MULT = 1.3
SUMMER_MULT = 0.7

days = np.arange(nticks, dtype=np.float64)
season = np.where(
    days < 90,
    WINTER_MULT,
    np.where(
        days < 150,
        WINTER_MULT + (SUMMER_MULT - WINTER_MULT) * (days - 90.0) / 60.0,
        np.where(
            days < 240,
            SUMMER_MULT,
            SUMMER_MULT + (WINTER_MULT - SUMMER_MULT) * (days - 240.0) / 125.0,
        ),
    ),
)

# Normalise to mean = 1.0  (should already be 1.0; guard against float drift)
season /= season.mean()

assert abs(season.mean() - 1.0) < 1e-10,    f"Seasonal mean = {season.mean():.6f}, expected 1.0"
assert np.all(season >= 0.0),                 "Seasonal multiplier contains negative values"
assert abs(season[:90].mean() - WINTER_MULT) < 0.01, "Winter window not at expected level"
assert abs(season[150:240].mean() - SUMMER_MULT) < 0.01, "Summer window not at expected level"

seasonality = ValuesMap.from_timeseries(season, n_patches)

print(f"\nSeasonal profile  (normalised mean = {season.mean():.6f})")
print(f"  Winter (days 0-90):    {season[:90].mean():.3f}× beta  "
      f"(effective beta = {params.beta * season[:90].mean():.4f})")
print(f"  Summer (days 150-240): {season[150:240].mean():.3f}× beta  "
      f"(effective beta = {params.beta * season[150:240].mean():.4f})")

# ── 5. Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.constant_int(4)    # latent period:     exactly 4 days
infdurdist = dists.constant_int(10)   # infectious period: exactly 10 days

# ── 6. Construct model and assign gravity network ─────────────────────────────
model = Model(scenario, params)
model.network = network   # replace default zero matrix with gravity coupling

# ── 7. Component list ─────────────────────────────────────────────────────────
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),  # ← seasonal
]

# ── 8. Run ────────────────────────────────────────────────────────────────────
print("\nRunning 4-patch SEIR with gravity coupling + seasonal forcing (365 days)...")
model.run("GravitySEIR_4patch_Seasonal")
print("Done.\n")

# ── 9. Epidemic summary ───────────────────────────────────────────────────────
print("=== Epidemic Summary ===")
for idx in range(n_patches):
    N         = int(populations[idx])
    I_trace   = model.nodes.I[:nticks, idx]
    peak_I    = int(I_trace.max())
    peak_day  = int(I_trace.argmax())
    S_final   = int(model.nodes.S[nticks - 1, idx])
    total_inf = int(model.nodes.newly_infected[:nticks, idx].sum())
    attack_rt = (N - S_final) / N * 100
    print(
        f"  Patch {idx} (N={N:,}): "
        f"peak I = {peak_I:,} on day {peak_day} | "
        f"total infected = {total_inf:,} | "
        f"attack rate = {attack_rt:.1f}%"
    )

for name, arr in [("S", model.nodes.S), ("E", model.nodes.E),
                  ("I", model.nodes.I), ("R", model.nodes.R)]:
    assert np.all(arr[:nticks, :] >= 0), f"Negative values in {name} compartment"
print("\nSanity check passed: all compartments non-negative")

# ── 10. Figures ───────────────────────────────────────────────────────────────
# Layout: 3 rows × 2 cols
#   Row 1: migration heatmap  |  seasonal forcing profile
#   Row 2: Patch 0 SEIR curves  |  Patch 1 SEIR curves
#   Row 3: Patch 2 SEIR curves  |  Patch 3 SEIR curves

PATCH_COLORS = ["steelblue", "darkorange", "forestgreen", "crimson"]

fig = plt.figure(figsize=(14, 11))
fig.suptitle(
    "4-Patch SEIR — Gravity Coupling + Seasonal Forcing\n"
    r"R$_0$=5, Latent=4d, Infectious=10d  |  k=0.01, a=1, b=1, c=1.5  |  "
    "Winter 1.3×, Summer 0.7×",
    fontsize=13,
)

# Panel 1 — migration network heatmap (origin = row, destination = column)
ax_net = fig.add_subplot(3, 2, 1)
im = ax_net.imshow(network * 100, cmap="YlOrRd", aspect="auto", vmin=0)
plt.colorbar(im, ax=ax_net, label="% of pop / day")
tick_labels = [f"P{j}\n({populations[j]//1000}k)" for j in range(n_patches)]
ax_net.set_xticks(range(n_patches))
ax_net.set_yticks(range(n_patches))
ax_net.set_xticklabels(tick_labels, fontsize=8)
ax_net.set_yticklabels([f"P{i}" for i in range(n_patches)], fontsize=8)
ax_net.set_xlabel("Destination", fontsize=8)
ax_net.set_ylabel("Origin", fontsize=8)
ax_net.set_title("Daily migration rates\n(row = origin patch)", fontsize=9)
for i in range(n_patches):
    for j in range(n_patches):
        txt_color = "white" if network[i, j] > 0.09 else "black"
        ax_net.text(j, i, f"{network[i, j]*100:.2f}%",
                    ha="center", va="center", fontsize=7, color=txt_color)

# Panel 2 — seasonal forcing profile
ax_sea = fig.add_subplot(3, 2, 2)
t = np.arange(nticks)
ax_sea.plot(t, season, color="indigo", linewidth=2.0, label="β multiplier")
ax_sea.axhline(1.0, color="grey", linestyle=":", linewidth=1.0, label="Baseline (1.0×)")
ax_sea.axhline(WINTER_MULT, color="cornflowerblue", linestyle="--",
               linewidth=1.0, label=f"Winter peak ({WINTER_MULT}×)")
ax_sea.axhline(SUMMER_MULT, color="tomato", linestyle="--",
               linewidth=1.0, label=f"Summer trough ({SUMMER_MULT}×)")
# Shade the fixed windows
ax_sea.axvspan(0,   90,  alpha=0.12, color="cornflowerblue", label="Winter window (days 0–90)")
ax_sea.axvspan(150, 240, alpha=0.12, color="tomato",         label="Summer window (days 150–240)")
ax_sea.set_xlim(0, nticks - 1)
ax_sea.set_ylim(0.55, 1.45)
ax_sea.set_xlabel("Day")
ax_sea.set_ylabel("β multiplier")
ax_sea.set_title(
    f"Seasonal forcing profile\n"
    f"mean={season.mean():.4f}, "
    f"eff. β: winter={params.beta*WINTER_MULT:.3f}, summer={params.beta*SUMMER_MULT:.3f}",
    fontsize=9,
)
ax_sea.legend(fontsize=7, loc="lower right", ncol=2)
ax_sea.grid(alpha=0.3)

# Panels 3–6 — SEIR curves per patch  (with seasonal window shading)
for idx in range(n_patches):
    ax = fig.add_subplot(3, 2, idx + 3)
    N   = populations[idx]
    pct = 100.0 / N

    # Background shading to indicate seasonal windows
    ax.axvspan(0,   90,  alpha=0.07, color="cornflowerblue")   # winter
    ax.axvspan(150, 240, alpha=0.07, color="tomato")            # summer

    ax.plot(t, model.nodes.S[:nticks, idx] * pct,
            label="S", color="steelblue",   linewidth=1.8)
    ax.plot(t, model.nodes.E[:nticks, idx] * pct,
            label="E", color="darkorange",  linewidth=1.8)
    ax.plot(t, model.nodes.I[:nticks, idx] * pct,
            label="I", color="crimson",     linewidth=1.8)
    ax.plot(t, model.nodes.R[:nticks, idx] * pct,
            label="R", color="forestgreen", linewidth=1.8)

    peak_day = int(model.nodes.I[:nticks, idx].argmax())
    peak_I   = int(model.nodes.I[:nticks, idx].max())
    ax.axvline(peak_day, color="crimson", linestyle="--", linewidth=1, alpha=0.45)

    ax.set_title(
        f"Patch {idx}  (N={N:,})\nPeak I = {peak_I:,} on day {peak_day}",
        fontsize=9,
    )
    ax.set_xlabel("Day")
    ax.set_ylabel("% of population")
    ax.set_xlim(0, nticks - 1)
    ax.set_ylim(0, 100)
    ax.legend(loc="center right", fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("seir_4patch_gravity_seasonal.png", dpi=150, bbox_inches="tight")
plt.show()
print("Plot saved → seir_4patch_gravity_seasonal.png")
```

**What changed from the previous version:**

1. **Import** — `from laser.generic.utils import ValuesMap` added.

2. **Section 4b (new)** — Builds a piecewise-linear 365-day seasonal multiplier array: constant at 1.3× for days 0–90, linear ramp down to 0.7× by day 150, constant at 0.7× through day 240, then linear ramp back to 1.3× by day 365. The mean is exactly 1.0 by construction (both 60-day and 125-day ramps average (1.3+0.7)/2 = 1.0), so `beta=0.5` still represents the annual-average transmission rate. The array is wrapped into a `ValuesMap` via `from_timeseries`.

3. **Section 7** — `SEIR.Transmission(..., seasonality=None)` → `SEIR.Transmission(..., seasonality=seasonality)`.

4. **Figure panel 2** — Replaced the outflow bar chart with a seasonal profile plot showing the multiplier curve, the fixed-window shading, and the effective `beta` at winter/summer extremes. The SEIR panels (3–6) also get faint blue/red background bands marking the winter and summer windows so the epidemic timing can be read against the forcing.
