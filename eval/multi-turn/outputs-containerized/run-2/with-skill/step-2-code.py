"""
4-Patch Gravity SEIR Model
==========================
Disease: respiratory (R0=5, latent=4d, infectious=10d)
Patches: pop 100k / 200k / 150k / 80k, arranged in a line 75 km apart
Gravity: k=0.01, a=1, b=1, c=1.5; row-normalised, max 15% export per patch
Run:     365 days — gravity spatial coupling, no vital dynamics, no seasonal forcing
Package: laser-generic v1.0.0
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
import matplotlib.pyplot as plt

# ── Simulation settings ────────────────────────────────────────────────────
NTICKS = 365
POPS   = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)
NNODES = len(POPS)
SEED   = 42

# ── Disease parameters ─────────────────────────────────────────────────────
R0       = 5.0
LATENT_D = 4
INFECT_D = 10
BETA     = R0 / INFECT_D  # = 0.5 per day

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
# Approximate lon offsets: ~0.675° per 75 km at equator (dummy positions)
PATCH_LON = np.array([0.000, 0.675, 1.351, 2.027])
scenario = gpd.GeoDataFrame({
    "nodeid":     np.arange(NNODES),
    "name":       [f"patch_{i}" for i in range(NNODES)],
    "population": POPS,
    "S": S_init, "E": E_init, "I": I_init, "R": R_init,
    "geometry":   [Point(lon, 0.0) for lon in PATCH_LON],
}, crs="EPSG:4326")

# ── Model parameters ───────────────────────────────────────────────────────
params = PropertySet({
    "prng_seed": SEED,
    "nticks":    NTICKS,
    "beta":      BETA,
})

# ── Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.gamma(shape=16, scale=0.25)   # E→I: mean = 4 d
infdurdist = dists.gamma(shape=40, scale=0.25)   # I→R: mean = 10 d

# ── Gravity coupling network ───────────────────────────────────────────────
# Patches are equally spaced along a line, 75 km apart.
# Gravity formula: W_ij = k * N_i^a * N_j^b / d_ij^c
# Dividing by N_i converts W to outflow fractions per capita.
# Each row is then capped so total export ≤ MAX_EXPORT.

PATCH_X_KM = np.array([0.0, 75.0, 150.0, 225.0])
GRAV_K, GRAV_A, GRAV_B, GRAV_C = 0.01, 1.0, 1.0, 1.5
MAX_EXPORT = 0.15

# Pairwise distances (km); inf on diagonal prevents self-coupling
DIST_KM = np.abs(PATCH_X_KM[:, None] - PATCH_X_KM[None, :])
np.fill_diagonal(DIST_KM, np.inf)

# Raw gravity weights: W_ij = k * N_i^a * N_j^b / d_ij^c
W = (GRAV_K
     * np.power(POPS[:, None].astype(np.float64), GRAV_A)
     * np.power(POPS[None, :].astype(np.float64), GRAV_B)
     / np.power(DIST_KM, GRAV_C))
np.fill_diagonal(W, 0.0)

# Convert to per-capita outflow fractions: frac_ij = W_ij / N_i
network = W / POPS[:, None].astype(np.float64)

# Row-normalise: scale each row so total export fraction ≤ MAX_EXPORT
row_sums = network.sum(axis=1)
scale = np.minimum(1.0, MAX_EXPORT / row_sums)
network = network * scale[:, None]

# ── Validate network ───────────────────────────────────────────────────────
assert np.all(network >= 0),                       "Negative network weight"
assert np.all(network.sum(axis=1) <= MAX_EXPORT + 1e-9), \
    "Export fraction exceeds 15% cap"
assert network.sum() > 0,                          "Network is all zeros"

print("Gravity network row sums (export fractions):")
for i in range(NNODES):
    row = network[i]
    print(f"  Patch {i} → {[f'{row[j]:.4f}' for j in range(NNODES)]}  "
          f"(total export = {row.sum():.4f})")
print()

# ── Construct model ────────────────────────────────────────────────────────
model = Model(scenario, params, birthrates=None)
model.network = network

# ── Assemble components ────────────────────────────────────────────────────
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=None),
]

# ── Run ────────────────────────────────────────────────────────────────────
print(f"Running gravity 4-patch SEIR model")
print(f"  β={BETA:.3f} d⁻¹,  R₀={R0},  latent={LATENT_D}d,  infectious={INFECT_D}d")
print(f"  Gravity: k={GRAV_K}, a={GRAV_A}, b={GRAV_B}, c={GRAV_C}")
print(f"  Spacing: 75 km;  max export fraction: {MAX_EXPORT:.0%}")
print(f"  Patches (pop): {POPS.tolist()}")
print(f"  Init: {F_S:.0%} S / {F_E:.0%} E / {F_I:.0%} I / {F_R:.0%} R")
model.run("Gravity 4-Patch SEIR")
print("Simulation complete.\n")

# ── Extract compartment time series ───────────────────────────────────────
days = np.arange(NTICKS)
S = model.nodes.S[:NTICKS, :]
E = model.nodes.E[:NTICKS, :]
I = model.nodes.I[:NTICKS, :]
R = model.nodes.R[:NTICKS, :]
N = S + E + I + R

newly_infected = model.nodes.newly_infected[:NTICKS, :].sum(axis=1)

# ── Validation ─────────────────────────────────────────────────────────────
assert np.all(S >= 0), "Negative S"
assert np.all(E >= 0), "Negative E"
assert np.all(I >= 0), "Negative I"
assert np.all(R >= 0), "Negative R"

for p in range(NNODES):
    assert I[:, p].max() > I_init[p], \
        f"No epidemic growth in patch {p} (pop={POPS[p]:,})"

# No vital dynamics → per-patch population must be constant
for p in range(NNODES):
    assert np.all(N[:, p] == POPS[p]), \
        f"Population changed in patch {p}"

print("All validation checks passed.\n")

# ── Summary statistics ─────────────────────────────────────────────────────
print("Summary:")
print(f"  Peak total infectious : {I.sum(axis=1).max():,}")
overall_ar = (S_init.sum() - S[-1, :].sum()) / POPS.sum()
print(f"  Overall attack rate   : {overall_ar:.1%}")
print()
for p in range(NNODES):
    peak_ip  = I[:, p].max()
    day_peak = int(I[:, p].argmax())
    ar_p     = (S_init[p] - S[-1, p]) / POPS[p]
    print(f"  Patch {p} (N={POPS[p]//1000}k): "
          f"peak I = {peak_ip:>6,}  day {day_peak:3d},  "
          f"attack rate = {ar_p:.1%}")

# ── Figure 1: SEIR compartments per patch ─────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
fig.suptitle(
    f"4-Patch SEIR — Gravity Spatial Coupling\n"
    f"R₀ = {R0},  β = {BETA:.2f} d⁻¹,  "
    f"Latent = {LATENT_D} d,  Infectious = {INFECT_D} d\n"
    f"Gravity k={GRAV_K}, a={GRAV_A}, b={GRAV_B}, c={GRAV_C};  "
    f"spacing = 75 km;  max export = {MAX_EXPORT:.0%}",
    fontsize=12,
)

STYLES = {
    "S": dict(color="steelblue",  lw=2.0, ls="-",  label="S"),
    "E": dict(color="darkorange", lw=1.5, ls="--", label="E"),
    "I": dict(color="firebrick",  lw=2.0, ls="-",  label="I"),
    "R": dict(color="seagreen",   lw=1.5, ls=":",  label="R"),
}

for p, ax in enumerate(axes.flat):
    for key, arr in (("S", S), ("E", E), ("I", I), ("R", R)):
        ax.plot(days, arr[:, p], **STYLES[key])
    ax.set_title(f"Patch {p}  (N = {POPS[p]:,})", fontsize=11)
    ax.set_ylabel("Count")
    ax.legend(fontsize=9, loc="center right")
    ax.grid(alpha=0.3)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):,}")
    )

for ax in axes[1]:
    ax.set_xlabel("Day")

plt.tight_layout()
plt.savefig("seir_4patch_gravity.png", dpi=120)
print("\nPlot saved → seir_4patch_gravity.png")
plt.show()

# ── Figure 2: Aggregate daily incidence ───────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.fill_between(days, newly_infected, alpha=0.25, color="firebrick")
ax2.plot(days, newly_infected, color="firebrick", lw=1.8)
ax2.set_title(
    f"Daily New Infections — All 4 Patches Combined  "
    f"(R₀ = {R0}, Gravity Coupling)",
    fontsize=12,
)
ax2.set_xlabel("Day")
ax2.set_ylabel("New infections per day")
ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("seir_4patch_gravity_incidence.png", dpi=120)
print("Plot saved → seir_4patch_gravity_incidence.png")
plt.show()

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
                 fontsize=9, color="black" if network[i, j] < MAX_EXPORT * 0.7 else "white")
plt.tight_layout()
plt.savefig("seir_4patch_gravity_network.png", dpi=120)
print("Plot saved → seir_4patch_gravity_network.png")
plt.show()