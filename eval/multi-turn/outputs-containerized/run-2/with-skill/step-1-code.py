"""
4-Patch Basic SEIR Model
========================
Disease: respiratory (R0=5, latent=4d, infectious=10d)
Patches: pop 100k / 200k / 150k / 80k
Init:    90% S, 0% E, 1% I, 9% R
Run:     365 days — no spatial coupling, no vital dynamics, no seasonal forcing
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
R0       = 5.0          # basic reproduction number
LATENT_D = 4            # mean exposed (E) period, days
INFECT_D = 10           # mean infectious (I) period, days
BETA     = R0 / INFECT_D  # = 0.5 per day

# ── Initial compartment fractions ──────────────────────────────────────────
F_S, F_E, F_I, F_R = 0.90, 0.00, 0.01, 0.09
assert abs(F_S + F_E + F_I + F_R - 1.0) < 1e-9, "Fractions must sum to 1"

# ── Compute initial counts (integer; rounding absorbed into S) ─────────────
I_init = np.maximum(1, np.round(F_I * POPS).astype(np.int32))
R_init = np.round(F_R * POPS).astype(np.int32)
E_init = np.zeros(NNODES, dtype=np.int32)
S_init = POPS - E_init - I_init - R_init

assert np.all(S_init + E_init + I_init + R_init == POPS), \
    "S+E+I+R must equal population in every patch"
assert np.all(I_init > 0),   "At least one infectious agent required per patch"
assert np.all(S_init >= 0),  "Negative S — check initial fractions"

# ── Scenario GeoDataFrame ──────────────────────────────────────────────────
# Dummy coordinates — spatial coupling is disabled, positions are irrelevant
scenario = gpd.GeoDataFrame({
    "nodeid":     np.arange(NNODES),
    "name":       [f"patch_{i}" for i in range(NNODES)],
    "population": POPS,
    "S": S_init, "E": E_init, "I": I_init, "R": R_init,
    "geometry":   [Point(float(i), 0.0) for i in range(NNODES)],
}, crs="EPSG:4326")

# ── Model parameters ───────────────────────────────────────────────────────
params = PropertySet({
    "prng_seed": SEED,
    "nticks":    NTICKS,
    "beta":      BETA,
})

# ── Duration distributions ─────────────────────────────────────────────────
# gamma(shape, scale): mean = shape × scale
#   Latent period  4 d  →  shape=16, scale=0.25  (mean = 16 × 0.25 = 4 d)
#   Infect period 10 d  →  shape=40, scale=0.25  (mean = 40 × 0.25 = 10 d)
# High shape → low CV → near-deterministic durations
expdurdist = dists.gamma(shape=16, scale=0.25)   # E→I transition
infdurdist = dists.gamma(shape=40, scale=0.25)   # I→R transition

# ── Construct model ────────────────────────────────────────────────────────
model = Model(scenario, params, birthrates=None)

# Explicitly zero the coupling network — no spatial coupling in this step
model.network = np.zeros((NNODES, NNODES), dtype=np.float64)

# ── Assemble components (order matters for S+E+I+R = N invariant) ──────────
# Susceptible and Recovered wrap the transition steps each tick
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist, seasonality=None),  # no seasonal forcing
]

# ── Run ────────────────────────────────────────────────────────────────────
print(f"Running basic 4-patch SEIR model")
print(f"  β={BETA:.3f} d⁻¹,  R₀={R0},  latent={LATENT_D}d,  infectious={INFECT_D}d")
print(f"  Patches (pop): {POPS.tolist()}")
print(f"  Init: {F_S:.0%} S / {F_E:.0%} E / {F_I:.0%} I / {F_R:.0%} R")
model.run("Basic 4-Patch SEIR")
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
assert np.all(S >= 0), "Negative S detected — compartment depletion bug"
assert np.all(E >= 0), "Negative E detected"
assert np.all(I >= 0), "Negative I detected"
assert np.all(R >= 0), "Negative R detected"

# Every patch should show epidemic growth
for p in range(NNODES):
    assert I[:, p].max() > I_init[p], \
        f"No epidemic growth in patch {p} (pop={POPS[p]:,})"

# With no births/deaths the per-patch population must be constant
for p in range(NNODES):
    assert np.all(N[:, p] == POPS[p]), \
        f"Population changed in patch {p} — unexpected coupling or births"

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
    f"4-Patch SEIR — Basic (No Spatial Coupling)\n"
    f"R₀ = {R0},  β = {BETA:.2f} d⁻¹,  "
    f"Latent = {LATENT_D} d,  Infectious = {INFECT_D} d",
    fontsize=13,
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
plt.savefig("seir_4patch_basic.png", dpi=120)
print("\nPlot saved → seir_4patch_basic.png")
plt.show()

# ── Figure 2: Aggregate daily incidence ───────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.fill_between(days, newly_infected, alpha=0.25, color="firebrick")
ax2.plot(days, newly_infected, color="firebrick", lw=1.8)
ax2.set_title(
    f"Daily New Infections — All 4 Patches Combined  (R₀ = {R0})",
    fontsize=12,
)
ax2.set_xlabel("Day")
ax2.set_ylabel("New infections per day")
ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("seir_4patch_incidence.png", dpi=120)
print("Plot saved → seir_4patch_incidence.png")
plt.show()