"""
Basic 4-patch SEIR model using the LASER framework (laser-generic v1.0.0).

Disease parameters:
  R0 = 5, latent period = 4 days, infectious period = 10 days

Initial conditions (per patch):
  90% Susceptible | 1% Infectious | 9% Recovered | 0% Exposed

No spatial coupling, no seasonal forcing, no demographics.
Simulation duration: 1 year (365 days).
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model
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

S_init = np.round(0.90 * populations).astype(np.int32)
I_init = np.round(0.01 * populations).astype(np.int32)
R_init = np.round(0.09 * populations).astype(np.int32)
E_init = populations - S_init - I_init - R_init   # = 0 (90+1+9 = 100%)

# Verify
assert np.all(S_init + E_init + I_init + R_init == populations), \
    "S + E + I + R must equal N in every patch"
assert np.all(I_init > 0), \
    "Every patch must have at least one initial infectious agent"

print(
    f"Initial conditions verified:\n"
    f"  Total N={populations.sum():,}, "
    f"S={S_init.sum():,} (90%), "
    f"I={I_init.sum():,} (1%), "
    f"R={R_init.sum():,} (9%), "
    f"E={E_init.sum():,} (0%)"
)

# ── 3. Scenario GeoDataFrame ───────────────────────────────────────────────────
# Geometry is required by Model, but arbitrary here — spatial coupling is off.
lats = np.array([30.0, 31.0, 32.0, 33.0])
lons = np.array([70.0, 71.0, 72.0, 73.0])

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

# ── 4. Model parameters ────────────────────────────────────────────────────────
# No gravity_k/b/c → model.network will not be auto-populated; we zero it below.
params = PropertySet(
    {
        "prng_seed": 42,
        "nticks":    nticks,
        "beta":      beta,
    }
)

# ── 5. Duration distributions ──────────────────────────────────────────────────
# gamma(shape, scale) → mean = shape * scale, variance = shape * scale²
#
# Exposed:    shape=16, scale=0.25 → mean = 4.0 days, std ≈ 1.0 day
# Infectious: shape=40, scale=0.25 → mean = 10.0 days, std ≈ 1.6 days
expdurdist = dists.gamma(shape=16, scale=0.25)   # latent period
infdurdist = dists.gamma(shape=40, scale=0.25)   # infectious period

# ── 6. Build model ─────────────────────────────────────────────────────────────
# birthrates=None → static population (no births/deaths).
model = Model(scenario, params, birthrates=None)

# Disable spatial coupling: set migration matrix to all zeros.
# TransmissionSE uses model.network in its FOI calculation; zeroing it makes
# each patch evolve independently.
model.network = np.zeros((nnodes, nnodes), dtype=np.float64)
assert model.network.sum() == 0.0, "Network must be zero for no spatial coupling"

# ── 7. Assemble components ─────────────────────────────────────────────────────
# Ordering: Susceptible first and Recovered last ensure S+E+I+R = N invariant.
# Transmission runs after Recovered so it reads the fully-updated compartments.
model.components = [
    SEIR.Susceptible(model),                          # propagates S[t] → S[t+1]
    SEIR.Exposed(model, expdurdist, infdurdist),      # E→I transitions
    SEIR.Infectious(model, infdurdist),               # I→R transitions
    SEIR.Recovered(model),                            # propagates R[t] → R[t+1]
    SEIR.Transmission(model, expdurdist, seasonality=None),  # S→E infections
]

# ── 8. Run ─────────────────────────────────────────────────────────────────────
print("\nRunning simulation ...")
model.run("Basic 4-Patch SEIR")
print("Done.\n")

# ── 9. Post-run verification ───────────────────────────────────────────────────
S = model.nodes.S[:nticks, :]    # shape (nticks, nnodes)
E = model.nodes.E[:nticks, :]
I = model.nodes.I[:nticks, :]
R = model.nodes.R[:nticks, :]
N = S + E + I + R                # should equal populations every tick

# 9a. Population conservation (no demographics → N constant per patch)
for i in range(nnodes):
    N_i = N[:, i]
    assert N_i.min() == N_i.max(), (
        f"Population not conserved in Patch {i}: "
        f"min={N_i.min():,}, max={N_i.max():,}"
    )

# 9b. Non-negativity
assert np.all(S >= 0) and np.all(E >= 0) and np.all(I >= 0) and np.all(R >= 0), \
    "Negative compartment count detected"

# 9c. Epidemic occurred in every patch
cumulative_infections = model.nodes.newly_infected[:nticks, :].sum(axis=0)
assert np.all(cumulative_infections > 0), (
    f"No infections in patch(es): {np.where(cumulative_infections == 0)[0]}"
)

print("Verification passed.")
print(f"{'Patch':<10} {'Population':>12} {'Total Infected':>16} {'Attack Rate':>12}")
print("-" * 52)
for i in range(nnodes):
    ar = cumulative_infections[i] / populations[i]
    print(
        f"Patch {i:<4}  {populations[i]:>12,}  "
        f"{cumulative_infections[i]:>14,.0f}  {ar:>11.1%}"
    )
total_pop = populations.sum()
total_inf = cumulative_infections.sum()
print("-" * 52)
print(f"{'Total':<10} {total_pop:>12,}  {total_inf:>14,.0f}  {total_inf/total_pop:>11.1%}")

# ── 10. Diagnostic plots ───────────────────────────────────────────────────────
days         = np.arange(nticks)
patch_labels = [f"Patch {i}  (N={populations[i]//1000}k)" for i in range(nnodes)]
colors       = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle(
    f"Basic 4-Patch SEIR  |  R₀={R0}, latent={latent_days}d, "
    f"infectious={infectious_days}d, β={beta}",
    fontsize=13,
)

# Panel A: daily new infections by patch
ax = axes[0, 0]
for i in range(nnodes):
    ax.plot(days, model.nodes.newly_infected[:nticks, i],
            label=patch_labels[i], color=colors[i])
ax.set_title("Daily New Infections by Patch")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel B: susceptible fraction S/N by patch
ax = axes[0, 1]
for i in range(nnodes):
    ax.plot(days, S[:, i] / populations[i],
            label=patch_labels[i], color=colors[i])
ax.set_title("Susceptible Fraction  S/N  by Patch")
ax.set_xlabel("Day")
ax.set_ylabel("S / N")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel C: aggregate SEIR compartments (all patches combined)
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

# Panel D: total daily incidence across all patches
ax = axes[1, 1]
ax.plot(days, model.nodes.newly_infected[:nticks, :].sum(axis=1), color="black")
ax.set_title("Total Daily New Infections  (All Patches)")
ax.set_xlabel("Day")
ax.set_ylabel("New Infections")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("seir_4patch_basic.png", dpi=100, bbox_inches="tight")
plt.show()
print("\nPlot saved: seir_4patch_basic.png")