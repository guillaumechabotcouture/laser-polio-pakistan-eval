Now I have the exact details I need. Key finding: `model.network[i,j]` is the **fraction of patch i's FOI exported to patch j**, so row sums should equal the export fraction (≤ 0.15) — **not** 1.0 like a row-stochastic matrix. The FOI formula is:

```
ft[j] = local_FOI[j] + Σ_i(ft[i]·network[i,j]) − Σ_j(ft[i]·network[i,j])
```

The gravity flow fraction is `F_ij / N_i` (raw flow divided by source pop), then rows are capped at 15%.

```python
"""
4-patch SEIR model with gravity-model spatial coupling (laser-generic v1.0.0).

Disease parameters:
  R0 = 5, latent period = 4 days, infectious period = 10 days

Initial conditions (per patch):
  90% Susceptible | 1% Infectious | 9% Recovered | 0% Exposed

Spatial coupling: gravity model  F_ij = k · N_i^a · N_j^b / d_ij^c
  k=0.01, a=1, b=1, c=1.5
  Patches in a line, 75 km apart.
  Row-normalized to cap each patch's total FOI export at 15%.

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
#
# With a=1 the N_i in the numerator and denominator cancel:
#   f_ij = k · N_j^b / d_ij^c
# (kept in the general form so changing a later just works).
export_frac_matrix = F_raw / populations[:, None]

# Row-cap: if a patch's total export exceeds 15%, scale that row down uniformly.
max_export_frac = 0.15
export_total    = export_frac_matrix.sum(axis=1)          # total export per patch
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

# ── 4. Scenario GeoDataFrame ───────────────────────────────────────────────────
# Geometry reflects the 1-D spatial arrangement (lon fixed, lat offset by km).
# 1 degree latitude ≈ 111 km.
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

# ── 5. Model parameters ────────────────────────────────────────────────────────
params = PropertySet(
    {
        "prng_seed": 42,
        "nticks":    nticks,
        "beta":      beta,
    }
)

# ── 6. Duration distributions ──────────────────────────────────────────────────
# gamma(shape, scale) → mean = shape * scale, variance = shape * scale²
#
# Exposed:    shape=16, scale=0.25 → mean = 4.0 days, std ≈ 1.0 day
# Infectious: shape=40, scale=0.25 → mean = 10.0 days, std ≈ 1.6 days
expdurdist = dists.gamma(shape=16, scale=0.25)   # latent period
infdurdist = dists.gamma(shape=40, scale=0.25)   # infectious period

# ── 7. Build model and install gravity network ─────────────────────────────────
# birthrates=None → static population (no births/deaths).
model = Model(scenario, params, birthrates=None)

# Install gravity-model network.
# network[i,j] = fraction of patch i's FOI exported to patch j each tick.
# Row sums are the export fractions (≤ 0.15 by construction above).
# The TransmissionSE FOI update is:
#   ft[j] += Σ_i ft[i]·network[i,j]   (incoming)
#   ft[i] -= Σ_j ft[i]·network[i,j]   (outgoing)
model.network = network
assert model.network.sum() > 0, "Network must be non-zero for spatial coupling"

# ── 8. Assemble components ─────────────────────────────────────────────────────
# Ordering: Susceptible first and Recovered last ensure S+E+I+R = N invariant.
# Transmission runs after Recovered so it reads the fully-updated compartments.
model.components = [
    SEIR.Susceptible(model),                          # propagates S[t] → S[t+1]
    SEIR.Exposed(model, expdurdist, infdurdist),      # E→I transitions
    SEIR.Infectious(model, infdurdist),               # I→R transitions
    SEIR.Recovered(model),                            # propagates R[t] → R[t+1]
    SEIR.Transmission(model, expdurdist, seasonality=None),  # S→E infections
]

# ── 9. Run ─────────────────────────────────────────────────────────────────────
print("\nRunning simulation ...")
model.run("4-Patch SEIR with Gravity Coupling")
print("Done.\n")

# ── 10. Post-run verification ──────────────────────────────────────────────────
S = model.nodes.S[:nticks, :]    # shape (nticks, nnodes)
E = model.nodes.E[:nticks, :]
I = model.nodes.I[:nticks, :]
R = model.nodes.R[:nticks, :]
N = S + E + I + R                # should equal populations every tick

# 10a. Population conservation (no demographics → N constant per patch)
for i in range(nnodes):
    N_i = N[:, i]
    assert N_i.min() == N_i.max(), (
        f"Population not conserved in Patch {i}: "
        f"min={N_i.min():,}, max={N_i.max():,}"
    )

# 10b. Non-negativity
assert np.all(S >= 0) and np.all(E >= 0) and np.all(I >= 0) and np.all(R >= 0), \
    "Negative compartment count detected"

# 10c. Epidemic occurred in every patch
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

# ── 11. Diagnostic plots ───────────────────────────────────────────────────────
days         = np.arange(nticks)
patch_labels = [f"Patch {i}  (N={populations[i]//1000}k)" for i in range(nnodes)]
colors       = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle(
    f"4-Patch SEIR with Gravity Coupling  |  "
    f"R₀={R0}, β={beta}  |  "
    f"k={gravity_k}, a={gravity_a}, b={gravity_b}, c={gravity_c}  |  "
    f"max export={max_export_frac:.0%}",
    fontsize=11,
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
plt.savefig("seir_4patch_gravity.png", dpi=100, bbox_inches="tight")
plt.show()
print("\nPlot saved: seir_4patch_gravity.png")
```

**Key changes from the baseline:**

**New section 3 — gravity network:**
- `positions_km = [0, 75, 150, 225]` — patches in a line, 75 km apart
- Raw flows: `F_ij = k · N_i^a · N_j^b / d_ij^c` with your parameters
- Divided by `N_i` to get dimensionless export fractions (`f_ij = F_ij / N_i`); with `a=1` this simplifies to `k · N_j / d_ij^1.5`
- `cap_scale` clips any row exceeding 15% down to exactly 15% — with these parameters all 4 patches hit the cap (raw exports would be ~370–430% of population)

**Section 7 — network assignment:**
- `model.network = network` replaces `model.network = np.zeros(...)` and the zero-assertion
- Added a comment explaining the FOI update formula that LASER's `TransmissionSE` applies

**Section 4 — geometry:** lats now track `positions_km / 111.0` so the GeoDataFrame geometry is consistent with the stated 75 km spacing.
