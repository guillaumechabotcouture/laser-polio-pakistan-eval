The workspace file system is read-only, so I cannot write to disk. Here is the complete, self-contained model code matching your specification exactly:

```python
#!/usr/bin/env python3
"""
Republic of Eastland — Spatial SEIR Model (Respiratory Illness)
4-patch spatial simulation using the LASER framework (laser-generic).

Specification
─────────────
  Patches:       4 (100k, 200k, 150k, 80k population)
  Disease:       R0 ~ 5
                 Latent     gamma(shape=4, scale=1)  → mean 4 days
                 Infectious gamma(shape=5, scale=2)  → mean 10 days
  Init:          90% S, 0% E, 1% I, 9% R per patch
  Demographics:  BirthsByCBR  CBR = 30 per 1000/year
                 MortalityByCDR CDR = 10 per 1000/year
  Capacity:      calc_capacity pre-allocates for 10 years of births
  Duration:      10 years (nticks = 3650)
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.core.utils import calc_capacity
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Constants ─────────────────────────────────────────────────────────────────

NTICKS = 3650           # 10 years (days)
NNODES = 4

# Disease: R0 = beta * mean_infectious_period = beta * (5 * 2) = beta * 10
# → beta = R0 / 10 = 5 / 10 = 0.5
BETA      = 0.5
EXP_SHAPE = 4           # latent gamma shape  (mean = shape * scale = 4 days)
EXP_SCALE = 1.0         # latent gamma scale
INF_SHAPE = 5           # infectious gamma shape (mean = 5 * 2 = 10 days)
INF_SCALE = 2.0         # infectious gamma scale

CBR = 30.0              # crude birth rate  (per 1000/year)
CDR = 10.0              # crude death rate  (per 1000/year)

PATCH_NAMES = ["Eastport", "Northton", "Centerburg", "Westham"]


# ── Patch configuration ───────────────────────────────────────────────────────

populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int32)

# Synthetic coordinates inside the Republic of Eastland (decimal degrees)
lats = np.array([30.0, 30.5, 31.0, 31.5])
lons = np.array([50.0, 50.5, 51.0, 51.5])


# ── Initial compartments: S=90%, E=0%, I=1%, R=9% ───────────────────────────

I_init = np.round(0.01 * populations).astype(np.int32)
R_init = np.round(0.09 * populations).astype(np.int32)
E_init = np.zeros(NNODES, dtype=np.int32)
S_init = populations - E_init - I_init - R_init

assert (S_init + E_init + I_init + R_init == populations).all(), \
    "S+E+I+R must equal population in every patch"
assert (I_init > 0).any(), "At least one patch needs initial infections"

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(NNODES, dtype=np.int32),
        "name":       PATCH_NAMES,
        "population": populations,
        "S":          S_init,
        "E":          E_init,
        "I":          I_init,
        "R":          R_init,
        "geometry":   [Point(lon, lat) for lat, lon in zip(lats, lons)],
    },
    crs="EPSG:4326",
)

print("Initial conditions:")
for _, row in scenario.iterrows():
    print(f"  {row['name']:<12}: pop={row.population:>7,}  "
          f"S={row.S:>6,}  E={row.E:>4}  I={row.I:>4,}  R={row.R:>5,}")


# ── Vital-dynamics rate arrays (shape: nticks × nnodes, per 1000/year) ───────
# LASER divides by 1000 internally. Passing daily per-capita is a silent error.

birthrate_array = np.full((NTICKS, NNODES), CBR, dtype=np.float32)
deathrate_array = np.full((NTICKS, NNODES), CDR, dtype=np.float32)

# Guard against the most common unit error
assert np.all(birthrate_array >= 1) and np.all(birthrate_array <= 60), \
    f"Birthrates must be per-1000/year; got {birthrate_array.min():.4f}"
assert np.all(deathrate_array >= 1) and np.all(deathrate_array <= 60), \
    f"Death rates must be per-1000/year; got {deathrate_array.min():.4f}"


# ── Pre-allocate capacity with calc_capacity ──────────────────────────────────
# Projects population growth over NTICKS at CBR=30 and applies a 2× safety
# factor so that LaserFrame.add() never exhausts free slots during the run.
# Model.__init__ performs this same call internally when birthrates are passed;
# we call it explicitly here to inspect the resulting per-patch slots.

capacity_per_patch = calc_capacity(birthrate_array, populations, safety_factor=2.0)
print(f"\ncalc_capacity pre-allocated slots (safety_factor=2.0):")
for name, cap in zip(PATCH_NAMES, capacity_per_patch):
    print(f"  {name:<12}: {cap:>10,}")
print(f"  Total         : {int(capacity_per_patch.sum()):>10,}")


# ── Age pyramid for BirthsByCBR ───────────────────────────────────────────────
# Exponential stable-age distribution over 0–99 years at the given CDR.

stable_age_wts = 1000.0 * np.exp(-(CDR / 1000.0) * np.arange(100))
pyramid = AliasedDistribution(stable_age_wts)


# ── Disease-duration distributions ───────────────────────────────────────────

expdurdist = dists.gamma(shape=EXP_SHAPE, scale=EXP_SCALE)   # latent:     mean 4 d
infdurdist = dists.gamma(shape=INF_SHAPE, scale=INF_SCALE)   # infectious: mean 10 d


# ── Model parameters ──────────────────────────────────────────────────────────
# Including gravity_k/a/b/c causes Model.__init__ to auto-build the coupling
# network from scenario centroids (Haversine distances + gravity formula).

params = PropertySet({
    "prng_seed":              42,
    "nticks":                 NTICKS,
    "beta":                   BETA,
    # Gravity network — auto-computed by Model.__init__ from scenario centroids
    "gravity_k":              0.02,
    "gravity_a":              0.0,    # source population exponent (LASER convention)
    "gravity_b":              1.0,    # destination population exponent
    "gravity_c":              2.0,    # distance decay exponent
    # capacity_safety_factor is forwarded to Model's internal calc_capacity call
    "capacity_safety_factor": 2.0,
})


# ── Build model ───────────────────────────────────────────────────────────────
# Passing birthrates triggers Model's internal calc_capacity call.
# Gravity params trigger automatic network construction from centroids.

model = Model(scenario, params, birthrates=birthrate_array)

# Verify the auto-built gravity network is non-trivial
assert model.network.sum() > 0, "Network is all zeros — check gravity params"
row_sums = model.network.sum(axis=1)
assert row_sums.max() < 0.3, \
    f"Network max row sum {row_sums.max():.3f} too high; agents may exceed capacity"

print(f"\nGravity network row sums (fraction of FOI exported per patch):")
for name, rs in zip(PATCH_NAMES, row_sums):
    print(f"  {name:<12}: {rs:.6f}")


# ── Assemble components in correct execution order ────────────────────────────
#
# Each tick, components run in list order:
#
#   1. Susceptible  — propagates S[t] → S[t+1] (bookkeeping)
#   2. Exposed      — ages exposed agents; fires E→I when etimer expires
#   3. Infectious   — ages infectious agents; fires I→R when itimer expires
#   4. Recovered    — propagates R[t] → R[t+1] (bookkeeping)
#   5. Transmission — computes spatial force of infection, draws S→E events,
#                     assigns etimer from expdurdist
#   6. BirthsByCBR  — Poisson births; calls on_birth on each component
#   7. MortalityByCDR — per-agent Bernoulli death trials; sets state = -1
#
# Susceptible/Recovered *wrap* the transition steps so that the invariant
# S + E + I + R = N holds at every tick.

model.components = [
    SEIR.Susceptible(model),                                          # 1
    SEIR.Exposed(model, expdurdist, infdurdist),                      # 2
    SEIR.Infectious(model, infdurdist),                               # 3
    SEIR.Recovered(model),                                            # 4
    SEIR.Transmission(model, expdurdist),                             # 5
    BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),  # 6
    MortalityByCDR(model, mortalityrates=deathrate_array),            # 7
]


# ── Run ───────────────────────────────────────────────────────────────────────

print(f"\nRunning {NTICKS}-tick ({NTICKS // 365}-year) Eastland SEIR simulation…")
print(f"  Agent capacity : {model.people.capacity:,}")
print(f"  R0=5, beta=0.5, latent=gamma(4,1), infectious=gamma(5,2)")
print(f"  CBR={CBR}/1000/yr, CDR={CDR}/1000/yr")

model.run("Eastland Respiratory SEIR")
print("Simulation complete.")


# ── Post-run verification ─────────────────────────────────────────────────────

print("\n─── Verification ────────────────────────────────────────────────────────")

S = model.nodes.S[:NTICKS]
E = model.nodes.E[:NTICKS]
I = model.nodes.I[:NTICKS]
R = model.nodes.R[:NTICKS]
N = S + E + I + R

# 1. Compartment non-negativity
assert np.all(S >= 0) and np.all(E >= 0) and np.all(I >= 0) and np.all(R >= 0), \
    "Negative compartment counts detected"
print("  Compartment non-negativity : PASS")

# 2. Population trajectory — must grow since CBR=30 > CDR=10
pop_0   = int(N[0].sum())
pop_end = int(N[-1].sum())
growth  = (pop_end / pop_0 - 1) * 100
status  = "PASS" if pop_end > pop_0 else "FAIL — check birth rate units!"
print(f"  Population growth          : {pop_0:,} → {pop_end:,} ({growth:+.1f}%)  {status}")

# 3. Epidemic dynamics
total_inf = int(model.nodes.newly_infected[:NTICKS].sum())
assert total_inf > 0, "No infections recorded — check beta and initial I counts"
print(f"  Total infections (10 yr)   : {total_inf:,}  PASS")

# 4. Spatial coupling — all patches should have cases
patch_totals   = model.nodes.newly_infected[:NTICKS].sum(axis=0)
patches_active = int(np.sum(patch_totals > 0))
status = "PASS" if patches_active == NNODES else "WARN — some patches had zero cases"
print(f"  Patches with cases         : {patches_active}/{NNODES}  {status}")

# 5. Compartment snapshot table
print(f"\n  {'Tick':<8} {'S':>10} {'E':>8} {'I':>8} {'R':>10} {'N':>10}")
print(f"  {'-' * 58}")
for tick in [0, NTICKS // 2, NTICKS - 1]:
    s = int(S[tick].sum())
    e = int(E[tick].sum())
    i = int(I[tick].sum())
    r = int(R[tick].sum())
    print(f"  {tick:<8} {s:>10,} {e:>8,} {i:>8,} {r:>10,} {s+e+i+r:>10,}")


# ── Diagnostic plots ──────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(14, 9))
colors = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12"]
days   = np.arange(NTICKS)

# (A) National epidemic curve (7-day rolling mean)
ax = axes[0, 0]
total_daily = model.nodes.newly_infected[:NTICKS].sum(axis=1).astype(np.float64)
kernel   = np.ones(7) / 7
smoothed = np.convolve(total_daily, kernel, mode="valid")
ax.plot(np.arange(len(smoothed)) / 365, smoothed, color="crimson", lw=1.2)
ax.set_xlabel("Year")
ax.set_ylabel("New infections (7-day avg)")
ax.set_title("(A) National Epidemic Curve")
ax.grid(True, alpha=0.3)

# (B) Susceptible fraction per patch
ax = axes[0, 1]
sample = np.arange(0, NTICKS, 7)
for i in range(NNODES):
    s = S[sample, i].astype(np.float64)
    n = np.maximum(N[sample, i].astype(np.float64), 1.0)
    ax.plot(sample / 365, s / n, color=colors[i], lw=1, label=PATCH_NAMES[i])
ax.axhline(1.0 / 5.0, color="black", ls=":", lw=1.5, alpha=0.6, label="1/R0 = 0.20")
ax.set_xlabel("Year")
ax.set_ylabel("S / N")
ax.set_title("(B) Susceptible Fraction per Patch")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (C) Population trajectory per patch
ax = axes[1, 0]
for i in range(NNODES):
    ax.plot(days / 365, N[:, i].astype(np.float64),
            color=colors[i], lw=1, label=PATCH_NAMES[i])
ax.set_xlabel("Year")
ax.set_ylabel("Population")
ax.set_title("(C) Population Trajectory (CBR=30, CDR=10 per 1000/yr)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (D) Weekly incidence heatmap (patches × time)
ax = axes[1, 1]
nweeks = NTICKS // 7
weekly = model.nodes.newly_infected[:nweeks * 7].reshape(nweeks, 7, NNODES).sum(axis=1)
im = ax.imshow(weekly.T, aspect="auto", cmap="hot",
               extent=[0, nweeks / 52, -0.5, NNODES - 0.5], origin="lower")
ax.set_yticks(range(NNODES))
ax.set_yticklabels(PATCH_NAMES)
ax.set_xlabel("Year")
ax.set_title("(D) Weekly Incidence Heatmap")
fig.colorbar(im, ax=ax, shrink=0.8, label="Cases / week")

plt.suptitle("Republic of Eastland — Respiratory SEIR  (R0=5, 10 yr)",
             fontsize=13, fontweight="bold")
plt.tight_layout()

outdir  = Path(__file__).parent / "eval" / "outputs"
outdir.mkdir(parents=True, exist_ok=True)
outpath = outdir / "eastland_respiratory_seir.png"
fig.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nPlot saved → {outpath}")
plt.close(fig)
```

---

**Key design decisions explained:**

**Beta derivation:** `R0 = beta × mean_infectious_period = beta × (5×2) = beta × 10`, so `beta = 0.5`.

**`calc_capacity` (explicit call, line ~80):** Called directly with `safety_factor=2.0` to inspect and log the per-patch slot allocations. `Model.__init__` performs the same call internally when `birthrates=` is passed — the `capacity_safety_factor=2.0` in `params` controls its internal call.

**Component order** — the 7-step sequence preserves `S+E+I+R=N` every tick:
| # | Component | Role |
|---|-----------|------|
| 1 | `Susceptible` | propagate S counts forward |
| 2 | `Exposed` | age E timers, fire E→I |
| 3 | `Infectious` | age I timers, fire I→R |
| 4 | `Recovered` | propagate R counts forward |
| 5 | `Transmission` | compute FOI, draw S→E, assign `etimer` |
| 6 | `BirthsByCBR` | Poisson births at CBR=30/1000/yr |
| 7 | `MortalityByCDR` | per-agent death trials at CDR=10/1000/yr |

**`MortalityByCDR` parameter name:** `mortalityrates=` (not `deathrates=`) — the skill flags this as a common silent error.

**Gravity network:** Including `gravity_k/a/b/c` in `params` triggers automatic network construction from scenario centroids inside `Model.__init__`. `gravity_a=0.0` follows the LASER convention that source population does not drive outward flow.
