Here is the complete simulation code written to `/tmp/polio_seir_pakistan.py`:

---

```python
"""
Spatial SEIR model for polio transmission dynamics across 10 Pakistan districts.
Bill & Melinda Gates Foundation – Global Health Division
LASER framework (laser-generic package).

Epidemiological configuration
------------------------------
Patches           : 10 districts, 100,000 population each
R0                : ≈ 6
Latent period     : 3 days   (sigma = 1/3 day⁻¹)
Infectious period : 28 days  (gamma = 1/28 day⁻¹)
beta              : R0 * gamma = 6/28 ≈ 0.2143 day⁻¹
Initial state     : 95 % recovered (immune), 5 infectious per patch
Simulation        : 10 years (3,650 daily ticks)
"""

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

from laser.core import PropertySet
import laser.core.distributions as dists
from laser.generic import SEIR, Model

# ── Epidemiological parameters ────────────────────────────────────────────────
N_PATCHES      = 10
POP_PER_PATCH  = 100_000
R0             = 6.0
LATENT_DAYS    = 3        # 1 / sigma
INFECT_DAYS    = 28       # 1 / gamma

gamma = 1.0 / INFECT_DAYS          # 0.03571 day⁻¹
beta  = R0 * gamma                  # 0.21429 day⁻¹  (frequency-dependent)

# Per-patch initial compartments
INIT_I = 5
INIT_R = int(POP_PER_PATCH * 0.95)          # 95,000 recovered / immune
INIT_E = 0
INIT_S = POP_PER_PATCH - INIT_R - INIT_I   # 4,995 susceptible

assert INIT_S + INIT_E + INIT_I + INIT_R == POP_PER_PATCH

# Simulation length
SIM_YEARS = 10
NTICKS    = SIM_YEARS * 365         # 3,650 daily ticks

# ── Build scenario GeoDataFrame ───────────────────────────────────────────────
# 10 patches in a 2×5 grid representing Pakistan districts.
# WGS-84 (EPSG:4326); each cell ≈ 1° × 1° (~100 km).
CELL_DEG  = 1.0
ORIGIN_X  = 65.0    # ~65°E (western Pakistan)
ORIGIN_Y  = 30.0    # ~30°N
GRID_COLS = 5

records = []
for patch_id in range(N_PATCHES):
    row, col = divmod(patch_id, GRID_COLS)
    x0, y0 = ORIGIN_X + col * CELL_DEG, ORIGIN_Y + row * CELL_DEG
    x1, y1 = x0 + CELL_DEG, y0 + CELL_DEG
    geom = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
    records.append({
        "nodeid":     patch_id,
        "population": POP_PER_PATCH,
        "geometry":   geom,
        "S": INIT_S,
        "E": INIT_E,
        "I": INIT_I,
        "R": INIT_R,
    })

scenario = gpd.GeoDataFrame(records, crs="EPSG:4326")

# ── Duration distributions ────────────────────────────────────────────────────
# Latent duration: Gamma(mean = 3 days)  →  shape=3, scale=1
exp_shape = 3.0
exp_scale = LATENT_DAYS / exp_shape     # 1.0 day
expdist   = dists.gamma(shape=exp_shape, scale=exp_scale)

# Infectious duration: Normal(mean = 28 days, std = 3 days)
infdist = dists.normal(loc=float(INFECT_DAYS), scale=3.0)

# ── Model configuration ───────────────────────────────────────────────────────
params = PropertySet({"nticks": NTICKS, "beta": beta})
model  = Model(scenario, params)

# Component pipeline (order: state trackers → Transmission)
#   Exposed      – decrements etimer; E→I, assigns itimer from infdist
#   Infectious   – decrements itimer; I→R (permanent immunity)
#   Transmission – computes FOI = beta*(I/N), draws S→E, assigns etimer from expdist;
#                  spatial coupling via internal gravity network
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdist, infdist),
    SEIR.Infectious(model, infdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdist),
]

# ── Run ───────────────────────────────────────────────────────────────────────
print("=" * 65)
print("Polio SEIR – 10 Pakistan districts, LASER framework")
print("=" * 65)
print(f"  R0 = {R0},  beta = {beta:.5f} day⁻¹")
print(f"  Latent {LATENT_DAYS} d, Infectious {INFECT_DAYS} d")
print(f"  Per patch: S={INIT_S:,}  E={INIT_E}  I={INIT_I}  R={INIT_R:,}")
print(f"  Duration: {SIM_YEARS} years ({NTICKS:,} ticks)")
print("-" * 65)

model.run("Polio Pakistan Districts")

# ── Results ───────────────────────────────────────────────────────────────────
# model.nodes.{S,E,I,R} → shape (T, N_PATCHES)
S_t   = model.nodes.S.sum(axis=1)
E_t   = model.nodes.E.sum(axis=1)
I_t   = model.nodes.I.sum(axis=1)
R_t   = model.nodes.R.sum(axis=1)
N_t   = S_t + E_t + I_t + R_t
years = np.arange(len(S_t)) / 365.0

assert abs(N_t[-1] - N_t[0]) < 1.0, "Population not conserved"

peak_day = int(np.argmax(I_t))
print(f"Peak infectious: {I_t[peak_day]:,.0f}  on day {peak_day} (year {peak_day/365:.2f})")
print(f"Final:  S={S_t[-1]:,.0f}  E={E_t[-1]:,.0f}  "
      f"I={I_t[-1]:,.0f}  R={R_t[-1]:,.0f}")

for p in range(N_PATCHES):
    pk = int(np.argmax(model.nodes.I[:, p]))
    print(f"  District {p:2d}: peak {model.nodes.I[pk,p]:6,.0f}  day {pk:4d}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

ax = axes[0]
ax.plot(years, S_t / 1e3, label="S – Susceptible",        linewidth=1.8)
ax.plot(years, E_t / 1e3, label="E – Exposed (latent)",   linewidth=1.8)
ax.plot(years, I_t / 1e3, label="I – Infectious",         linewidth=1.8)
ax.plot(years, R_t / 1e3, label="R – Recovered / Immune", linewidth=1.8)
ax.axvline(peak_day / 365, color="grey", linestyle="--",
           linewidth=0.9, label=f"Peak I (day {peak_day})")
ax.set_ylabel("Population (thousands)")
ax.set_title(f"Aggregate SEIR – 10 Pakistan districts  (polio, R0={R0}, {SIM_YEARS} yr)")
ax.legend(loc="center right", fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
colors = plt.cm.tab10(np.linspace(0, 1, N_PATCHES))
for p in range(N_PATCHES):
    ax.plot(years, model.nodes.I[:, p],
            color=colors[p], linewidth=0.9, alpha=0.85, label=f"District {p}")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Infectious count")
ax.set_title("Per-district infectious prevalence")
ax.legend(fontsize=7, ncol=5, loc="upper right")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("polio_seir_pakistan.png", dpi=150, bbox_inches="tight")
print("Figure saved -> polio_seir_pakistan.png")
plt.show()
```

---

**Design decisions and their rationale:**

| Parameter | Value | Derivation |
|---|---|---|
| `beta` | 6/28 ≈ 0.2143 day⁻¹ | Frequency-dependent: β = R0 × γ |
| `expdist` | Gamma(shape=3, scale=1) | Mean = 3 days; shape=3 gives a realistic right-skewed incubation distribution |
| `infdist` | Normal(mean=28, std=3) | 28-day infectious period; std=3 avoids near-zero draws |
| Initial immune fraction | 95 % | Below the 83.3 % herd-immunity threshold, so an outbreak occurs |
| Spatial layout | 2×5 grid in WGS-84 | ~65–70°E, 30–31°N (Pakistan); LASER derives a gravity coupling network from centroids |

**Component pipeline order** matters: `Susceptible → Exposed → Infectious → Recovered → Transmission`. The `Transmission` component runs last so the force of infection computed at tick `t` uses the infectious count from the same tick.

**Post-run access:** `model.nodes.S/E/I/R` are 2D arrays of shape `(T, N_PATCHES)`. Summing `axis=1` gives the aggregate time-series; indexing `[:, p]` gives the per-district curve.
