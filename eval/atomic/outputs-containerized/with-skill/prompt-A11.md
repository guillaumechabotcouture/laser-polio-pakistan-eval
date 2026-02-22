Here is the complete code. The key unit fact is prominently documented in the header.

```python
"""
Birth process demo using BirthsByCBR.

BirthsByCBR unit convention
────────────────────────────
The `birthrates` parameter must be in per-1,000 population per year.
  ✓ Correct:   44.0           (44 births per 1,000 per year)
  ✗ Wrong:     0.044          (fraction)
  ✗ Wrong:     44/365/1000    (daily per-capita ≈ 0.000121)

The framework converts internally each tick:
  p_birth_per_tick = (1 + CBR/1000)^(1/365) - 1

This same per-1000/year value must also be passed to Model() at
construction so calc_capacity() can pre-allocate enough agent slots.
If you omit it or pass the wrong units, capacity ≈ initial population,
LaserFrame.add() returns no free slots, and no births occur — silently.
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


# ── 1. Simulation parameters ──────────────────────────────────────────────────
CBR    = 44.0        # crude birth rate  — per 1,000 population per year
CDR    = 10.0        # crude death rate  — per 1,000 population per year
NTICKS = 5 * 365     # 5-year simulation

params = PropertySet({
    "prng_seed":             42,
    "nticks":                NTICKS,
    "beta":                  0.0,    # no transmission — isolates birth dynamics
    "exp_shape":             40,
    "exp_scale":             0.25,
    "inf_mean":              8.0,
    "inf_sigma":             2.0,
    "capacity_safety_factor": 3.0,  # extra headroom for population growth
})


# ── 2. Scenario (3 patches) ───────────────────────────────────────────────────
pops   = np.array([100_000, 50_000, 25_000], dtype=np.int64)
nnodes = len(pops)

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(nnodes),
        "name":       ["PatchA", "PatchB", "PatchC"],
        "population": pops,
        "S":          pops.copy(),
        "E":          np.zeros(nnodes, dtype=np.int64),
        "I":          np.zeros(nnodes, dtype=np.int64),
        "R":          np.zeros(nnodes, dtype=np.int64),
        "geometry":   [Point(73.0, 33.0), Point(74.0, 34.0), Point(75.0, 35.0)],
    },
    crs="EPSG:4326",
)

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "Initial compartments must sum to population in every patch"


# ── 3. Birth rate array ───────────────────────────────────────────────────────
# ValuesMap.from_scalar produces shape (NTICKS, nnodes) filled with CBR = 44.0.
# Every element must be in per-1,000/year — the assertion below catches wrong units.
birthrate_map = ValuesMap.from_scalar(CBR, NTICKS, nnodes)

assert np.all(birthrate_map.values >= 1) and np.all(birthrate_map.values <= 60), (
    f"Birthrates must be per-1000/year (plausible range 1–60). "
    f"Got min={birthrate_map.values.min():.6f}. "
    f"Do not convert to daily fractions before passing."
)


# ── 4. Age pyramid for newborns ───────────────────────────────────────────────
# BirthsByCBR samples ages for newly created agents from this distribution.
# Place all probability mass at index 0 so every newborn starts at age 0 days.
newborn_ages      = np.zeros(365 * 90, dtype=float)
newborn_ages[0]   = 1.0
pyramid           = AliasedDistribution(newborn_ages)


# ── 5. Death rate array ───────────────────────────────────────────────────────
# MortalityByCDR also expects per-1,000/year.
# NOTE: parameter keyword is `mortalityrates=`, not `deathrates=`.
deathrate_map = ValuesMap.from_scalar(CDR, NTICKS, nnodes)


# ── 6. Build model ────────────────────────────────────────────────────────────
# Pass birthrates to Model() so calc_capacity() pre-allocates enough agent slots.
# Use the same per-1000/year array that BirthsByCBR will consume each tick.
model = Model(scenario, params, birthrates=birthrate_map.values)


# ── 7. Duration distributions ─────────────────────────────────────────────────
expdurdist = dists.gamma(shape=params.exp_shape, scale=params.exp_scale)
infdurdist = dists.normal(loc=params.inf_mean,  scale=params.inf_sigma)


# ── 8. Component list ─────────────────────────────────────────────────────────
# BirthsByCBR placement: after disease components so newborns (state=SUSCEPTIBLE)
# are not exposed to transmission in the same tick they are created.
#
# BirthsByCBR calls on_birth(istart, iend, tick) on any other component that
# defines that method, allowing custom per-agent initialization at birth.
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist),
    BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=deathrate_map.values),
]


# ── 9. Run ────────────────────────────────────────────────────────────────────
initial_pop = int(pops.sum())
print(f"Initial population: {initial_pop:,}")

model.run("BirthDemo")


# ── 10. Verify births occurred ────────────────────────────────────────────────
N = (
    model.nodes.S[:NTICKS]
    + model.nodes.E[:NTICKS]
    + model.nodes.I[:NTICKS]
    + model.nodes.R[:NTICKS]
)

final_pop      = int(N[-1].sum())
net_rate       = (CBR - CDR) / 1000            # 0.034 → 3.4 %/yr
expected_final = initial_pop * (1 + net_rate) ** 5

print(f"Final population:   {final_pop:,}")
print(f"Expected (~{net_rate*100:.1f}%/yr × 5 yr): {expected_final:,.0f}")
print(f"Population change:  {final_pop - initial_pop:+,}")

assert final_pop > initial_pop, (
    "Population did not grow. Most likely cause: birthrates passed in wrong units. "
    "BirthsByCBR requires per-1000/year (e.g. 44.0), not a daily fraction."
)
print("PASSED: births occurred and population grew.")
```

---

**Unit summary for `BirthsByCBR`:**

| What to pass | Example value | What NOT to pass |
|---|---|---|
| Per-1,000 population **per year** | `44.0` | `0.044` (fraction) or `0.000121` (daily per-capita) |

The framework does the internal conversion: `p_birth_per_tick = (1 + CBR/1000)^(1/365) - 1`. The same array — same units — must go to both `Model(birthrates=...)` and `BirthsByCBR(birthrates=...)`: the former is for capacity pre-allocation, the latter drives the per-tick birth draws.
