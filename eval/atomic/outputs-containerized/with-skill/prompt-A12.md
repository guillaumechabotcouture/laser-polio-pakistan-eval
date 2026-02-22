Here is a complete, self-contained example. The critical unit fact is front and center.

**`MortalityByCDR` expects rates in per-1,000/year.** The framework converts internally each tick using: `p_death = 1 - (1 - CDR/1000)^(1/365)`. Do **not** pass a daily per-capita fraction (e.g., 0.0000356) — that is the most common silent failure.

```python
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR

# ---------------------------------------------------------------------------
# 1. Scenario: 3 patches
# ---------------------------------------------------------------------------
pops = np.array([500_000, 200_000, 80_000], dtype=np.int32)
lats = np.array([30.0, 31.5, 32.0])
lons = np.array([70.0, 71.0, 71.5])

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(len(pops)),
        "name":       ["PatchA", "PatchB", "PatchC"],
        "population": pops,
        # Seed 3 infectious, seed 70% recovered, rest susceptible
        "S": (pops * 0.30).astype(np.int32),
        "E": np.zeros(len(pops), dtype=np.int32),
        "I": np.full(len(pops), 3, dtype=np.int32),
        "R": (pops * 0.70 - 3).astype(np.int32),
        "geometry": [Point(lon, lat) for lat, lon in zip(lats, lons)],
    },
    crs="EPSG:4326",
)

assert (scenario.S + scenario.E + scenario.I + scenario.R == scenario.population).all(), \
    "S+E+I+R must equal population in every patch"

# ---------------------------------------------------------------------------
# 2. Parameters
# ---------------------------------------------------------------------------
nticks = 10 * 365   # 10-year simulation
nnodes = len(scenario)

parameters = PropertySet({
    "prng_seed":              42,
    "nticks":                 nticks,
    # Disease
    "beta":                   14.0,
    "exp_shape":              40,
    "exp_scale":              0.25,   # gamma exposed duration ~10 days
    "inf_mean":               8,
    "inf_sigma":              2,
    # Vital dynamics (per-1000/year)
    "cbr":                    30.0,   # crude birth rate
    "cdr":                    13.0,   # crude death rate  <-- 13 per 1000/year
    # Gravity network
    "gravity_k":              0.01,
    "gravity_b":              0.5,
    "gravity_c":              1.5,
    # Capacity
    "capacity_safety_factor": 3.0,
})

# ---------------------------------------------------------------------------
# 3. Build mortalityrates array (nticks × nnodes) in per-1000/year
#
#    MortalityByCDR EXPECTS: per-1000/year (e.g., 13.0)
#    Do NOT pass:            daily per-capita (e.g., 0.0000356)
#                            annual fraction  (e.g., 0.013)
#
#    The component converts internally each tick:
#      p_death = 1 - (1 - CDR/1000)^(1/365)
# ---------------------------------------------------------------------------
CDR_PER_1000_PER_YEAR = 13.0   # <-- this is the value you pass

# Uniform CDR across all patches and all ticks
mortalityrates = ValuesMap.from_scalar(CDR_PER_1000_PER_YEAR, nticks, nnodes)

# Unit sanity check — catches the "passed a fraction by mistake" failure
assert np.all(mortalityrates.values >= 1) and np.all(mortalityrates.values <= 60), (
    f"mortalityrates must be per-1000/year (expected 1–60), "
    f"got min={mortalityrates.values.min():.6f}"
)

# ---------------------------------------------------------------------------
# 4. Birth rates (also per-1000/year — same convention as mortality)
# ---------------------------------------------------------------------------
CBR_PER_1000_PER_YEAR = 30.0
birthrates = ValuesMap.from_scalar(CBR_PER_1000_PER_YEAR, nticks, nnodes)

assert np.all(birthrates.values >= 1) and np.all(birthrates.values <= 60), \
    f"birthrates must be per-1000/year, got {birthrates.values.min():.6f}"

# ---------------------------------------------------------------------------
# 5. Age pyramid for births (stable exponential)
# ---------------------------------------------------------------------------
stable_ages = np.exp(-0.02 * np.arange(90))
pyramid = AliasedDistribution(stable_ages)

# ---------------------------------------------------------------------------
# 6. Disease duration distributions
# ---------------------------------------------------------------------------
expdurdist = dists.gamma(shape=parameters.exp_shape, scale=parameters.exp_scale)
infdurdist = dists.normal(loc=parameters.inf_mean, scale=parameters.inf_sigma)

# ---------------------------------------------------------------------------
# 7. Build model — birthrates passed here for capacity pre-allocation
# ---------------------------------------------------------------------------
model = Model(scenario, parameters, birthrates=birthrates.values)

# ---------------------------------------------------------------------------
# 8. Assemble components
#    MortalityByCDR keyword is `mortalityrates=` (NOT `deathrates=`)
# ---------------------------------------------------------------------------
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist),
    BirthsByCBR(model, birthrates=birthrates.values, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=mortalityrates.values),   # <-- CDR in per-1000/year
]

# ---------------------------------------------------------------------------
# 9. Run
# ---------------------------------------------------------------------------
model.run("CDR Example")

# ---------------------------------------------------------------------------
# 10. Quick verification: population should not be static
# ---------------------------------------------------------------------------
initial_pop = pops.sum()
final_pop = (
    model.nodes.S[nticks - 1].sum() +
    model.nodes.E[nticks - 1].sum() +
    model.nodes.I[nticks - 1].sum() +
    model.nodes.R[nticks - 1].sum()
)
print(f"Initial population : {initial_pop:,}")
print(f"Final population   : {final_pop:,}")
print(f"Net change         : {final_pop - initial_pop:+,}")

# With CBR=30 > CDR=13, expect modest population growth
assert final_pop != initial_pop, \
    "Population is static — check that rates are per-1000/year, not fractions"
```

---

**Key facts about `MortalityByCDR`:**

| | Detail |
|---|---|
| **Keyword** | `mortalityrates=` — not `deathrates=` |
| **Units** | Per-1,000/year. Pass `13.0`, not `0.013` or `0.0000356` |
| **Shape** | `(nticks, nnodes)` ndarray — use `ValuesMap.from_scalar()` for uniform rates |
| **Internal conversion** | `p_death = 1 - (1 - CDR/1000)^(1/365)` per tick per agent |
| **On death** | Sets `state = -1` (DECEASED) and decrements the agent's compartment count in `model.nodes` |
| **Valid range** | 1–60 per 1000/year for the assertion; typical CDR is 5–20 |
