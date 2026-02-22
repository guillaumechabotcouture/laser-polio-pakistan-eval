Now I have all the information needed. Here is the complete code example with explicit documentation of units.

---

## `BirthsByCBR` in the LASER Framework

**Package location:** `laser.generic.vitaldynamics` (part of the `laser-generic` package, built on top of `laser-core`)

**CBR parameter units:** The `birthrates` array must contain values in **births per 1,000 people per year**. So a CBR of 44 means pass `44.0`. The component converts this internally to a per-day probability using:

```
daily_rate = (1 + CBR/1000)^(1/365) - 1
births_today = Poisson(N × daily_rate)
```

---

```python
import numpy as np
import pandas as pd

from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution
from laser.generic import Model, SEIR
from laser.generic.vitaldynamics import BirthsByCBR
from laser.generic.utils import ValuesMap

# ---------------------------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------------------------
NTICKS = 365 * 10   # 10 years of daily time steps

# ---------------------------------------------------------------------------
# Scenario: spatial patches (nodes) with initial populations
# ---------------------------------------------------------------------------
scenario = pd.DataFrame({
    "name":       ["Node_A", "Node_B"],
    "population": [500_000,  200_000],
    "latitude":   [9.0,      12.0],
    "longitude":  [7.5,      13.5],
})

# Initialise compartments: everyone starts susceptible
scenario["S"] = scenario["population"]
scenario["E"] = 0
scenario["I"] = 0
scenario["R"] = 0

# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------
params = PropertySet({
    "nticks": NTICKS,
    "seed":   42,
})

# ---------------------------------------------------------------------------
# Age pyramid (required by BirthsByCBR to assign dates of birth to newborns)
#
# A real model should load empirical age-distribution data, e.g.:
#   age_weights = np.loadtxt("Nigeria-AgeDistribution-2020.csv", delimiter=",")
#
# Here a synthetic uniform pyramid is used for illustration.
# ---------------------------------------------------------------------------
age_weights = np.ones(89)          # cohorts for ages 0–88
pyramid = AliasedDistribution(age_weights)

# ---------------------------------------------------------------------------
# Birth rates
#
# BirthsByCBR expects the crude birth rate (CBR) in units of:
#
#     births per 1,000 people per year
#
# A CBR of 44 means 44 live births per 1,000 population per year.
# The component converts this to a daily probability internally:
#
#     daily_rate  =  (1 + CBR/1000)^(1/365) - 1
#     births_today = Poisson(N_node × daily_rate)
#
# ValuesMap.from_scalar(value, nticks, nnodes) broadcasts a single scalar
# into a constant [nticks × nnodes] rate array.
# ---------------------------------------------------------------------------
N_NODES = len(scenario)
CBR = 44.0   # births per 1,000 people per year  <-- units are explicit here

birthrates = ValuesMap.from_scalar(CBR, NTICKS, N_NODES)

# ---------------------------------------------------------------------------
# Build model and attach components
#
# IMPORTANT: BirthsByCBR must be placed LAST in model.components so that
# end-of-tick population counts are fully updated before births are computed.
# ---------------------------------------------------------------------------
model = Model(scenario, params, birthrates=birthrates)

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model),
    SEIR.Infectious(model),
    SEIR.Recovered(model),
    SEIR.Transmission(model),
    BirthsByCBR(           # ← must be last
        model,
        birthrates,
        pyramid,
        track=True,        # records dob on each agent; required by MortalityByEstimator
    ),
]

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
model.run()

# ---------------------------------------------------------------------------
# Inspect results
# ---------------------------------------------------------------------------
total_births = model.nodes.births[:NTICKS].sum()
print(f"Total births over simulation:    {total_births:,}")
print(f"Final S (susceptible) per node:  {model.nodes.S[NTICKS]}")
print(f"Dates of birth (first 5 agents): {model.people.dob[:5]}")
```

---

### Key points

| Detail | Value |
|---|---|
| Import path | `from laser.generic.vitaldynamics import BirthsByCBR` |
| CBR units | **births per 1,000 people per year** |
| Value to pass for CBR = 44 | `44.0` (not 0.044, not 44/365) |
| Internal daily conversion | `(1 + 44/1000)^(1/365) - 1 ≈ 0.0001178` per person per day |
| Placement in `model.components` | **Last** — uses end-of-tick population counts |
| `track=True` | Adds a `dob` property to each agent; required if using `MortalityByEstimator` |

For spatially varying rates (different CBR per node), use `ValuesMap.from_nodes(cbr_array, nticks)`. For time-varying rates, use `ValuesMap.from_timeseries(cbr_timeseries, nnodes)`.

Sources:
- [LASER GitHub — InstituteforDiseaseModeling/laser](https://github.com/InstituteforDiseaseModeling/laser)
- [Vital Dynamics Model — LASER 0.4.0 documentation](https://docs.idmod.org/projects/laser/en/v0.4.1/vdexample.html)
- [laser-core on PyPI](https://pypi.org/project/laser-core/)
- [laser-generic on PyPI](https://pypi.org/project/laser-generic/)
