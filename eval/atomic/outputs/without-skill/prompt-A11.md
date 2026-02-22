Here's the complete code adding a birth process using `BirthsByCBR`:

```python
import numpy as np
from laser.generic import Model
from laser.generic.vitaldynamics import BirthsByCBR
from laser.generic.utils import ValuesMap
from laser.core.demographics import AliasedDistribution
from laser.core import PropertySet

# --- Configuration ---
NTICKS = 365          # 1 year of daily ticks
NNODES = 1            # single node/patch
INITIAL_POP = 100_000

# Crude birth rate: 44 per 1,000 population per year
CBR = 44.0

# --- Birth rate setup ---
# BirthsByCBR expects CBR in units of: births per 1,000 population per year.
# Internally it converts to a daily per-capita rate via:
#   daily_rate = (1 + CBR/1000)^(1/365) - 1
# Then draws: births = Poisson(daily_rate * N)
birthrates = ValuesMap.from_scalar(CBR, NTICKS, NNODES)

# --- Age pyramid (uniform placeholder) ---
# In practice, load real demographic data, e.g.:
#   age_data = np.loadtxt("age_distribution.csv", delimiter=",", usecols=0)
age_data = np.ones(89)  # flat distribution across 89 age cohorts (0-88)
pyramid = AliasedDistribution(age_data)

# --- Build and run the model ---
params = PropertySet({"nticks": NTICKS})
model = Model(params=params)

# Initialize population across nodes
model.nodes.population = np.array([INITIAL_POP])

# Add BirthsByCBR as the LAST component (so end-of-tick counts are correct)
model.components = [
    # ... other components (e.g., SEIR.Susceptible, Transmission, etc.) go here ...
    BirthsByCBR(model, birthrates, pyramid),
]

model.run()

# --- Inspect results ---
total_births = model.nodes.births[:NTICKS].sum()
print(f"CBR: {CBR} per 1,000 per year")
print(f"Population: {INITIAL_POP:,}")
print(f"Simulated days: {NTICKS}")
print(f"Total births: {total_births:,}")
print(f"Expected births (approx): {int(INITIAL_POP * CBR / 1000 * NTICKS / 365):,}")
```

## Units for `BirthsByCBR`

**The `birthrates` parameter expects the crude birth rate in births per 1,000 population per year** — the standard demographic convention.

So for a CBR of 44 per 1,000/year, you pass `44.0` directly. Internally, the class converts this to a daily per-capita probability:

```
daily_rate = (1 + 44/1000)^(1/365) - 1 ≈ 0.000118
```

Then at each daily tick it draws births from a Poisson distribution: `Poisson(daily_rate * N)`.

### Key notes

- **Component ordering**: `BirthsByCBR` should be the **last** component so population counts at the end of the tick are accurate.
- **Newborns enter as Susceptible** and trigger `on_birth()` callbacks on all other components.
- **Spatial variation**: Use `ValuesMap.from_nodes()` or `ValuesMap.from_array()` instead of `from_scalar()` if different nodes have different CBRs.
