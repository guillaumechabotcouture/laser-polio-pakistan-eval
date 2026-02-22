Here's the complete code for adding births via `BirthsByCBR`:

```python
import numpy as np
from laser.generic.vitaldynamics import BirthsByCBR
from laser.generic.utils import ValuesMap

# --- Units ---
# BirthsByCBR expects the crude birth rate in **per-1000 population per year**.
# It divides by 1000 internally to convert to a per-capita rate, then
# divides by 365 to get a daily probability.
#
# A CBR of 44 means: 44 births per 1000 people per year.
# Do NOT convert to daily per-capita (e.g., 44/1000/365 = 0.00012) — that
# would be silently treated as 0.00012 per 1000/year, producing zero births.

cbr = 44  # per-1000/year — pass this value directly

# Build the birthrates array: shape (nticks, nnodes)
# ValuesMap.from_scalar broadcasts a single value across all ticks and nodes.
nticks = parameters.nticks       # e.g., 20 * 365
nnodes = len(scenario)           # number of spatial patches
birthrate_map = ValuesMap.from_scalar(cbr, nticks, nnodes)

# Validate units before passing to LASER
assert np.all(birthrate_map.values >= 1) and np.all(birthrate_map.values <= 60), \
    f"Birthrates must be per-1000/year, got {birthrate_map.values.min():.4f}–{birthrate_map.values.max():.4f}"

# --- Model construction ---
# IMPORTANT: pass birthrates to Model.__init__ so calc_capacity() can
# pre-allocate enough agent slots for population growth.
model = Model(scenario, parameters, birthrates=birthrate_map.values)

# --- Add BirthsByCBR as a component ---
model.components = [
    # ... other components (Susceptible, Transmission, etc.) ...
    BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
    # ... other components ...
]
```

**The critical unit point:** `BirthsByCBR` expects the birth rate as **per-1000 population per year** — a value like `44`, not `0.00012`. Passing a pre-converted daily per-capita rate is the #1 silent failure: `calc_capacity` sees near-zero growth, allocates no extra agent slots, and births silently never happen. The model runs fine but the population stays flat.

The same value must also be passed to `Model(...)` via `birthrates=` so that `calc_capacity()` correctly pre-allocates space for the newborn agents over the simulation.
