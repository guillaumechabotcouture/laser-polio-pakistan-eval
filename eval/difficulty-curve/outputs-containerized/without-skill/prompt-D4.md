Here is the complete code at `/tmp/eastland_seir.py`. Let me walk through every design decision:

---

**`eastland_seir.py`** — full annotated listing:

```python
# ── Imports ───────────────────────────────────────────────────────────────────
from laser.core import PropertySet
from laser.core.demographics import AliasedDistribution    # Walker-alias sampler
from laser.core.distributions import gamma as gamma_dist   # Numba-wrapped Γ dist
from laser.core.utils import calc_capacity, grid, initialize_population
from laser.generic import Model
from laser.generic.SEIR import (
    Exposed, Infectious, Recovered, Susceptible, Transmission,
)
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
```

---

**Key design decisions:**

| Decision | Detail |
|---|---|
| **β calculation** | `β = R0 / mean_I = 5 / (5×2) = 0.5 day⁻¹` |
| **Latent dist** | `gamma_dist(4, 1)` → Gamma(shape=4, scale=1), mean=4 d |
| **Infectious dist** | `gamma_dist(5, 2)` → Gamma(shape=5, scale=2), mean=10 d |
| **Scenario** | `grid(M=2, N=2, population_fn=pop_fn)` → 2×2 GeoDataFrame with real geometry for spatial coupling |
| **Initial conditions** | `initialize_population(scenario, [0.90, 0.00, 0.01, 0.09])` applied per patch |
| **Capacity pre-allocation** | `birthrates` array passed to `Model()` → triggers `calc_capacity(birthrates, population)` internally, sizing the `LaserFrame` for all births over 10 yr |
| **Birthrates shape** | `np.full((3650, 4), 30.0)` — (nticks × nnodes), CBR=30/1000/yr |
| **Pyramid** | `AliasedDistribution(np.ones(20))` — uniform across 20 five-year bins; `track=False` since no age-based interventions |
| **Component order** | `Susceptible → Transmission → Exposed → Infectious → Recovered → BirthsByCBR → MortalityByCDR` |

**Component ordering rationale:**  
`Transmission` (2nd) reads `I[tick+1]` which `_initialize_flows` has already set to `I[tick]`, so it sees the correct tick-start infectious count before any recoveries are processed. `Exposed` (3rd) then handles E→I, `Infectious` (4th) handles I→R. Demographics come last so they act on the fully resolved disease state.
