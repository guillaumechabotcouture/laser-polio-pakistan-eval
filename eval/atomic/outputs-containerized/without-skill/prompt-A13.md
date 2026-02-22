Now I have everything I need. Let me write the complete code.

```python
import numpy as np
from laser_measles.utils import calc_capacity
from laser_core.laserframe import LaserFrame

# ── Simulation parameters ────────────────────────────────────────────────────
INITIAL_POPULATION = 500_000   # agents present at t = 0
CBR                = 44.0      # crude birth rate: births per 1,000 people per year
YEARS              = 10
NTICKS             = YEARS * 365  # calc_capacity works in days (ticks), not years

# ── Pre-allocate storage ─────────────────────────────────────────────────────
# calc_capacity projects the maximum population size the simulation can reach,
# assuming births accumulate continuously at the given CBR.
# A 1 % safety buffer is standard practice to cover stochastic over-runs.
raw_capacity = calc_capacity(
    np.uint32(INITIAL_POPULATION),
    np.uint32(NTICKS),
    np.float32(CBR),
    verbose=True,               # prints projected growth to stdout
)
capacity = int(1.01 * raw_capacity)

print(f"Initial population  : {INITIAL_POPULATION:>10,}")
print(f"Simulation length   : {YEARS} years ({NTICKS:,} daily ticks)")
print(f"CBR                 : {CBR} per 1,000 per year")
print(f"Projected max pop   : {raw_capacity:>10,}")
print(f"Pre-alloc capacity  : {capacity:>10,}  (projected × 1.01 buffer)")

# ── Build the agent store ────────────────────────────────────────────────────
# LaserFrame allocates all agent arrays up front at `capacity` slots.
# Only `INITIAL_POPULATION` slots are "active" at t = 0; the rest sit
# dormant and are activated one-by-one as births occur during the run.
population = LaserFrame(capacity)
population.add(INITIAL_POPULATION)  # activate the starting cohort

# Attach per-agent properties (arrays are sized to `capacity`, not `count`)
population.add_scalar_property("date_of_birth", dtype=np.int32,  default=-1)
population.add_scalar_property("date_of_death", dtype=np.int32,  default=-1)
population.add_scalar_property("alive",         dtype=np.bool_,  default=True)
population.add_scalar_property("nodeid",        dtype=np.uint16, default=0)

print(f"\nLaserFrame ready:")
print(f"  .capacity = {population.capacity:>10,}  — total pre-allocated slots")
print(f"  .count    = {population.count:>10,}  — agents active at t = 0")
print(f"  headroom  = {population.capacity - population.count:>10,}  — slots reserved for future births")
```

---

## What `calc_capacity` expects — unit reference

| Parameter | Type | Units |
|---|---|---|
| `population` | `np.uint32` | Agents (integer count) |
| `nticks` | `np.uint32` | **Days** — convert years × 365 before passing |
| `cbr` | `np.float32` | **Births per 1,000 people per year** — pass 44, not 0.044 |
| `verbose` | `bool` | Optional; prints projected growth |

The key thing to get right is the **day-based tick count**. Internally the function converts:

```
daily_rate = (cbr / 1000) / 365.0
capacity   = population × (1 + daily_rate) ** nticks
```

So `nticks` must be in days (3,650 for a 10-year run), and `cbr` stays in the conventional epidemiological unit of *per 1,000 per year* — no unit conversion needed for CBR itself.

For 500,000 people, CBR = 44, over 10 years the projected capacity works out to roughly **776,000**, and with the 1 % buffer you pre-allocate about **784,000** slots. All agent property arrays (dates, node IDs, etc.) are allocated at that full size immediately, so no reallocation ever happens during the simulation loop.

Sources:
- [laser_measles.model — laser_measles 0.5.0 documentation](https://docs.idmod.org/projects/laser-measles/en/latest/_modules/laser_measles/model.html)
- [Vital Dynamics Model — LASER 0.4.0 documentation](https://docs.idmod.org/projects/laser/en/v0.4.1/vdexample.html)
- [Software Architecture & Design — LASER 0.6.0 documentation](https://docs.idmod.org/projects/laser/en/latest/architecture.html)
