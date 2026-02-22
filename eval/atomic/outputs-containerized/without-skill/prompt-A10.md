The filesystem is read-only, so I'll present the answer directly.

---

## SEIR + Transmission in a LASER Model

### Component registration

```python
from laser_core.model import Model
from laser_core.propertyset import PropertySet

# SEIR state-tracking / progression components
from laser_measles.seir import (
    Susceptible,   # S pool – births, waning immunity
    Exposed,       # E pool – incubation timer countdown, E → I
    Infectious,    # I pool – infection timer countdown, I → R
    Recovered,     # R pool – immunity duration tracking
    Transmission,  # Force-of-infection engine – S → E
)

params = PropertySet({
    "population": 100_000,
    "nodes":      1,
    "nticks":     365,        # daily time steps
    "beta":       0.5,        # transmission rate
    "sigma":      1 / 8,      # 1 / incubation period (days⁻¹)
    "gamma":      1 / 5,      # 1 / infectious period (days⁻¹)
    "seed":       42,
})

model = Model(params)

# Assign components in execution order
model.components = [
    Susceptible,   # 1
    Exposed,       # 2
    Infectious,    # 3
    Recovered,     # 4
    Transmission,  # 5  ← always last
]

model.run()
```

`model.components` setter instantiates each class (passing `model` to its constructor) and appends a callable to an internal `phases` list. `model.run()` loops over ticks and calls each phase in order.

---

### Why this order is mandatory

Each component's `.step(model, tick)` is called **once per tick, sequentially**. All state reads within a tick reflect the **start-of-tick snapshot**.

| Step | Component | Transition | Reads state from |
|------|-----------|------------|-----------------|
| 1 | `Susceptible` | births / waning → S | start of tick |
| 2 | `Exposed` | E → I (incubation expiry) | start of tick |
| 3 | `Infectious` | I → R (recovery) | start of tick |
| 4 | `Recovered` | R bookkeeping | start of tick |
| 5 | `Transmission` | **S → E** (new exposures) | start of tick |

**`Transmission` must run last** for two related reasons:

1. **Prevents double-processing.** If `Transmission` ran before `Exposed`, agents freshly moved from S → E this tick would be immediately seen by `Exposed.step()` and advanced toward I in the *same tick*. That collapses the incubation period by one day for every newly exposed agent.

2. **Preserves correct discrete-time semantics.** All progression steps (E → I, I → R) should operate on the cohort that was *already* in that state at the start of the tick. New exposures created by `Transmission` belong to the *next* tick's cohort.

The official LASER SIR example confirms this pattern — `IntrahostProgression` (I → R) is registered before `Transmission` (S → I) for exactly the same reason.

---

Sources:
- [SIR Model Example — LASER 0.3.0 documentation](https://docs.idmod.org/projects/laser/en/v0.3.0/example.html)
- [Software Architecture & Design — LASER 0.6.0](https://docs.idmod.org/projects/laser/en/latest/architecture.html)
- [laser_cholera metapop model source](https://docs.idmod.org/projects/laser-cholera/en/v0.8.0/_modules/laser_cholera/metapop/model.html)
- [GitHub — InstituteforDiseaseModeling/laser](https://github.com/InstituteforDiseaseModeling/laser)
