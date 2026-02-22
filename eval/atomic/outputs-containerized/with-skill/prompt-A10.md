Here is the correct component assembly with explanation:

```python
import laser.core.distributions as dists
from laser.generic import SEIR, Model

# Distribution objects required by the transition components
expdurdist = dists.gamma(shape=40, scale=0.25)   # Exposed duration ~10 days
infdurdist = dists.normal(loc=8, scale=2)         # Infectious duration ~8 days

model = Model(scenario, parameters)

model.components = [
    SEIR.Susceptible(model),                             # 1
    SEIR.Exposed(model, expdurdist, infdurdist),         # 2
    SEIR.Infectious(model, infdurdist),                  # 3
    SEIR.Recovered(model),                               # 4
    SEIR.Transmission(model, expdurdist),                # 5
]
```

**Why this order matters**

Each tick, components execute sequentially. The ordering is designed around a single invariant: `S[t+1] + E[t+1] + I[t+1] + R[t+1] = N` must hold after every tick.

| Position | Component | What it does at tick `t` |
|---|---|---|
| 1 | `Susceptible` | Copies `S[t]` → `S[t+1]` as a baseline for this tick |
| 2 | `Exposed` | Decrements `etimer` per agent; agents at 0 transition E→I, setting `itimer` via `infdurdist` |
| 3 | `Infectious` | Decrements `itimer` per agent; agents at 0 transition I→R |
| 4 | `Recovered` | Copies `R[t]` → `R[t+1]` as a baseline |
| 5 | `Transmission` | Reads `I[t]` to compute force of infection; new S→E exposures subtract from `S[t+1]`, add to `E[t+1]` |

**The wrapping pattern**

`Susceptible` (step 1) and `Recovered` (step 4) establish the `t+1` baselines before any compartment counts are modified. The transition components (`Exposed`, `Infectious`) then operate on the `t+1` arrays in between. `Transmission` runs last so it uses the fully-settled `I[t]` count — before any new exposures have been added — to compute an unbiased force of infection, and then modifies `S[t+1]` and `E[t+1]`.

Placing `Transmission` before `Exposed` or `Infectious` would mean new exposures are counted before timer-based transitions resolve, corrupting the `I[t]` count used in the FOI calculation and breaking the population invariant.

**Component aliases**

For reference, the `SEIR.*` names map to the underlying classes:

| Alias | Class |
|---|---|
| `SEIR.Susceptible` | `Susceptible` |
| `SEIR.Exposed` | `Exposed` |
| `SEIR.Infectious` | `InfectiousIR` |
| `SEIR.Recovered` | `Recovered` |
| `SEIR.Transmission` | `TransmissionSE` |

Both import styles work:
```python
from laser.generic import SEIR
SEIR.Susceptible(model)  # via alias

from laser.generic.components import Susceptible, Exposed, InfectiousIR, Recovered, TransmissionSE
Susceptible(model)       # direct
```
