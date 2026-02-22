The filesystem is read-only, so I'll present the code directly.

---

## How to correctly vaccinate agents in LASER with `TransmissionSE`

### Why the state field — not susceptibility — is what matters

`TransmissionSE` uses a Numba-compiled inner loop that looks like this (from `laser.generic.components`):

```python
for i in nb.prange(len(states)):
    if states[i] == State.SUSCEPTIBLE.value:   # checks state, not susceptibility
        draw = np.random.rand()
        if draw < ft[nid]:
            states[i] = State.EXPOSED.value    # infect
```

Setting `model.people.susceptibility[i] = 0` has **no effect** here — that field is never read. The only way to make an agent immune is to give them a state value that is not `State.SUSCEPTIBLE` (0). The conventional immune state is `State.RECOVERED` (3), which is also what LASER's own `RoutineImmunizationEx` component uses.

---

### State enum (from `laser.generic.shared`)

```python
class State(Enum):
    DECEASED    = -1   # int8
    SUSCEPTIBLE =  0   # int8  ← TransmissionSE infects agents HERE
    EXPOSED     =  1   # int8
    INFECTIOUS  =  2   # int8
    RECOVERED   =  3   # int8  ← vaccinated agents must land HERE
```

---

### Reusable `VaccinationCampaign` component

```python
import numpy as np
from laser.generic.shared import State


class VaccinationCampaign:
    """
    One-time mass vaccination campaign.

    Sets state = RECOVERED on a random `coverage` fraction of SUSCEPTIBLE
    agents so that TransmissionSE skips them entirely.

    Wrong approach: model.people.susceptibility[chosen] = 0
      → TransmissionSE never reads susceptibility; agents still get infected.

    Correct approach: model.people.state[chosen] = State.RECOVERED.value
      → TransmissionSE's `if states[i] == State.SUSCEPTIBLE.value` gate
        excludes them; they cannot be infected.
    """

    def __init__(self, model, coverage: float = 0.80, tick: int = 0, seed: int = None):
        if not (0.0 <= coverage <= 1.0):
            raise ValueError(f"coverage must be in [0, 1], got {coverage}")
        self.model    = model
        self.coverage = coverage
        self.tick     = tick
        self.rng      = np.random.default_rng(seed)

    def step(self, tick: int) -> None:
        if tick != self.tick:
            return

        states = self.model.people.state                          # np.int8 array

        susceptible_idx = np.nonzero(states == State.SUSCEPTIBLE.value)[0]
        n_susceptible   = len(susceptible_idx)

        if n_susceptible == 0:
            return

        n_vaccinate = int(np.round(n_susceptible * self.coverage))
        chosen      = self.rng.choice(susceptible_idx, size=n_vaccinate, replace=False)

        # ---------------------------------------------------------------
        # THE KEY LINE: set state to RECOVERED, not susceptibility to 0.
        # TransmissionSE checks states[i] == State.SUSCEPTIBLE.value (0).
        # Any other value — including RECOVERED (3) — passes the gate check
        # and the agent is never exposed.
        # ---------------------------------------------------------------
        states[chosen] = State.RECOVERED.value                    # integer 3

        print(
            f"[tick {tick}] Vaccinated {n_vaccinate:,} / {n_susceptible:,} "
            f"susceptible agents ({100*n_vaccinate/n_susceptible:.1f}%)"
        )
```

---

### Wiring it into a model

```python
from laser.generic.model  import LaserModel
from laser.generic        import components as comp

params = {
    "nticks"     : 365,
    "nodes"      : 1,
    "population" : [10_000],
    "prevalence" : [0.01],
    "beta"       : [0.3],
    "exp_mean"   : 5.0,  "exp_std" : 1.0,
    "inf_mean"   : 10.0, "inf_std" : 2.0,
}

model = LaserModel(params)

# Standard SEIR components
model.components.append(comp.Susceptible(model))
model.components.append(comp.Exposed(model))
model.components.append(comp.InfectiousIR(model))
model.components.append(comp.TransmissionSE(model))

# Vaccination campaign: 80 % coverage on day 0, before any transmission
model.components.append(VaccinationCampaign(model, coverage=0.80, tick=0, seed=42))

model.run()
```

---

### Demonstrating the difference

```python
import numpy as np
from laser.generic.shared import State

n        = 1_000
states   = np.zeros(n, dtype=np.int8)   # all SUSCEPTIBLE
foi      = 0.5                          # force of infection

# ---- Wrong: susceptibility field (TransmissionSE ignores it) ----
susceptibility = np.ones(n, dtype=np.float32)
susceptibility[:800] = 0.0              # "vaccinate" 800 agents

infected = sum(
    1 for i in range(n)
    if states[i] == State.SUSCEPTIBLE.value   # all 1000 still pass!
    and np.random.rand() < foi
)
print(f"Wrong  (susceptibility=0): {infected} infected  ← all 1000 were eligible")

# ---- Correct: state = RECOVERED ----
states_v = np.zeros(n, dtype=np.int8)
vax_idx  = np.random.choice(n, size=800, replace=False)
states_v[vax_idx] = State.RECOVERED.value       # only 200 remain SUSCEPTIBLE

infected_v = sum(
    1 for i in range(n)
    if states_v[i] == State.SUSCEPTIBLE.value   # only 200 pass
    and np.random.rand() < foi
)
print(f"Correct (state=RECOVERED): {infected_v} infected  ← only 200 were eligible")
```

---

### Summary

| Approach | Field modified | TransmissionSE sees it? | Result |
|---|---|---|---|
| `pop.susceptibility[i] = 0` | `susceptibility` | No | Agent still infected |
| `pop.state[i] = State.RECOVERED.value` | `state` | **Yes** | Agent fully immune |

The `VaccinationCampaign.step()` component can be appended to `model.components` and will fire exactly once at the specified tick, converting 80 % of susceptibles to `RECOVERED` before `TransmissionSE` has a chance to expose them.
