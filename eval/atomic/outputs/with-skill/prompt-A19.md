Here's the core issue and the correct solution:

## The LASER Vaccination Trap

LASER's `TransmissionSE` (and all Transmission kernels) determines who can be infected by checking **`state == SUSCEPTIBLE`** (int8 value 0). It does **not** check the `susceptibility` property. This means:

```python
# WRONG — has zero effect on transmission
people.susceptibility[indices] = 0.0   # TransmissionSE ignores this field
```

```python
# CORRECT — prevents infection by removing agents from the SUSCEPTIBLE state
people.state[indices] = SEIR.State.RECOVERED.value  # int8 == 3
```

## Minimal Correct Example: Vaccinate 80% of Susceptibles

```python
import numpy as np
from laser.generic import SEIR

def vaccinate_80pct(model, tick):
    """Vaccinate 80% of susceptible agents by setting state = RECOVERED.
    
    LASER's TransmissionSE only checks `state == SUSCEPTIBLE` (int8 == 0)
    to select agents for infection. It does NOT check `susceptibility`.
    Setting susceptibility = 0 has zero effect. You must change the state.
    """
    people = model.people
    nodes = model.nodes
    count = people.count

    # 1. Find all susceptible agents
    susceptible_idx = np.nonzero(
        people.state[:count] == SEIR.State.SUSCEPTIBLE.value  # int8 == 0
    )[0]

    if len(susceptible_idx) == 0:
        return

    # 2. Randomly select 80% of them
    n_to_vaccinate = int(0.80 * len(susceptible_idx))
    chosen = np.random.choice(susceptible_idx, size=n_to_vaccinate, replace=False)

    # 3. Set state = RECOVERED (value 3) — this is what makes them immune
    people.state[chosen] = SEIR.State.RECOVERED.value

    # 4. Update node-level compartment counts (S decreases, R increases)
    vax_by_node = np.bincount(
        people.nodeid[chosen], minlength=nodes.count
    ).astype(nodes.S.dtype)

    nodes.S[tick + 1] -= vax_by_node
    np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])  # guard against negatives
    nodes.R[tick + 1] += vax_by_node
```

## Why Both Steps (3 and 4) Matter

- **Step 3** (`people.state = RECOVERED`) is what actually prevents `TransmissionSE` from infecting the agent. This operates at the **agent level**.
- **Step 4** (`nodes.S -= ...`, `nodes.R += ...`) keeps the **node-level compartment counts** consistent. `TransmissionSE` uses `nodes.I` and `nodes.S` to compute force of infection, and the `Susceptible`/`Recovered` components propagate these counts forward each tick. If you skip this step, the node-level bookkeeping drifts out of sync with the actual agent states.

## As a Reusable Component

To plug this into a LASER model's component list:

```python
class SimpleVaccination:
    """One-time 80% vaccination of all susceptibles at a given tick."""

    def __init__(self, model, coverage=0.80, vaccinate_at_tick=1):
        self.model = model
        self.coverage = coverage
        self.vaccinate_at_tick = vaccinate_at_tick

    def step(self, tick):
        if tick != self.vaccinate_at_tick:
            return

        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        susceptible_idx = np.nonzero(
            people.state[:count] == SEIR.State.SUSCEPTIBLE.value
        )[0]

        if len(susceptible_idx) == 0:
            return

        n_vax = int(self.coverage * len(susceptible_idx))
        chosen = np.random.choice(susceptible_idx, size=n_vax, replace=False)

        # Change STATE, not susceptibility — this is the critical line
        people.state[chosen] = SEIR.State.RECOVERED.value

        vax_by_node = np.bincount(
            people.nodeid[chosen], minlength=nodes.count
        ).astype(nodes.S.dtype)
        nodes.S[tick + 1] -= vax_by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.R[tick + 1] += vax_by_node

        print(f"Tick {tick}: vaccinated {len(chosen):,} agents "
              f"({100*len(chosen)/len(susceptible_idx):.0f}% of susceptibles)")
```

Then add it to the model's component list:

```python
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SimpleVaccination(model, coverage=0.80, vaccinate_at_tick=1),
    SEIR.Transmission(model, expdurdist),
    # ... other components
]
```

The existing codebase at `scripts/custom_components.py` has production-grade examples (e.g., `PerPatchVaccination` at line 290) that follow this same pattern — setting `state = RECOVERED` and updating node counts — plus features like per-patch coverage arrays and correlated missedness for hard-to-reach populations.
