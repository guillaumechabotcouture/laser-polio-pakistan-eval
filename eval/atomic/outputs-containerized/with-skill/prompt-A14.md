`PulseImportation` already exists in `scripts/custom_components.py` at line 410, implementing exactly what you described. Here's what it does:

```python
class PulseImportation:
    def __init__(self, model, infdurdist, period=90, count=5):
        self.model = model
        self.infdurdist = infdurdist
        self.period = period
        self.count = count

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0:
            return
        ...
```

**Behavior summary** (`scripts/custom_components.py:410-458`):

- `__init__` stores the model reference and a distribution for sampling infectious durations. `period=90` and `count=5` are the defaults matching your spec.
- `step` skips tick 0 and any tick not divisible by `period`, then:
  1. Finds all susceptible agents (`state == SUSCEPTIBLE`) currently in patch 0 via boolean mask + `np.nonzero`.
  2. Selects `min(count, available)` of them with `np.random.choice(..., replace=False)`.
  3. Sets `people.state[chosen] = INFECTIOUS` and assigns sampled `itimer` durations (floored to ≥1).
  4. Decrements `nodes.S[tick+1, 0]` and increments `nodes.I[tick+1, 0]` to keep compartment counts consistent.

To use it, pass an infectious duration distribution at construction:

```python
from custom_components import PulseImportation

importation = PulseImportation(model, infdurdist=my_infdur_dist)
model.components.append(importation)
```

No changes were needed — the implementation already satisfies the LASER component protocol and your requirements.
