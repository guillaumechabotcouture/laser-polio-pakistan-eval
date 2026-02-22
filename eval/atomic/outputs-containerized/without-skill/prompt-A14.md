The `/eval-work` filesystem is read-only, so I've written the file to `/tmp/pulse_importation.py`. The component passes a smoke test. Here is the full implementation:

```python
import numpy as np


class PulseImportation:
    """
    LASER component that seeds 5 infectious agents into patch 0
    at tick 0 and again at every subsequent multiple of 90 ticks.
    """

    _SUSCEPTIBLE: int = 0   # laser.generic.shared.State.SUSCEPTIBLE
    _INFECTIOUS:  int = 2   # laser.generic.shared.State.INFECTIOUS

    _PULSE_INTERVAL: int = 90
    _N_INFECTIONS:   int = 5
    _TARGET_PATCH:   int = 0

    def __init__(self, model) -> None:
        self.model = model

    def step(self, tick: int) -> None:
        if tick % self._PULSE_INTERVAL != 0:
            return

        people = self.model.people

        susceptible_in_patch = (
            (people.nodeid == self._TARGET_PATCH)
            & (people.state == self._SUSCEPTIBLE)
        )
        candidates = np.nonzero(susceptible_in_patch)[0]

        n = min(self._N_INFECTIONS, len(candidates))
        if n == 0:
            return

        rng = getattr(self.model, "prng", np.random.default_rng())
        chosen = rng.choice(candidates, size=n, replace=False)

        people.state[chosen] = self._INFECTIOUS

        if hasattr(people, "_itimer"):
            params = getattr(self.model, "params", None)
            infectious_period = int(getattr(params, "infectious_period", 10))
            people.itimer[chosen] = infectious_period
```

**How it follows the LASER component protocol:**

| Requirement | Implementation |
|---|---|
| `__init__(self, model)` | Stores the model reference; no side effects on `model.people` at construction time |
| `step(self, tick)` | Called once per tick by `model.run()`; applies the pulse logic |
| Agent state access | `people.nodeid` and `people.state` — standard LaserFrame scalar properties returned as NumPy views |
| State constants | `0` = SUSCEPTIBLE, `2` = INFECTIOUS — matching `laser.generic.shared.State` values |

**Pulse logic:**
- Fires at tick 0, 90, 180, … via `tick % 90 == 0`
- Finds all susceptible agents in patch 0 with a boolean NumPy mask, selects up to 5 at random
- Uses `model.prng` (the model's seeded RNG) when available, for reproducibility
- Initialises `itimer` when that property exists on the frame (SIR/SEIRS models), gracefully skipping it otherwise (SI models)
