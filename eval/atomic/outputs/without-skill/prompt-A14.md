The component is written to `pulse_importation.py`. Here's how it works:

**`__init__(self, model)`** — Stores a reference to the model and configures the pulse parameters (every 90 ticks, 5 infections, targeting patch 0).

**`step(self, tick)`** — On each tick divisible by 90, it:
1. Finds susceptible agents (`disease_state == 0`) located in patch 0 (`node_id == 0`)
2. Randomly selects up to 5 of them
3. Sets their `disease_state` to `1` (infected)

If fewer than 5 susceptible agents remain in patch 0, it infects all that are available. If none remain, it's a no-op.

**Usage** with a LASER model:

```python
model.add_component(PulseImportation(model))
```
