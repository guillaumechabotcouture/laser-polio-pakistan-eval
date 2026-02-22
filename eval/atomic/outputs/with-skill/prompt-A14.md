The `PulseImportation` class has been added to `scripts/custom_components.py` (line 408). Here's what it does:

- **`__init__(self, model, infdurdist, period=90, count=5)`** — stores the model reference, infectious duration distribution, pulse interval (90 ticks), and infection count (5).
- **`step(self, tick)`** — every 90 ticks, finds susceptible agents in patch 0, selects up to 5, and transitions them S→I with sampled infectious durations. Updates node-level `S` and `I` counts accordingly.

It follows the same patterns as the existing `PatchImportation` component: targets only susceptible agents (avoiding double-infection), uses `dists.sample_floats` for infectious duration, and correctly adjusts both agent state and node-level compartment counts.

Usage in a model's component list:

```python
PulseImportation(model, infdurdist)  # defaults: 5 infections into patch 0 every 90 ticks
```
