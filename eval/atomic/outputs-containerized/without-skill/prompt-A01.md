Now I have all the information needed. Here is the complete Python code:

```python
import numpy as np
from laser.core import PropertySet
from laser.core.utils import grid
from laser.generic import Model

# ── Scenario: 4-patch (2×2) GeoDataFrame ─────────────────────────────────────
# grid() returns a GeoDataFrame with columns: nodeid, population, geometry,
# S, E, I, R — with S already set to the full population (all susceptible).
scenario = grid(
    M=2,
    N=2,
    population_fn=lambda row, col: 50_000,
)

# Explicit: all individuals susceptible, none exposed/infectious/recovered.
scenario["S"] = scenario["population"]
scenario["E"] = 0
scenario["I"] = 0
scenario["R"] = 0

# ── Parameters ────────────────────────────────────────────────────────────────
params = PropertySet({
    "nticks":    365,
    "prng_seed": 20250101,
})

# ── Model ─────────────────────────────────────────────────────────────────────
model = Model(scenario, params)

# ── Inspect ───────────────────────────────────────────────────────────────────
print(f"Patches : {len(model.scenario)}")
print(f"Agents  : {model.people.count:,}")
print(f"Nodes   : {model.nodes.count}")
print(scenario[["nodeid", "population", "S", "E", "I", "R"]])
```

**Key points:**

- **`grid(M=2, N=2, ...)`** — creates a 2×2 grid of 4 rectangular patches in geographic (EPSG:4326) coordinates, each with a unique `nodeid` (0–3). Any `M×N` combination giving 4 patches works (`M=1, N=4`, etc.). The function sets `S = population` and `E = I = R = 0` by default.
- **`population_fn=lambda row, col: 50_000`** — assigns 50,000 to every patch regardless of position.
- **`PropertySet({"nticks": 365, "prng_seed": 20250101})`** — `Model.__init__` searches for `prng_seed`, then `prngseed`, then `seed` in params (in that order) to seed the PRNG via `laser.core.random.seed`.
- **`Model(scenario, params)`** — constructs the agent pool (`model.people` as a `LaserFrame`), node data (`model.nodes`), and a gravity-based inter-patch migration network. No `beta` or other disease parameters are required at construction time; those are consumed by components when added.

Sources:
- [laser-base/laser-generic (GitHub)](https://github.com/laser-base/laser-generic)
- [laser-generic on PyPI](https://pypi.org/project/laser.generic/)
