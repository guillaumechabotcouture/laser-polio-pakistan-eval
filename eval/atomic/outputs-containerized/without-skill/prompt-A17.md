The LASER framework does not export named integer constants like `SUSCEPTIBLE = 0`. Instead, states are defined as a positional string list and the integer indices are derived from it. Here is the code showing exactly where they come from and what each value is:

```python
from laser_measles.abm.params import STATES

# STATES = ["S", "E", "I", "R"]
# Integer values are list indices — derive them explicitly:

SUSCEPTIBLE = STATES.index("S")
EXPOSED     = STATES.index("E")
INFECTIOUS  = STATES.index("I")
RECOVERED   = STATES.index("R")

print(f"SUSCEPTIBLE = {SUSCEPTIBLE}")  # 0
print(f"EXPOSED     = {EXPOSED}")      # 1
print(f"INFECTIOUS  = {INFECTIOUS}")   # 2
print(f"RECOVERED   = {RECOVERED}")    # 3
```

**Output:**
```
SUSCEPTIBLE = 0
EXPOSED     = 1
INFECTIOUS  = 2
RECOVERED   = 3
```

**Key points:**

| State | Value | Source |
|---|---|---|
| `SUSCEPTIBLE` | `0` | `STATES.index("S")` |
| `EXPOSED` | `1` | `STATES.index("E")` |
| `INFECTIOUS` | `2` | `STATES.index("I")` |
| `RECOVERED` | `3` | `STATES.index("R")` |

- Import from: `laser_measles.abm.params`
- The framework does **not** export pre-named integer constants — you derive them from `STATES`, which is `["S", "E", "I", "R"]`
- The framework's own internal Numba/NumPy kernels use the raw literals `0`, `1`, `2`, `3` directly (e.g., `states == 0` for susceptible checks in `process_transmission.py`)
