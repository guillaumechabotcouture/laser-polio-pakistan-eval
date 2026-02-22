The states are defined as a numpy `int8` enum in `laser.generic.shared` and re-exported from `laser.generic.SEIR`.

```python
from laser.generic.SEIR import State  # or: from laser.generic.shared import State

# Integer values (stored as numpy.int8)
print(State.SUSCEPTIBLE.value)  # 0
print(State.EXPOSED.value)      # 1
print(State.INFECTIOUS.value)   # 2
print(State.RECOVERED.value)    # 3

# Bonus: DECEASED is also in the enum
print(State.DECEASED.value)     # -1
```

Summary:

| State | Value |
|-------|-------|
| `SUSCEPTIBLE` | `0` |
| `EXPOSED` | `1` |
| `INFECTIOUS` | `2` |
| `RECOVERED` | `3` |
| `DECEASED` | `-1` |

The canonical source is `laser.generic.shared`, but `laser.generic.SEIR` re-exports `State` and is what the skill docs and most model code use. Both import paths work. Values are `numpy.int8`, so use `.value` to get the raw integer when needed (e.g., for array comparisons: `population.states == State.SUSCEPTIBLE.value`).
