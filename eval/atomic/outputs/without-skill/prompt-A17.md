Here's the result. The `State` enum is defined in `laser.generic.shared` and re-exported from `laser.generic`:

```python
from laser.generic import State
```

The four SEIR states and their integer values (stored as `numpy.int8`):

| State | Value |
|---|---|
| `State.SUSCEPTIBLE` | **0** |
| `State.EXPOSED` | **1** |
| `State.INFECTIOUS` | **2** |
| `State.RECOVERED` | **3** |

There is also a `State.DECEASED = -1` defined in the same enum but it's outside the standard SEIR compartments.

The underlying enum (from `laser/generic/shared.py`) stores values as `np.int8` via a custom `__new__` method, which matches the `int8` dtype used for the `state` property on the agent population (e.g., `model.people.state`). In model code you'll typically see comparisons like `model.people.state == State.INFECTIOUS.value`.
