Here are the SEIR state values from the LASER framework:

**Import path:** `from laser.generic.SEIR import State`
(or directly: `from laser.generic.shared import State`)

**Integer values (stored as `np.int8`):**

| State | Value |
|-------|-------|
| `State.DECEASED` | -1 |
| `State.SUSCEPTIBLE` | 0 |
| `State.EXPOSED` | 1 |
| `State.INFECTIOUS` | 2 |
| `State.RECOVERED` | 3 |

Defined in `/opt/anaconda3/lib/python3.12/site-packages/laser/generic/shared.py:30-40`.

Usage example:

```python
from laser.generic.SEIR import State

# Access values
State.SUSCEPTIBLE.value  # np.int8(0)
State.EXPOSED.value      # np.int8(1)
State.INFECTIOUS.value   # np.int8(2)
State.RECOVERED.value    # np.int8(3)
State.DECEASED.value     # np.int8(-1)
```

Note: In project code you'll also see the older import style `from laser.generic.SEIR import SEIR` followed by `SEIR.State.SUSCEPTIBLE.value` — both work, but the direct `State` import shown above is cleaner since `State` is re-exported from the `SEIR` module.
