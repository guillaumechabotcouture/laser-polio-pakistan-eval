Here's how to add mortality using `MortalityByCDR`:

## Units

`MortalityByCDR` expects the crude death rate in **per-1000 population per year**. It divides by 1000 internally to get a per-capita rate, then converts from annual to daily. So you pass `13` directly — do **not** convert to a daily per-capita fraction.

## Complete code

```python
import numpy as np
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import MortalityByCDR

# CDR = 13 per 1000 population per year — pass as-is
CDR = 13

# Build a ValuesMap of shape (nticks, nnodes)
# ValuesMap.from_scalar broadcasts a single value across all ticks and nodes
nticks = 20 * 365
nnodes = 13  # number of patches in your scenario
deathrate_map = ValuesMap.from_scalar(CDR, nticks, nnodes)

# Validate units (should be in the 5-60 range, never a tiny decimal)
assert np.all(deathrate_map.values >= 1) and np.all(deathrate_map.values <= 60), \
    f"CDR must be per-1000/year, got {deathrate_map.values.min():.4f}-{deathrate_map.values.max():.4f}"

# Add to model components list (after births, before or after transmission)
model.components = [
    # ... other components ...
    BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=deathrate_map.values),
    # ... transmission, vaccination, etc. ...
]
```

**Two things to watch:**

1. **Parameter name is `mortalityrates=`**, not `deathrates=` or `cdr=`.
2. **Pass the raw number `13`**, not `0.013` (per-capita) or `0.0000356` (daily per-capita). LASER handles the conversion internally.
