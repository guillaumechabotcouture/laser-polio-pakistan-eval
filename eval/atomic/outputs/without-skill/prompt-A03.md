The script runs successfully. Here's a summary of what it demonstrates:

**Creating a `PropertySet`** — Import from `laser.core` and pass dictionaries to the constructor. Each dictionary becomes an accessible parameter:

```python
from laser.core import PropertySet

params = PropertySet(
    {"nticks": 3650},
    {"beta": 0.25},
    {"gamma": 1 / 14},
    {"sigma": 1 / 5},
    {"prng_seed": 42},
    {"seasonal_amplitude": 1.3},
)
```

**Accessing individual parameters** — Two equivalent ways:
- **Dot notation:** `params.beta` → `0.25`
- **Bracket notation:** `params['beta']` → `0.25`

PropertySet also supports the `<<` operator (override existing keys only) and `|` operator (add or override keys), useful for merging parameter sets — e.g., overriding defaults with scenario-specific values.

The derived R0 = beta/gamma = 0.25 / (1/14) = **3.50**, consistent with a moderately transmissible disease like measles in a partially vaccinated population.
