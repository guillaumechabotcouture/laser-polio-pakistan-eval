The filesystem is read-only, so here is the complete, runnable code with explanation:

```python
from laser_core.propertyset import PropertySet

# Create a PropertySet with SEIR model parameters plus a custom seasonal parameter
params = PropertySet({
    "nticks":             3650,
    "beta":               0.25,
    "gamma":              1 / 14,
    "sigma":              1 / 5,
    "prng_seed":          42,
    "seasonal_amplitude": 1.3,
})
```

`PropertySet` accepts a plain dict in its constructor. Values are stored internally — `gamma` and `sigma` are evaluated immediately to their float equivalents (`0.07142...` and `0.2`).

---

**Accessing individual parameters**

```python
# Dot notation (preferred in LASER model code)
params.nticks             # → 3650
params.beta               # → 0.25
params.gamma              # → 0.07142857142857142
params.sigma              # → 0.2
params.prng_seed          # → 42
params.seasonal_amplitude # → 1.3

# Dictionary-style (interchangeable)
params["beta"]            # → 0.25
params["seasonal_amplitude"]  # → 1.3
```

---

**Other useful operations**

```python
# Membership test
"seasonal_amplitude" in params   # → True
"mu" in params                   # → False

# Number of stored parameters
len(params)                      # → 6

# Print all params as formatted JSON
print(params)

# Add a new parameter after construction
params.mu = 1 / 70 / 365        # birth/death rate

# Merge with another PropertySet (raises ValueError on duplicate keys)
overrides = PropertySet({"beta": 0.30})
merged = params + overrides      # combines into a new PropertySet
```

---

**Key points:**

| Feature | Detail |
|---|---|
| Dot access | `params.gamma` — clean, readable in model code |
| Bracket access | `params["gamma"]` — useful when key is a variable |
| Custom params | Any key works (`seasonal_amplitude`, etc.) — no schema enforced |
| Computed values | `1/14` and `1/5` are evaluated at assignment time to floats |
| `prng_seed` | Passed to `np.random.default_rng(params.prng_seed)` in LASER components |
