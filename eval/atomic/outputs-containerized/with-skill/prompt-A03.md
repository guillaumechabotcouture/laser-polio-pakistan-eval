The API reference confirms the `PropertySet` interface. Here's a self-contained example:

```python
from laser.core.propertyset import PropertySet

# Create PropertySet with standard and custom parameters
params = PropertySet({
    "nticks":             3650,
    "beta":               0.25,
    "gamma":              1 / 14,
    "sigma":              1 / 5,
    "prng_seed":          42,
    "seasonal_amplitude": 1.3,   # custom parameter — PropertySet accepts any key
})

# Attribute-style access
print(params.nticks)              # 3650
print(params.beta)                # 0.25
print(params.gamma)               # 0.07142857142857142
print(params.sigma)               # 0.2
print(params.prng_seed)           # 42
print(params.seasonal_amplitude)  # 1.3

# Dictionary-style access (equivalent)
print(params["beta"])             # 0.25

# Inspect all parameters
print(params.to_dict())
```

A few notes:

- `PropertySet` accepts any key/value pair — `seasonal_amplitude` is treated identically to built-in params. LASER only reads the keys it recognizes (like `gravity_k`, `nticks`); custom keys are silently ignored by the framework and available only to your code.
- `gamma` and `sigma` are stored as Python floats (`1/14 ≈ 0.0714`, `1/5 = 0.2`) since `PropertySet` doesn't impose types.
- LASER's `TransmissionSE` uses `beta` directly as a per-day rate; `gamma`/`sigma` are not consumed by built-in components — you'd pass them to distribution constructors (e.g., `dists.constant(1/params.gamma)` for mean infectious duration).
- Both access styles (`params.key` and `params["key"]`) are identical.
