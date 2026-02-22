Here's how to create and use a LASER `PropertySet`:

```python
from laser.core.propertyset import PropertySet

# Create a PropertySet with standard and custom parameters
params = PropertySet({
    "nticks": 3650,
    "beta": 0.25,
    "gamma": 1 / 14,
    "sigma": 1 / 5,
    "prng_seed": 42,
    "seasonal_amplitude": 1.3,
})

# Access individual parameters as attributes
print(f"nticks:              {params.nticks}")
print(f"beta:                {params.beta}")
print(f"gamma:               {params.gamma:.6f}")
print(f"sigma:               {params.sigma:.4f}")
print(f"prng_seed:           {params.prng_seed}")
print(f"seasonal_amplitude:  {params.seasonal_amplitude}")
```

`PropertySet` wraps a dictionary and exposes every key as a dot-accessible attribute — both built-in LASER parameters (like `nticks`, `prng_seed`) and any custom ones you define (like `seasonal_amplitude`) work the same way.
