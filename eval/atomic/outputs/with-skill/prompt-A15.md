Now I have what I need. `AliasedDistribution` takes an array of weights (counts/probabilities) where each index represents that age (in days, typically). Let me write the solution.

Here's how to create and sample from an `AliasedDistribution` representing your age pyramid:

```python
import numpy as np
from laser.core.demographics import AliasedDistribution

# Define age bins and their probabilities
# Bins: 0-4 (30%), 5-14 (25%), 15-29 (20%), 30-49 (15%), 50+ (10%)
# AliasedDistribution assigns one weight per index, so we spread
# each bin's probability uniformly across its constituent ages (in years).

max_age = 80  # upper bound for 50+ bin
weights = np.zeros(max_age)

weights[0:5]   = 0.30 / 5    # 30% spread over 5 ages
weights[5:15]  = 0.25 / 10   # 25% spread over 10 ages
weights[15:30] = 0.20 / 15   # 20% spread over 15 ages
weights[30:50] = 0.15 / 20   # 15% spread over 20 ages
weights[50:80] = 0.10 / 30   # 10% spread over 30 ages

# Create the aliased distribution (O(1) sampling via alias method)
pyramid = AliasedDistribution(weights)

# Sample 1000 ages (returns integer indices = ages in years)
ages = pyramid.sample(count=1000, dtype=np.int32)

# Verify the distribution matches expectations
print(f"Sampled {len(ages)} ages")
print(f"  0-4:   {np.mean((ages >= 0) & (ages < 5)):.1%}  (target 30%)")
print(f"  5-14:  {np.mean((ages >= 5) & (ages < 15)):.1%}  (target 25%)")
print(f"  15-29: {np.mean((ages >= 15) & (ages < 30)):.1%}  (target 20%)")
print(f"  30-49: {np.mean((ages >= 30) & (ages < 50)):.1%}  (target 15%)")
print(f"  50+:   {np.mean(ages >= 50):.1%}  (target 10%)")
```

**How it works:**

- `AliasedDistribution` takes an array of weights — one per index. When you call `.sample()`, it returns random indices drawn proportionally to those weights using the [alias method](https://en.wikipedia.org/wiki/Alias_method) (O(n) setup, O(1) per sample).
- Each index corresponds to an age. We spread each bin's total probability uniformly across its constituent ages so that, e.g., a 2-year-old and a 4-year-old are equally likely within the 0-4 bin.
- In LASER models, ages are typically tracked in **days**. If you need day-resolution ages, multiply `max_age` by 365 and adjust the bin boundaries accordingly (e.g., `weights[0:5*365] = 0.30 / (5*365)`).
