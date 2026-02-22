The API reference shows that `AliasedDistribution` takes an array of weights (unnormalized is fine) and returns sampled indices. For an age pyramid, each index represents one year of age. Here's the complete example:

```python
import numpy as np
from laser.core.demographics import AliasedDistribution

# Age pyramid: weights per year-of-age, proportional to each bin's share.
# Sampled index = age in years, so we need weights[0..max_age-1].
#
# 50+ is modeled as 50–79 (30 years), a common demographic ceiling.
max_age = 80  # ages 0–79

weights = np.zeros(max_age)
weights[0:5]   = 0.30 / 5    # 0–4:   30% spread over 5 years
weights[5:15]  = 0.25 / 10   # 5–14:  25% spread over 10 years
weights[15:30] = 0.20 / 15   # 15–29: 20% spread over 15 years
weights[30:50] = 0.15 / 20   # 30–49: 15% spread over 20 years
weights[50:80] = 0.10 / 30   # 50–79: 10% spread over 30 years

# AliasedDistribution accepts unnormalized weights, but this is already ~1.0
assert abs(weights.sum() - 1.0) < 1e-10, \
    f"Weights sum to {weights.sum():.8f}, expected 1.0"

pyramid = AliasedDistribution(weights)

# Sample 1000 ages — returned values are integer year-of-age indices
ages = pyramid.sample(count=1000, dtype=np.int32)

# --- Verification ---
print(f"Sampled {len(ages)} ages  |  min={ages.min()}  max={ages.max()}  mean={ages.mean():.1f}\n")

groups = [
    ("0–4",   0.30, (0,  5)),
    ("5–14",  0.25, (5,  15)),
    ("15–29", 0.20, (15, 30)),
    ("30–49", 0.15, (30, 50)),
    ("50+",   0.10, (50, 80)),
]

print(f"{'Age group':<10} {'Target':>8} {'Observed':>10}")
print("-" * 30)
for label, target, (lo, hi) in groups:
    observed = ((ages >= lo) & (ages < hi)).sum() / len(ages)
    flag = "  ✓" if abs(observed - target) < 0.05 else "  !"
    print(f"{label:<10} {target:>8.0%} {observed:>10.1%}{flag}")
```

**Key points:**

- `AliasedDistribution(weights)` — weights are per-year-of-age; unnormalized is fine (it normalizes internally). Spreading each bin's proportion uniformly across its years gives correct sampling.
- `pyramid.sample(count=N, dtype=np.int32)` — returns an integer array of sampled year-of-age values.
- The "50+" bin requires a concrete upper bound. 80 years (`50:80`) is standard; adjust if your model uses a different demographic ceiling.
- This `pyramid` object is what you pass to `BirthsByCBR(model, birthrates=..., pyramid=pyramid)` to set the age distribution of newborns.

With 1000 samples the observed proportions will be close to targets but will have ~1–3% sampling noise — that's expected. At 10,000+ samples the match tightens significantly.
