Here's how to create and use a gamma distribution with LASER's distribution utilities:

## Creating the Distribution

```python
import laser.core.distributions as dists

# Gamma with mean=10, shape=5, scale=2 (mean = shape × scale = 5 × 2 = 10)
infectious_duration = dists.gamma(shape=5, scale=2)
```

That's it — `dists.gamma(shape, scale)` returns a Numba-compatible callable.

## How It's Used

You then pass this distribution to LASER components that need stochastic durations. In an SEIR model, the two main duration distributions are:

1. **Exposed (latent) period** — time from infection to becoming infectious (E→I)
2. **Infectious period** — time from becoming infectious to recovery (I→R)

For example, to use it as the infectious period distribution:

```python
infdurdist = dists.gamma(shape=5, scale=2)   # mean 10 days
expdurdist = dists.gamma(shape=40, scale=0.25)  # mean 10 days, narrower

model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    SEIR.Transmission(model, expdurdist),
    # ...
]
```

## Sampling

The distribution can be consumed two ways:

```python
# Batch: fill an array of samples
samples = dists.sample_floats(infectious_duration, offsets_array)
int_samples = dists.sample_ints(infectious_duration, offsets_array)

# Inside a Numba @nb.njit loop: call directly with (tick, node)
# value = infectious_duration(tick, node)
```

## Why Gamma for Durations?

The **shape parameter controls the variance** around the mean. With shape=5, scale=2:
- **Mean** = 5 × 2 = 10 days
- **Variance** = 5 × 2² = 20, so **std** ≈ 4.5 days

Compare to the skill's example `dists.gamma(shape=40, scale=0.25)` — same 10-day mean, but shape=40 gives a much tighter distribution (std ≈ 1.6 days). Higher shape → more concentrated around the mean → more realistic for well-characterized disease stages where individual variation is small.

An exponential distribution (shape=1) would be memoryless — biologically unrealistic for most disease stages. The gamma lets you tune how "peaked" the duration distribution is, which affects epidemic dynamics (sharper durations produce sharper epidemic peaks).
