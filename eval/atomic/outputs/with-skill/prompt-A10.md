Here's the correct setup with the SEIR components and Transmission, and why the order matters:

```python
import laser.core.distributions as dists
from laser.generic import SEIR, Model

# Duration distributions for exposed and infectious periods
expdurdist = dists.gamma(shape=40, scale=0.25)   # ~10 day latent period
infdurdist = dists.normal(loc=8, scale=2)         # ~8 day infectious period

model = Model(scenario, parameters)

model.components = [
    SEIR.Susceptible(model),                                    # 1. Propagate S counts
    SEIR.Exposed(model, expdurdist, infdurdist),                # 2. E→I transitions
    SEIR.Infectious(model, infdurdist),                         # 3. I→R transitions
    SEIR.Recovered(model),                                      # 4. Propagate R counts
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),  # 5. S→E transitions
]
```

## Why the order matters

Components execute **in list order every tick**. The ordering above follows LASER's convention:

1. **`Susceptible` first** — Copies the S count from the previous tick forward (`S[t] = S[t-1]`). This establishes the baseline susceptible pool before any transitions modify it.

2. **`Exposed` second** — Processes E→I transitions: agents whose latent period expires move to Infectious. This decrements `E[t]` and increments `I[t]`. It needs to run before Transmission so that the current tick's exposed count reflects completed transitions, not new infections.

3. **`Infectious` third** — Processes I→R transitions: agents whose infectious period expires move to Recovered. This decrements `I[t]` and increments `R[t]`. Running this before Transmission ensures the force of infection is computed from the *current* number of infectious agents (after recoveries), not an inflated count.

4. **`Recovered` fourth** — Copies the R count forward (`R[t] = R[t-1]`), then absorbs increments from I→R transitions above. Like Susceptible, it's a bookkeeping step.

5. **`Transmission` last** — Computes the force of infection using the now-accurate `I[t]` and `S[t]` counts, then moves newly-infected susceptibles from S→E. Running it last ensures it sees the correctly updated compartment counts after all other transitions have resolved.

**The key invariant**: `S + E + I + R = N` must hold at every tick. `Susceptible` and `Recovered` act as bookkeeping "bookends" that propagate counts forward, while `Exposed` and `Infectious` handle the internal transitions. `Transmission` comes last because it needs accurate compartment counts to compute the force of infection — if it ran before I→R transitions, it would overestimate infectiousness; if it ran before E→I, it would miscalculate the exposed pool.
