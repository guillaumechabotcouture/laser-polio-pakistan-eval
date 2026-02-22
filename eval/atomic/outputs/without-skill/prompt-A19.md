Here's the file at `vaccinate.py`. The key insight:

### Why `state`, not `susceptibility`?

LASER's `TransmissionSE` Numba kernel has this inner loop:

```python
for i in nb.prange(len(states)):
    if states[i] == State.SUSCEPTIBLE.value:  # == 0
        # ... attempt infection ...
```

It **only** checks whether `state == 0`. It never reads a `susceptibility` property. So:

| Approach | Effect on TransmissionSE |
|---|---|
| `model.people.susceptibility[idx] = 0.0` | **No effect** — agent is still state 0 and will be infected |
| `model.people.state[idx] = State.RECOVERED.value` | **Correct** — agent is state 3, kernel skips it entirely |

### What the code does

1. **Finds all susceptible agents** — `state == State.SUSCEPTIBLE.value` (0)
2. **Randomly selects 80%** of them with `np.random.default_rng().choice(..., replace=False)`
3. **Sets their state to `RECOVERED` (3)** — this is what makes them invisible to TransmissionSE
4. **Updates node-level S/R counters** with `np.subtract.at` / `np.add.at` — LASER tracks aggregate counts per node per tick, and these must stay consistent or validation breaks

### The common mistake

```python
# WRONG — TransmissionSE ignores this field entirely
model.people.susceptibility[chosen] = 0.0

# CORRECT — TransmissionSE checks state, not susceptibility
model.people.state[chosen] = State.RECOVERED.value
```

The demo at the bottom builds a minimal model with 10,000 agents across 4 nodes, vaccinates 80%, and asserts that both agent states and node-level bookkeeping are consistent.
