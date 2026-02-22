Here's the complete setup in `seir_model.py`. The five components are assigned to `model.components` as a list, and the model executes them **in list order** each tick.

## Execution order and why it matters

```
Tick start → _initialize_flows() copies state(t) → state(t+1)
         ↓
   1. Susceptible    — tallies S counts, initializes agent properties
   2. Transmission   — reads S and I counts → computes force of infection → S → E
   3. Exposed        — decrements incubation timers → E → I when timer hits 0
   4. Infectious     — decrements infectious timers → I → R when timer hits 0
   5. Recovered      — tallies final R counts
```

**Why this order is required:**

1. **Susceptible runs first** because it initializes per-agent properties (`nodeid`, `state`) and produces the node-level S counts that Transmission needs.

2. **Transmission runs second** because it must compute force of infection from the *beginning-of-tick* S and I counts (before any within-tick transitions alter them). It converts S agents to E and sets their incubation timer (`etimer`).

3. **Exposed runs third** — after Transmission has created new E agents, this component decrements `etimer` for all exposed agents. When `etimer` reaches 0, it transitions E → I and sets the infectious timer (`itimer`).

4. **Infectious runs fourth** — after Exposed may have added new I agents, this component decrements `itimer` and transitions I → R when the infectious period ends.

5. **Recovered runs last** because it passively tallies R counts. All transitions feeding into R (from Infectious) must have already occurred for the counts to be correct.

Reversing or shuffling the order would cause components to read stale or not-yet-computed counts, leading to incorrect force-of-infection calculations and broken state accounting.
