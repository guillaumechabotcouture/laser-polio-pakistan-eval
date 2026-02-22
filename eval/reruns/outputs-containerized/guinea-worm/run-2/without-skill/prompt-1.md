The file is at `/home/evaluser/guinea_worm_seis.py`. Here is a full walkthrough of every design decision:

---

## Model structure

### SEIS via component composition

`laser-generic` has no built-in SEIS class, but it ships the right building blocks:

| Component | Transition | Role in GW model |
|---|---|---|
| `TransmissionSE` | S → E | Susceptible host drinks copepod-contaminated water; sets `etimer` |
| `Exposed` | E → I | Pre-patent larval maturation over ~365 days; sets `itimer` on emergence |
| **`InfectiousIS`** | **I → S** | Worm expelled/dies; **returns host to S, not R** — this single swap is what makes the model SEIS rather than SEIR |
| `Susceptible` | — | Passive bookkeeping |

### Duration distributions

```
Pre-patent (E timer): gamma(shape=44, scale=8.3)  → mean=365 d, sd≈55 d, min=300 d
Infectious  (I timer): gamma(shape=4,  scale=5.25) → mean=21 d,  sd≈10 d, min=14 d
```

### R₀ calibration

```
β = R0 / mean_infectious_duration
Humans: β = 1.5 / 21 ≈ 0.071
Dogs:   β = 2.0 / 21 ≈ 0.095
```

### Annual case counting

Annual incidence is derived from daily prevalence via **Little's Law**:

```
new_cases/year = mean_daily_I × 365 / mean_infectious_duration
```

This is the number of worm-emergence events (the epidemiologically reported "case" definition for guinea worm).

---

## What the model does NOT yet do (and how to extend)

**Cross-species coupling** — the two models run independently. In reality, dogs and humans draw from the same water source. The comment in section 8 gives the coupling formula:

```python
λ_i = β × (I_human_i + α × I_dog_i) / (N_human_i + N_dog_i)
```

where α ≈ 19 encodes the 95% dog reservoir contribution. Implementing this requires a custom `TransmissionSE`-derived component that reads `model_d.nodes.I` at each tick. The two-model architecture here is already structured to accept that extension with minimal refactoring.

**Seasonal transmission** — Guinea worm peaks in the wet season (May–Oct in Chad). Pass a `seasonality=ValuesMap(...)` argument to `TransmissionSE` to add a sinusoidal forcing term.
