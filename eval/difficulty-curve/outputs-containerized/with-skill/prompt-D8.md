The code is written to `/tmp/eastland_respiratory_model.py` (438 lines). Here is the complete model with all key design decisions explained:

---

## `eastland_respiratory_model.py`

### Architecture at a glance

| Section | Decision |
|---|---|
| **Beta** | `R0 / mean_inf_period = 5 / (shape×scale) = 5/10 = 0.5` |
| **Distributions** | `expdurdist = dists.gamma(4, 1.0)`, `infdurdist = dists.gamma(5, 2.0)` |
| **Seasonal forcing** | Cosine peaked at day 15 (mid-Jan), amplitude 0.3 → `[0.7, 1.3]`, normalized to mean=1.0, tiled into `ValuesMap.from_timeseries` |
| **Gravity network** | **Critical fix:** `dist_matrix` initialized with `np.full(..., np.inf)` on the diagonal — a zero diagonal would give `k·p²/0 = inf`, which `row_normalizer` turns into `nan`, silently zeroing the network |
| **Vital dynamics** | `birthrates` shape `(NTICKS, nnodes)` = 30.0 per-1000/yr (asserted in range 1–60), same array passed to both `Model()` and `BirthsByCBR` |
| **Component order** | `Susceptible → Exposed → Infectious → Recovered → PulseImportation → Transmission → BirthsByCBR → MortalityByCDR` |

### `PulseImportation` component

```python
class PulseImportation:
    def __init__(self, model, infdurdist, period=60, count=3):
        self.model = model; self.infdurdist = infdurdist
        self.period = period; self.count = count

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0:
            return
        # find susceptibles in patch 0, infect min(3, available)
        # set state=INFECTIOUS, sample itimer from infdurdist
        # update nodes.S[tick+1, 0] -= n_infect; nodes.I[tick+1, 0] += n_infect
```

Fires at ticks 60, 120, 180, … (every 60 days). Targets susceptibles only, transitions directly to `INFECTIOUS` (same pattern as the reference `Importation` in `custom_components.py`), and updates node-level counts for the next tick slot.

### Key bug fixed vs. the pre-existing draft

```python
# WRONG (pre-existing): diagonal = 0 → gravity gives inf → row_normalizer gives nan
dist_matrix = np.zeros((nnodes, nnodes))

# CORRECT: diagonal = inf → gravity gives 0 → row_normalizer is clean
dist_matrix = np.full((nnodes, nnodes), np.inf)
```

Run with:
```bash
python /tmp/eastland_respiratory_model.py
```
