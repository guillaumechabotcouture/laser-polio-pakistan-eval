The model is written to `~/guinea_worm_model.py` (943 lines). The `/eval-work` filesystem is read-only, so this was the appropriate location. The simulation ran successfully and all outputs are consistent with Guinea worm epidemiology.

---

## Model architecture

The file implements all three LASER architectural patterns from scratch using only NumPy and GeoPandas (no laser-core/laser-generic required):

**`PropertySet`** — dot-notation parameter container, mirrors `laser.core.PropertySet`

**`LaserFrame`** — pre-allocated NumPy-backed state store with:
- `add_scalar_property(name)` → 1-D array `[n_nodes]`
- `add_vector_property(name, length)` → 2-D array `[length, n_nodes]` for time-series indexed as `nodes.I_d[tick, node]`

**Components** — each has `step(tick: int)` and is registered in `model.components`:
| Component | Transition | Notes |
|---|---|---|
| `SeasonalForcingComponent` | — | computes φ(t), caches as `model._phi` |
| `WaterTransmissionComponent` | S → E | reads `I[tick]`, writes `tick+1` deltas |
| `ExposedProgressionComponent` | E → I | simultaneous-update: reads `E[tick]` |
| `InfectiousProgressionComponent` | I → S | reads `I[tick]`, prevents same-tick double-transition |

`_initialize_tick(tick)` mirrors LASER's `_initialize_flows`: copies all six state arrays from `tick` to `tick+1` before components apply signed deltas.

---

## Features implemented

**1. Water-source mediated transmission** — effective copepod density at each district's shared water bodies:

```
W_local[i] = 0.95 · I_d[i]/N_d[i]  +  0.05 · I_h[i]/N_h[i]
W_ext[i]   = Σ_j  W[j,i] · W_local[j]     ← gravity spillover
W_eff[i]   = W_local[i] + W_ext[i]
λ_h[i]     = β_h · φ(t) · W_eff[i]
p_SE[i]    = 1 − exp(−λ)               ← Bernoulli
```

**2. Gravity network** (`build_gravity_network`): `W[i,j] = k · pop_i^1.0 · pop_j^0.3 / d_ij^2.0`, column-normalised with configurable inflow cap. Uses vectorised Haversine for great-circle distances.

**3. Seasonal forcing** — monthly multipliers:
- Jan–Mar: **1.8×** (peak; shrinking stagnant pools concentrate copepods)
- Apr: 1.3×, May: 0.7× (declining)
- Jun–Sep: **0.4×** (rainy season; flowing water disperses larvae)
- Oct: 0.7×, Nov: 1.2×, Dec: 1.6× (building back to peak)

**4. Cross-species asymmetric coupling** — dogs contribute 95% of water contamination; both species are exposed to the same `W_eff` but with different `β` values (`beta_d = 0.060`, `beta_h = 0.002`).

---

## Simulation output (10-year run, 9 districts)

- **Population conservation**: S + E + I = N exactly conserved every tick for both species
- **Dog endemic prevalence**: ~4.3% infectious, ~36% incubating (expected given 180-day incubation / 21-day infectious ratio)
- **Human prevalence**: ~0.31% (spill-over driven by dog contamination — consistent with Chad's dog-dominated epidemic)
- **Seasonal R0_dog**: 2.16 (Jan–Mar), 0.48 (Jun–Sep), mean = 1.25 — correctly epidemic in dry season, sub-threshold in rains
- **Spatial coupling**: strongest between geographically proximate districts (Logone-Occidental ↔ Logone-Oriental at W=0.24)
