The complete model is at `/home/evaluser/guinea_worm_water_transmission.py` and passes all 11 verification checks.

---

## What was built and why

### The gap in the existing model

`DualHostWaterTransmission` in `guinea_worm_components.py` computes cross-species FOI using **local prevalence only** — dogs in patch j infect only humans in patch j. The three requested features require that dogs in district j contaminate water bodies shared with humans in neighboring district i — the cross-species coupling needs to be spatially explicit.

### New component: `WaterSourceCrossHostTransmission`

This is the key extension. It computes gravity-coupled cross-species prevalence before applying the FOI:

```
dog_coupled_prev[i] = (1 - Σⱼ w[i,j]) * dog_prev[i]
                    + Σⱼ w[i,j] * dog_prev[j]
```

This is the same formula `SEIR.Transmission` uses for within-species spatial coupling — applied here to dog prevalence as the source for human infection. The shared `water_network` (built with c=2.0, b=0.3) represents the physical water bodies both species access.

### Three distinct networks

| Network | Purpose | `gravity_k` |
|---|---|---|
| `human_net` → `human_model.network` | Within-species human FOI via `SEIR.Transmission` | 0.005 |
| `dog_net` → `dog_model.network` | Within-species dog FOI via `SEIR.Transmission` | 0.003 |
| `water_net` → `WaterSourceCrossHostTransmission` | Cross-species FOI across shared water bodies | 0.005 |

### Verification output

```
Dog rate / Human rate = 22.1x  (dogs: 7,453 infections / 180k pop vs humans: 4,114 / 2.2M)
Peak-tick spread across districts: 355 days (humans), 427 days (dogs) — spatial heterogeneity confirmed
R at final tick = 0 — SEIS dynamics correct (no accumulated immunity)
```
