The simulation is written to `polio_seir_model.py`. Here's a summary of the model design:

**Model structure:**
- **Framework:** LASER (`laser-generic`) spatial agent-based SEIR
- **10 patches** with real geographic coordinates for Pakistani polio-priority districts (Karachi, Quetta, Peshawar, Lahore, Rawalpindi, Bannu, Tank, D.I. Khan, Chaman, Zhob)
- **Gravity-model spatial coupling** using Haversine inter-district distances

**Parameters:**
| Parameter | Value | Derivation |
|---|---|---|
| `beta` | 0.2143/day | R0 / infectious duration = 6 / 28 |
| Latent period | Gamma(shape=3, scale=1) | Mean = 3 days |
| Infectious period | Gamma(shape=14, scale=2) | Mean = 28 days |
| Duration | 3,650 ticks | 10 years at daily resolution |

**Initial conditions per patch:**
| Compartment | Count | Rationale |
|---|---|---|
| S | 4,995 | Remaining after R and I |
| E | 0 | — |
| I | 5 | Seed cases |
| R | 95,000 | 95% baseline immunity |

**Components registered** (in execution order): `Susceptible` → `Exposed` → `Infectious` → `Recovered` → `Transmission`

**Output:** The script prints a full summary including aggregate totals, per-district final state, annual incidence breakdown, and cumulative incidence by district.

Run with:
```bash
python polio_seir_model.py
```
