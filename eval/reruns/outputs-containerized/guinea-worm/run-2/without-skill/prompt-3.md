All four files are written to `/home/evaluser/`. Here is a full summary of the design.

---

## File layout

```
seis_components.py   – SEIS disease state + dog reservoir
interventions.py     – 4 prevention-based intervention components
guinea_worm_model.py – Model subclass that assembles everything
run_simulation.py    – Example scenario, params, runner, and report
```

---

## SEIS dynamics (`seis_components.py`)

The defining biological constraint — no lasting immunity — is encoded in `GuineaWormInfectious.step()`:

```python
# I → S: worm fully emerged; host returns to Susceptible (SEIS — no immunity)
people.state[trans_mask]  = State.SUSCEPTIBLE.value
model.nodes.S[tick + 1, :] += delta   # back to S, NOT R
```

**State machine per agent (human)**:

```
S ──(WaterTransmission, copepod ingestion)──► E
                                              │ etimer ~365 days (GuineaWormExposed)
                                              ▼
                                              I
                                              │ itimer ~30 days  (GuineaWormInfectious)
                                              ▼
                                              S  ← returns to susceptible, no R compartment
```

**`DogReservoir`** runs the same S→E→I→S cycle at node level using daily binomial draws, sharing the `water_contamination` scalar computed by `WaterTransmission`.

---

## Four intervention components (`interventions.py`)

Each writes a scalar modifier array; `WaterTransmission` reads all of them when building `W_i`:

| Component | Property written | Mechanism |
|---|---|---|
| `AbateLarvicide` | `nodes.abate_factor[:]` | `1 − (coverage × 0.80)` in season; 1.0 outside |
| `WaterFilter` | `people.has_filter[:]` | Assigns at init; `WaterTransmission` scales `λ_j × (1 − 0.95)` for owners |
| `CaseContainment` | `nodes.human_contamination_factor[:]` | Detects 70 % of I agents; reduces their water contribution by 90 % |
| `DogTethering` | `nodes.dog_contamination_factor[:]` | Flags 50 % of villages; `dcf_i = 1 − 0.85 = 0.15` when active |

The contamination formula inside `WaterTransmission.step()`:

```python
W_i = (I_human_i × hcf_i  +  I_dog_i × dcf_i) / (N_human_i + N_dog_i)
W_i *= abate_factor_i
λ_i = beta × W_i
λ_j = λ_i × (1 − 0.95)  if agent has filter,  else λ_i
```

---

## Component step order

```
CaseContainment      → human_contamination_factor
DogTethering         → dog_contamination_factor
AbateLarvicide       → abate_factor
WaterTransmission    → water_contamination  +  S→E transitions
DogReservoir         → dog S/E/I (reads water_contamination)
GuineaWormExposed    → E→I  (etimer countdown)
GuineaWormInfectious → I→S  (itimer countdown, SEIS return)
```

Interventions always run before `WaterTransmission` so the force-of-infection calculation sees current-tick modifiers. `DogReservoir` always runs after `WaterTransmission` so it uses the freshly computed contamination level.
