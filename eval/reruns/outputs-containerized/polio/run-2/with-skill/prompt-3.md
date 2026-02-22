Now I have everything I need. The components already exist in `custom_components.py` and are used in `polio_seir_10patch.py`. Let me explain the complete setup clearly, including the key architectural decisions.

---

## Vaccination Setup for the LASER Polio Model

### Why built-in immunization classes aren't used

The LASER API provides three immunization classes. None can satisfy the requirement to track V separately from R:

| Class | Sets | Effect on transmission | Usable here? |
|---|---|---|---|
| `ImmunizationCampaign` | `susceptibility = 0` | **None** — `TransmissionSE` only checks `state == SUSCEPTIBLE(0)`, ignores `susceptibility` | No |
| `RoutineImmunization` | `susceptibility = 0` | **None** — same bug | No |
| `RoutineImmunizationEx` | `state = RECOVERED(3)` | Correct — blocks infection | No — merges V into R, loses OPV/natural distinction. Also requires Numba callables, no SIA analog |

All three must be skipped. The V compartment distinction (OPV wanes → S, natural recovery is permanent) requires custom components.

---

### Component 1: `VaccinatedCompartment`

Manages the V state and the per-agent OPV waning timer. Must come **before** `BirthsByCBR` and the vaccination component in the component list.

```python
# In custom_components.py
VACCINATED = np.int8(4)  # Distinct from RECOVERED = 3

class VaccinatedCompartment:
    """Tracks vaccinated agents with waning OPV immunity (V→S).

    - Adds people.vtimer (int16): per-agent countdown in days
    - Adds nodes.V and nodes.newly_vax_waned (per-tick arrays)
    - Each tick: decrements vtimer for VACCINATED agents;
      when vtimer hits 0, agent returns to SUSCEPTIBLE
    """
    def __init__(self, model, wandurdist, wandurmin=365):
        ...
        model.people.add_scalar_property("vtimer", dtype=np.int16, default=0)
        # nodes.V is created by Model(additional_states=["V"]) — see model init below
```

Waning duration for OPV type 1 — gamma with mean ~3 years:

```python
vax_waning_dist = dists.gamma(shape=3, scale=365.0)  # mean = 3×365 = 1095 days
vax_compartment = VaccinatedCompartment(model, vax_waning_dist, wandurmin=365)
```

---

### Component 2: `PerPatchVaccinationSEIRV`

Implements both RI and SIA in a single component. Uses correlated missedness — each agent receives a permanent `reachable` flag at birth (via `on_birth` callback), so the same children are missed by every campaign round.

**Routine Immunization logic** (`scripts/custom_components.py:224-244`):
- Fires every `ri_period=7` days
- Finds susceptible, reachable agents whose age falls in `[ri_age, ri_age + ri_period)` — the 7-day window catches everyone turning 42 days old that week
- Draws from per-district `ri_coverage` → moves S → VACCINATED, assigns vtimer

**SIA logic** (`scripts/custom_components.py:246-287`):
- Fires when `tick % sia_period == 0`
- Susceptible children 0–5 years → VACCINATED (new vaccination)
- Already-VACCINATED children → vtimer reset (booster, renews waning clock)
- Reachable agents only, per-district `sia_coverage` draw

```python
vaccination = PerPatchVaccinationSEIRV(
    model,
    vaccinated_compartment=vax_compartment,
    unreachable_frac=np.array(scenario.unreachable_frac, dtype=np.float32),
    ri_coverage=0.80,    # 80% among reachable agents
    sia_coverage=0.90,   # 90% among reachable agents
    ri_age=42,           # 6 weeks = 42 days
    ri_period=7,         # Weekly RI check
    sia_period=180,      # Every 6 months
    sia_max_age=1825,    # 0–5 years (5×365)
)
```

---

### Complete model assembly

Three non-obvious requirements beyond the component list itself:

**1. `additional_states=["V"]`** — registers V with the node frame and auto-copies `V[t+1] = V[t]` each tick (same as S, E, I, R):

```python
model = Model(scenario, PARAMS, birthrates=birthrate_array,
              additional_states=["V"])
```

**2. `MortalityByCDR` must map the V state** — otherwise deaths of vaccinated agents don't decrement `nodes.V`, breaking `S+E+I+R+V = N`:

```python
MortalityByCDR(model, mortalityrates=mortality_array, mappings=[
    (SEIR.State.SUSCEPTIBLE.value, "S"),
    (SEIR.State.EXPOSED.value,     "E"),
    (SEIR.State.INFECTIOUS.value,  "I"),
    (SEIR.State.RECOVERED.value,   "R"),
    (VACCINATED,                   "V"),   # <-- required, non-default
])
```

**3. `initialize_ages(model, pyramid)` before `model.run()`** — without this, every agent has `dob=0`, creating a massive artificial vaccination spike at day 42 when RI first fires:

```python
# Set DOBs from age pyramid AFTER components are set (BirthsByCBR creates dob array)
# but BEFORE model.run()
count = model.people.count
ages_years = pyramid.sample(count=count, dtype=np.int32)
ages_days = np.minimum(ages_years, 99) * 365 + np.random.randint(0, 365, count)
model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)
```

**Full component list with ordering rationale** (`scripts/polio_seir_10patch.py:262-280`):

```python
model.components = [
    SEIR.Susceptible(model),                              # 1. Propagates S[t+1]=S[t]
    SEIR.Exposed(model, expdurdist, infdurdist),          # 2. E→I (decrements etimer)
    SEIR.Infectious(model, infdurdist),                   # 3. I→R (decrements itimer)
    SEIR.Recovered(model),                                # 4. Propagates R[t+1]=R[t]
    vax_compartment,                                      # 5. V waning: V→S (decrements vtimer)
    BirthsByCBR(model, birthrates=birthrate_array,        # 6. Births → triggers on_birth
                pyramid=pyramid),                         #    (sets reachable flag on newborns)
    MortalityByCDR(model, mortalityrates=mortality_array, # 7. Deaths (V mapping required)
                   mappings=[...]),
    vaccination,                                          # 8. RI (S→V at 6wk) + SIA (S→V, V boost)
    PatchImportation(model, infdurdist,                   # 9. Endemic corridor seeding
                     endemic_patch_indices, period=30, count=3),
    SEIR.Transmission(model, expdurdist,                  # 10. Force of infection + Bernoulli trials
                      seasonality=seasonality),
]
```

**Why this ordering matters:**
- `Susceptible` before `VaccinatedCompartment`: V→S waning must find stable S counts to increment
- `BirthsByCBR` after `VaccinatedCompartment`: the `on_birth` callback on `vaccination` sets the `reachable` flag on newborns — `BirthsByCBR` must be instantiated before `vaccination` is added to components, but `vaccination` must appear after `BirthsByCBR` in the list so births are committed before RI checks that tick
- `Transmission` last: operates on the S counts left after all other modifications that tick

---

### State transitions summary

```
Birth → S
S  →(RI at day 42, 80% of reachable)→ V  [vtimer ~ Gamma(3, 365), mean 3yr]
S  →(SIA, 90% of reachable, 0-5yr)→   V  [vtimer assigned]
V  →(SIA booster)→                     V  [vtimer reset]
V  →(vtimer = 0)→                      S  [waned, susceptible again]
S  →(force of infection)→              E → I → R  [permanent natural immunity]
```

---

### Outputs tracked

```python
model.nodes.V[tick]               # vaccinated count per patch per tick
model.nodes.R[tick]               # naturally recovered per patch per tick
model.nodes.ri_doses_given[tick]  # RI doses per patch per tick
model.nodes.sia_doses_given[tick] # SIA doses (new + booster) per patch per tick
model.nodes.newly_vax_waned[tick] # agents returning S→V per patch per tick
```

The complete implementation is in `scripts/custom_components.py` (both components) and `scripts/polio_seir_10patch.py` (model assembly at lines 208–280).
