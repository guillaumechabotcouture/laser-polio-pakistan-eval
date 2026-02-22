The vaccination system is already fully implemented in the codebase. Let me walk through the complete design, explaining what LASER's built-in classes provide versus what required custom components.

---

## Why custom components are needed

LASER provides two built-in immunization classes, but both fail for this use case:

- **`ImmunizationCampaign` / `RoutineImmunization`**: set `susceptibility = 0`, but `TransmissionSE` only checks `state == SUSCEPTIBLE`. **Zero effect on transmission.** (Documented in the skill as a known silent failure.)
- **`RoutineImmunizationEx`**: sets `state = RECOVERED` — correct for transmission, but merges vaccine-derived and natural immunity into the same compartment. No OPV waning.

For a polio model distinguishing OPV (mucosal, waning) from natural infection (strong, long-lasting), custom components are required for all three vaccination tasks.

---

## Component 1: `VaccinatedCompartment`

**File:** `scripts/custom_components.py:37`

Manages the V compartment with per-agent waning timers.

```python
# State values
VACCINATED = np.int8(4)   # distinct from RECOVERED = 3
# state=4 → blocked from infection (TransmissionSE checks state==0)
# state=3 → also blocked, but permanent (natural immunity)
```

**Setup:**
```python
# OPV mucosal immunity: gamma(shape=3, scale=365) → mean ~3 years
vax_waning_dist = dists.gamma(shape=3, scale=365.0)
vax_compartment = VaccinatedCompartment(model, vax_waning_dist, wandurmin=365)
```

What it adds to the model:
- `model.people.vtimer` — `int16` countdown per agent (days until V→S)
- `model.nodes.V` — node-level V count (via `additional_states=["V"]` on Model construction)
- `model.nodes.newly_vax_waned` — tracking array

**`step(tick)` logic** (`custom_components.py:78`):
1. Find all agents where `state == VACCINATED` and `vtimer > 0`
2. Decrement timers by 1
3. Agents where `vtimer` just hit 0: `state → SUSCEPTIBLE`, update `nodes.V` and `nodes.S`

Natural immunity (state=RECOVERED) is never touched — `VaccinatedCompartment` only processes state=4.

---

## Component 2: `PerPatchVaccinationSEIRV`

**File:** `scripts/custom_components.py:113`

Implements both RI and SIA with correlated missedness.

**Setup:**
```python
vaccination = PerPatchVaccinationSEIRV(
    model, vax_compartment, unreachable_frac,
    ri_coverage=0.80,    # 80% among reachable agents
    sia_coverage=0.90,   # 90% among reachable agents
    ri_age=42,           # 6 weeks = 42 days
    sia_period=180,      # every 6 months
    sia_max_age=1825,    # 0–5 years (5*365)
)
```

### Correlated missedness (`custom_components.py:183`)

Each agent receives a permanent `reachable` flag (0/1) at birth, based on district-level `unreachable_frac`. Agents with `reachable=0` are never vaccinated by any campaign.

```python
# At birth (or initialization):
draws = np.random.random(n)
people.reachable[istart:iend] = (draws >= self.unreachable_frac[nodeids])
```

This is critical for Pakistan: districts like N. Waziristan (40% unreachable) and Quetta (35%) have structural immunity gaps that independent Bernoulli draws would overestimate as filling over multiple SIA rounds. With correlated missedness, the same children are missed every time.

### Routine Immunization (`custom_components.py:224`)

Runs every `ri_period=7` days. Eligible agents must be:
- `state == SUSCEPTIBLE`
- `reachable == 1`
- Age between 42 and 49 days (6-week window, checked weekly to avoid double-dosing)

```python
ages = tick - people.dob[:count]
ri_eligible = np.nonzero(
    (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
    (people.reachable[:count] == 1) &
    (ages >= 42) & (ages < 49)
)[0]
# Stochastic draw against per-district coverage
draws = np.random.random(len(ri_eligible))
vaccinated = ri_eligible[draws < self.ri_coverage[people.nodeid[ri_eligible]]]
```

On vaccination: `state → VACCINATED(4)`, `vtimer` set from `vax_compartment.sample_waning_duration()`, node `S--`, `V++`.

### SIA Campaigns (`custom_components.py:246`)

Fires when `tick % 180 == 0`. Two sub-operations for ages 0–5, reachable only:

1. **Susceptible children** (`state == SUSCEPTIBLE`): `S → V`, assign waning timer
2. **Already-vaccinated children** (`state == VACCINATED`): reset `vtimer` (booster — extends immunity without double-counting)

Exposed (E), Infectious (I), and naturally Recovered (R) children are not targeted — SIAs don't provide additional protection beyond existing immunity, and natural immunity is stronger.

---

## Complete component assembly

**File:** `scripts/polio_seir_10patch.py:208`

The model must be constructed with `additional_states=["V"]` so LASER's internal state propagation includes V:

```python
model = Model(scenario, PARAMS, birthrates=birthrate_array,
              additional_states=["V"])
```

Component order (`polio_seir_10patch.py:262`):

```python
model.components = [
    SEIR.Susceptible(model),           # 1. Propagate S[t+1]=S[t], register nodeid/state
    SEIR.Exposed(model, ...),          # 2. E→I transitions (etimer countdown)
    SEIR.Infectious(model, ...),       # 3. I→R transitions (itimer countdown)
    SEIR.Recovered(model),             # 4. Propagate R[t+1]=R[t] — natural immunity
    vax_compartment,                   # 5. V waning: vtimer--, V→S when expired
    BirthsByCBR(model, ...),           # 6. New births → S, triggers on_birth callbacks
    MortalityByCDR(model, ...),        # 7. Deaths from all states including V
    vaccination,                       # 8. RI (weekly) + SIA (biannual): S→V
    PatchImportation(model, ...),      # 9. Seed infections in endemic corridors
    SEIR.Transmission(model, ...),     # 10. FOI: S→E based on current I and network
]
```

**Ordering rationale:**
- `vax_compartment` (5) runs before `BirthsByCBR` (6) so waning timers are decremented before newborns are added — newborns are born into S, not counted in this tick's waning
- `vaccination` (8) runs after births/deaths so the eligible population is current
- `Transmission` (10) runs last so it sees the post-vaccination susceptible pool

**MortalityByCDR mapping** must include V (`custom_components.py:269`):
```python
MortalityByCDR(model, mortalityrates=mortality_array, mappings=[
    (SEIR.State.SUSCEPTIBLE.value, "S"),
    (SEIR.State.EXPOSED.value,     "E"),
    (SEIR.State.INFECTIOUS.value,  "I"),
    (SEIR.State.RECOVERED.value,   "R"),
    (VACCINATED, "V"),              # must be explicit — not a standard SEIR state
])
```

---

## Compartment semantics summary

| State | Value | Transitions | Immunity |
|-------|-------|-------------|----------|
| S | 0 | → E (infection), → V (vaccination) | None |
| E | 1 | → I (etimer) | None |
| I | 2 | → R (itimer) | None |
| **R** | **3** | Permanent | **Natural: strong, lifelong** |
| **V** | **4** | → S (vtimer=0) | **OPV: mucosal, ~3yr mean** |

Both R and V are protected from infection (TransmissionSE checks `state == 0`), but only V can revert to S. The summary statistics at `polio_seir_10patch.py:458` report S/N, V/N, and R/N separately to verify the compartments are tracking correctly.
