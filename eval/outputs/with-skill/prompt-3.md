The model already has all three vaccination features fully implemented. Here's a walkthrough of the complete component setup:

## What's already in place

### 1. Custom state value (`custom_components.py:31`)

```python
VACCINATED = np.int8(4)  # distinct from RECOVERED=3
```

LASER's built-in `TransmissionSE` only infects agents with `state == SUSCEPTIBLE (0)`, so any non-zero state blocks infection. Using value 4 (instead of reusing 3) lets us track vaccine-derived vs natural immunity separately.

### 2. `VaccinatedCompartment` — waning OPV immunity (`custom_components.py:34-107`)

This component handles the V compartment lifecycle:

- Adds a per-agent `vtimer` (int16 countdown) and node-level `V[t]` array
- Each tick: decrements `vtimer` for all `state == VACCINATED` agents
- When `vtimer` hits 0: agent returns to `SUSCEPTIBLE`, `nodes.V` decreases, `nodes.S` increases
- Waning duration sampled from a gamma distribution (mean ~3 years for OPV type 1)

**Why custom:** LASER's built-in `RoutineImmunization` and `ImmunizationCampaign` only set `susceptibility = 0`, which `TransmissionSE` ignores. `RoutineImmunizationEx` sets state correctly but doesn't support waning or a separate V compartment.

### 3. `PerPatchVaccinationSEIRV` — RI + SIA (`custom_components.py:110-284`)

This handles both vaccination programs in a single component:

**Routine Immunization (RI):**
- Checked every 7 days (`ri_period=7`)
- Targets susceptible, reachable agents aged exactly 42 days (6 weeks)
- 80% coverage among reachable agents
- Sets `state = VACCINATED (4)` and assigns a waning timer

**SIAs:**
- Every 180 days (`sia_period=180`)
- Targets children 0-1825 days old (0-5 years)
- 90% coverage among reachable agents
- Susceptible children → `VACCINATED` (new vaccination)
- Already-vaccinated children → timer reset (booster effect)

**Correlated missedness:**
- Each agent gets a permanent `reachable` flag at birth based on district-level `unreachable_frac`
- Unreachable agents are never vaccinated by RI or SIA — the same children are missed every time
- This prevents the independent-Bernoulli trap where SIAs appear to cumulatively eliminate immunity gaps

### 4. Model assembly (`polio_seir_10patch.py:230-280`)

The complete component wiring:

```python
# OPV waning: gamma(shape=3, scale=365) → mean 1095 days (~3 years)
vax_waning_dist = dists.gamma(shape=3, scale=365.0)

vax_compartment = VaccinatedCompartment(model, vax_waning_dist, wandurmin=365)
vaccination = PerPatchVaccinationSEIRV(
    model, vax_compartment, unreachable_frac,
    ri_coverage=0.80,       # 80% RI coverage (OPV at 6 weeks)
    sia_coverage=0.90,      # 90% SIA coverage
    ri_age=42,              # 6 weeks = 42 days
    sia_period=180,          # Every 6 months
    sia_max_age=1825,        # 0-5 years
)

model.components = [
    susceptible,          # 1. Count S
    exposed,              # 2. E→I transitions
    infectious,           # 3. I→R transitions
    recovered,            # 4. Count R (natural immunity — permanent)
    vax_compartment,      # 5. V waning (V→S, ~3yr OPV decay)
    BirthsByCBR(...),     # 6. Births
    MortalityByCDR(..., mappings=[
        ...,
        (VACCINATED, "V"),  # Deaths from V reduce V count
    ]),
    vaccination,          # 8. RI + SIA (S→V)
    PatchImportation(...),# 9. Endemic importation
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
]
```

Key details:
- `Model(scenario, PARAMS, additional_states=["V"])` allocates the node-level `V` array and auto-copies `V[t+1] = V[t]` each tick
- `MortalityByCDR` includes a `(VACCINATED, "V")` mapping so deaths of vaccinated agents correctly decrement `nodes.V`
- Component order matters: `VaccinatedCompartment` (waning) runs before vaccination so newly waned agents can be re-vaccinated in the same tick

### Immunity distinction summary

| Property | OPV (Vaccinated) | Natural infection (Recovered) |
|---|---|---|
| State value | `4` (VACCINATED) | `3` (RECOVERED) |
| Node counter | `nodes.V[t]` | `nodes.R[t]` |
| Duration | ~3 years (gamma, wanes to S) | Permanent |
| Booster | SIA resets waning timer | N/A |
| Blocks infection? | Yes (`state != 0`) | Yes (`state != 0`) |

Everything is already wired up and running. The summary output at `polio_seir_10patch.py:458` reports both `V/N` and `R/N` fractions per district so you can verify the two pools are tracked independently.
