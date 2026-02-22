The intervention components and model are already fully implemented in this codebase. Here's a precise walkthrough of the complete design.

---

## What's already in the codebase

**`scripts/guinea_worm_components.py`** — all four intervention components plus two intervention-aware transmission components  
**`scripts/guinea_worm_chad_interventions.py`** — 10-year model using them

The baseline model (`guinea_worm_chad.py`) uses `SEIR.Transmission` + `DualHostWaterTransmission`. The intervention model replaces both with the intervention-aware variants.

---

## Component design

### 1. `ABATELarvicide` — patch-level FOI multiplier

`guinea_worm_components.py:362`

```python
# Initialization adds a vector property to nodes
model.nodes.add_vector_property("abate_factor", nticks + 1, dtype=np.float32)
model.nodes.abate_factor[:] = 1.0

# step() sets next tick's factor
if season_start <= day_of_year <= season_end:
    nodes.abate_factor[tick + 1] = 1.0 - coverage * efficacy  # e.g. 0.50 * 0.80 = 0.40 residual
else:
    nodes.abate_factor[tick + 1] = 1.0
```

Applied to **both** human and dog models — same water sources. Coverage varies per-patch via a 1D array.

---

### 2. `WaterFilterDistribution` — per-agent ingestion risk reduction

`guinea_worm_components.py:421`

```python
# Adds scalar property on people
model.people.add_scalar_property("has_filter", dtype=np.int8, default=0)

# _assign_filters() uses per-patch adoption rates
rates = self.adoption_rate[nodeids]
model.people.has_filter[istart:iend] = (random(n) < rates).astype(np.int8)

# on_birth() called by BirthsByCBR — newborns inherit household rate
def on_birth(self, istart, iend, tick):
    self._assign_filters(istart, iend)

# step() is a no-op — protection applied inside IntervenedWaterTransmission
```

Human model only. Dogs don't use filters.

---

### 3. `CaseContainment` — contained vs. escaped infectious humans

`guinea_worm_components.py:477`

Three-state flag on agents: 0 = pending, 1 = contained, 2 = uncontained.

```python
# Reset for agents who left INFECTIOUS state (SEIS cycle — they can be
# reinfected and need a fresh containment decision)
not_infectious = (people.state[:count] != SEIR.State.INFECTIOUS.value)
people.contained[:count][not_infectious] = 0

# One-time Bernoulli decision for newly infectious agents
pending = np.nonzero(
    (people.state[:count] == SEIR.State.INFECTIOUS.value) &
    (people.contained[:count] == 0)
)[0]
draws = np.random.random(len(pending))
people.contained[pending[draws < detection_rate]] = 1   # contained
people.contained[pending[draws >= detection_rate]] = 2  # escaped
```

The SEIS-reset is critical: without it, a re-infected agent who was previously contained would skip re-evaluation.

---

### 4. `DogTethering` — mirrors CaseContainment for dogs

`guinea_worm_components.py:553`

Identical structure, but `village_coverage` replaces `detection_rate` and the property is named `tethered` instead of `contained`. Dog model only.

---

## Intervention-aware transmission

### `IntervenedWaterTransmission` — replaces `SEIR.Transmission`

`guinea_worm_components.py:627`

The FOI pipeline with all four modifiers:

```python
# 1. Effective I per patch (containment/tethering reduce contribution)
contributions = np.ones(len(infectious_idx))
if hasattr(people, "contained"):
    contributions[contained_mask] *= (1.0 - containment_efficacy)  # 0.10 residual
if hasattr(people, "tethered"):
    contributions[tethered_mask] *= (1.0 - tether_efficacy)        # 0.15 residual

I_eff = np.zeros(nnodes)
np.add.at(I_eff, people.nodeid[infectious_idx], contributions)

# 2. Gravity-coupled effective prevalence
prev = I_eff / N
for i in range(nnodes):
    local_frac = 1.0 - network[i].sum()
    coupled_prev[i] = local_frac * prev[i] + dot(network[i], prev)

# 3. Patch-level FOI (ABATE applied here)
abate = nodes.abate_factor[tick + 1]  # set by ABATELarvicide this tick
foi = beta * season * coupled_prev * abate

# 4. Per-agent probability (filters applied here)
p_inf = 1.0 - exp(-foi[agent_patch])
if has_filter:
    p_inf *= (1.0 - filter_efficacy)  # 0.05 residual
```

`hasattr` checks let the same class work for human (containment + filters + ABATE) and dog (tethering + ABATE) models without subclassing.

### `IntervenedDualHostTransmission` — replaces `DualHostWaterTransmission`

`guinea_worm_components.py:782`

Same logic, but `_compute_effective_prevalence()` accounts for containment/tethering in the **source** species before computing cross-species FOI. ABATE and filters then apply on the **target** side.

---

## Component ordering in each model

```python
human_model.components = [
    SEIR.Susceptible(human_model),           # initializes state/nodeid
    SEIR.Exposed(human_model, ...),
    SEIR.Infectious(human_model, ...),
    SEIR.Recovered(human_model),             # counts R at tick+1
    SEISRecovery(human_model),               # immediately flips R→S (SEIS)
    BirthsByCBR(...),
    MortalityByCDR(...),
    ABATELarvicide(human_model, ...),        # writes abate_factor[tick+1]
    WaterFilterDistribution(human_model, ...), # passive; on_birth assigns filters
    CaseContainment(human_model, ...),       # decides contained/escaped
    IntervenedWaterTransmission(human_model, ...), # reads all of the above
]

dog_model.components = [
    SEIR.Susceptible, SEIR.Exposed, SEIR.Infectious,
    SEIR.Recovered, SEISRecovery,
    BirthsByCBR, MortalityByCDR,
    ABATELarvicide(dog_model, ...),          # same water sources as humans
    DogTethering(dog_model, ...),
    IntervenedWaterTransmission(dog_model, ...),
]

# After both component lists:
cross_host = IntervenedDualHostTransmission(human_model, dog_model, ...)
cross_host.step(tick)  # called in the manual tick loop
```

Ordering matters: `ABATELarvicide` and `CaseContainment`/`DogTethering` must run **before** `IntervenedWaterTransmission` within the same tick so the multipliers are set when transmission reads them.

---

## Time-varying coverage

`guinea_worm_chad_interventions.py:137`

Coverages ramp annually in the tick loop:

```python
if tick > 0 and tick % 365 == 0:
    year = tick // 365
    cov = get_intervention_coverages(year)   # e.g. ABATE 30%→90% over 7 yrs

    h_abate.coverage[:] = cov["abate"]       # patch array update in-place
    d_abate.coverage[:] = cov["abate"]
    update_filter_coverage(human_model, h_filters, cov["filter"])  # distribute new filters
    h_containment.detection_rate = cov["containment"]
    d_tethering.village_coverage = cov["tethering"]
```

`update_filter_coverage()` (`line 232`) finds agents currently without filters and gives them out to reach the new national rate, then updates `filter_component.adoption_rate` for future newborns.
