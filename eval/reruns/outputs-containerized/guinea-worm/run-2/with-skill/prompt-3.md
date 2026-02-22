The workspace filesystem is mounted read-only — I can't edit files in this environment. Let me instead give you a complete walkthrough of what already exists and the exact changes needed to fix the three spec discrepancies.

---

## What already exists

All four components are implemented in `scripts/guinea_worm_components.py` and wired into `scripts/guinea_worm_chad_interventions.py`. Here's the complete architecture.

---

## Component implementations (`guinea_worm_components.py`)

### 1. `ABATELarvicide` — patch-level copepod reduction

```python
class ABATELarvicide:
    """ABATE (temephos) larvicide applied monthly to known water sources.

    Sets nodes.abate_factor[tick+1] each tick:
        1.0                          — off season
        1 - (coverage * efficacy)    — dry season (Jan–May, days 0-150)

    Must be added to BOTH human and dog models (same water sources).
    """
    def __init__(self, model, coverage_by_patch, efficacy=0.80,
                 season_start_day=0, season_end_day=150):
        nnodes = model.nodes.count
        if np.isscalar(coverage_by_patch):
            self.coverage = np.full(nnodes, coverage_by_patch, dtype=np.float32)
        else:
            self.coverage = np.asarray(coverage_by_patch, dtype=np.float32)

        assert len(self.coverage) == nnodes
        assert np.all((self.coverage >= 0) & (self.coverage <= 1))

        self.efficacy = efficacy          # 0.80 kills 80% of copepods
        self.season_start = season_start_day
        self.season_end = season_end_day

        if not hasattr(model.nodes, "abate_factor"):
            model.nodes.add_vector_property(
                "abate_factor", model.params.nticks + 1, dtype=np.float32
            )
            model.nodes.abate_factor[:] = 1.0

    def step(self, tick):
        day_of_year = tick % 365
        if self.season_start <= day_of_year <= self.season_end:
            self.model.nodes.abate_factor[tick + 1] = \
                1.0 - self.coverage * self.efficacy  # vector broadcast
        else:
            self.model.nodes.abate_factor[tick + 1] = 1.0
```

### 2. `WaterFilterDistribution` — per-agent filter protection

```python
class WaterFilterDistribution:
    """Cloth/pipe water filter distribution to households.

    Assigns people.has_filter (0/1) at init based on per-patch adoption_rate.
    on_birth() initialises newborns. Passive during step() — protection is
    applied inside IntervenedWaterTransmission.

    Only added to the human model (dogs don't use filters).
    """
    def __init__(self, model, adoption_rate=0.60, filter_efficacy=0.95):
        nnodes = model.nodes.count
        if np.isscalar(adoption_rate):
            self.adoption_rate = np.full(nnodes, adoption_rate, dtype=np.float32)
        else:
            self.adoption_rate = np.asarray(adoption_rate, dtype=np.float32)
        self.filter_efficacy = filter_efficacy

        if not hasattr(model.people, "has_filter"):
            model.people.add_scalar_property("has_filter", dtype=np.int8, default=0)
        self._assign_filters(0, model.people.count)

    def _assign_filters(self, istart, iend):
        n = iend - istart
        if n <= 0:
            return
        nodeids = self.model.people.nodeid[istart:iend]
        rates = self.adoption_rate[nodeids]
        draws = np.random.random(n).astype(np.float32)
        self.model.people.has_filter[istart:iend] = (draws < rates).astype(np.int8)

    def on_birth(self, istart, iend, tick):
        self._assign_filters(istart, iend)

    def step(self, tick):
        pass  # protection applied by IntervenedWaterTransmission
```

### 3. `CaseContainment` — human case detection and water exclusion

```python
class CaseContainment:
    """Case containment for human guinea worm infections.

    When an agent enters INFECTIOUS state (worm emerging), it is detected with
    probability detection_rate and flagged as contained (people.contained == 1).
    Contained cases contribute (1 - containment_efficacy) to water FOI instead
    of 1.0.  Status resets when the agent leaves INFECTIOUS (SEIS cycle).

    contained == 0  pending (newly infectious, no decision yet)
    contained == 1  detected and water-access prevented
    contained == 2  escaped surveillance
    """
    def __init__(self, model, detection_rate=0.70, containment_efficacy=0.90):
        self.detection_rate = detection_rate
        self.containment_efficacy = containment_efficacy

        if not hasattr(model.people, "contained"):
            model.people.add_scalar_property("contained", dtype=np.int8, default=0)
        model.nodes.add_vector_property(
            "cases_contained", model.params.nticks + 1, dtype=np.int32
        )

    def on_birth(self, istart, iend, tick):
        self.model.people.contained[istart:iend] = 0

    def step(self, tick):
        people = self.model.people
        count = people.count

        # Reset for agents no longer infectious (SEIS: I→R→S resets flag)
        not_inf = (people.state[:count] != SEIR.State.INFECTIOUS.value)
        people.contained[:count][not_inf] = 0

        # Identify newly infectious with no decision yet
        pending = np.nonzero(
            (people.state[:count] == SEIR.State.INFECTIOUS.value) &
            (people.contained[:count] == 0)
        )[0]
        if len(pending) == 0:
            return

        draws = np.random.random(len(pending))
        detected = pending[draws < self.detection_rate]
        escaped  = pending[draws >= self.detection_rate]
        people.contained[detected] = 1
        people.contained[escaped]  = 2

        if len(detected) > 0:
            self.model.nodes.cases_contained[tick] += np.bincount(
                people.nodeid[detected], minlength=self.model.nodes.count
            ).astype(np.int32)
```

### 4. `DogTethering` — tethering infectious dogs away from water

```python
class DogTethering:
    """Dog tethering/management for infected dogs.

    Mirrors CaseContainment for the dog population. Infectious dogs are
    tethered with probability village_coverage; tethered dogs contribute
    (1 - tether_efficacy) to water FOI.

    tethered == 0  pending, 1 = tethered, 2 = not tethered
    Only added to the dog model.
    """
    def __init__(self, model, village_coverage=0.50, tether_efficacy=0.85):
        self.village_coverage = village_coverage
        self.tether_efficacy = tether_efficacy

        if not hasattr(model.people, "tethered"):
            model.people.add_scalar_property("tethered", dtype=np.int8, default=0)
        model.nodes.add_vector_property(
            "dogs_tethered", model.params.nticks + 1, dtype=np.int32
        )

    def on_birth(self, istart, iend, tick):
        self.model.people.tethered[istart:iend] = 0

    def step(self, tick):
        people = self.model.people
        count = people.count

        not_inf = (people.state[:count] != SEIR.State.INFECTIOUS.value)
        people.tethered[:count][not_inf] = 0

        pending = np.nonzero(
            (people.state[:count] == SEIR.State.INFECTIOUS.value) &
            (people.tethered[:count] == 0)
        )[0]
        if len(pending) == 0:
            return

        draws = np.random.random(len(pending))
        caught = pending[draws < self.village_coverage]
        free   = pending[draws >= self.village_coverage]
        people.tethered[caught] = 1
        people.tethered[free]   = 2

        if len(caught) > 0:
            self.model.nodes.dogs_tethered[tick] += np.bincount(
                people.nodeid[caught], minlength=self.model.nodes.count
            ).astype(np.int32)
```

---

## How FOI is modified (`IntervenedWaterTransmission`)

The within-species transmission replaces `SEIR.Transmission` and composes all four modifiers:

```
I_eff_i   = Σ_k  contribution_k
              where contribution_k = 1.0              if not contained/tethered
                                     1 - efficacy     if contained (human) or tethered (dog)

prev_i    = I_eff_i / N_i
coupled_i = (1 - Σ_j network[i,j]) * prev_i + Σ_j network[i,j] * prev_j
FOI_i     = beta * season(t) * coupled_i * abate_factor_i
P(infect j in patch i) = (1 - exp(-FOI_i)) * (1 - has_filter_j * filter_efficacy)
```

`IntervenedDualHostTransmission` applies the same formula for cross-species contributions (dog→human, human→dog), calling `_compute_effective_prevalence()` on the source model first.

---

## Component assembly (`guinea_worm_chad_interventions.py`)

```python
# ---- Human model ----
h_abate       = ABATELarvicide(human_model, coverage_by_patch=y0["abate"])
h_filters     = WaterFilterDistribution(human_model,
                    adoption_rate=y0["filter"],
                    filter_efficacy=FILTER_EFFICACY)     # 0.95
h_containment = CaseContainment(human_model,
                    detection_rate=y0["containment"],    # 0.70
                    containment_efficacy=CONTAINMENT_EFFICACY)  # 0.90
h_transmission = IntervenedWaterTransmission(
                    human_model, expdurdist,
                    seasonality_values=season_tiled,
                    containment_efficacy=CONTAINMENT_EFFICACY,
                    filter_efficacy=FILTER_EFFICACY)

human_model.components = [
    SEIR.Susceptible(human_model),
    SEIR.Exposed(human_model, expdurdist, infdurdist),
    SEIR.Infectious(human_model, infdurdist),
    SEIR.Recovered(human_model),
    SEISRecovery(human_model),      # I→R→S: no lasting immunity
    BirthsByCBR(human_model, ...),
    MortalityByCDR(human_model, ...),
    h_abate,        # sets nodes.abate_factor each tick
    h_filters,      # passive; on_birth assigns has_filter
    h_containment,  # containment decisions for newly-infectious
    h_transmission, # gravity-coupled FOI with all modifiers
]

# ---- Dog model ----
d_abate    = ABATELarvicide(dog_model, coverage_by_patch=y0["abate"])
d_tethering = DogTethering(dog_model,
                  village_coverage=y0["tethering"],   # 0.50
                  tether_efficacy=TETHER_EFFICACY)    # 0.85
d_transmission = IntervenedWaterTransmission(
                  dog_model, expdurdist,
                  seasonality_values=season_tiled,
                  tether_efficacy=TETHER_EFFICACY)

dog_model.components = [
    SEIR.Susceptible(dog_model),
    SEIR.Exposed(dog_model, expdurdist, infdurdist),
    SEIR.Infectious(dog_model, infdurdist),
    SEIR.Recovered(dog_model),
    SEISRecovery(dog_model),
    BirthsByCBR(dog_model, ...),
    MortalityByCDR(dog_model, ...),
    d_abate,
    d_tethering,
    d_transmission,
]

# ---- Cross-species coupling ----
cross_host = IntervenedDualHostTransmission(
    human_model, dog_model, expdurdist,
    dog_to_human=0.15, human_to_dog=0.05,
    seasonality_values=season_tiled,
    containment_efficacy=CONTAINMENT_EFFICACY,
    tether_efficacy=TETHER_EFFICACY,
    filter_efficacy=FILTER_EFFICACY,
)
```

In the tick loop, cross-species coupling steps **after** both model component lists:

```python
for tick in range(NTICKS):
    human_model._initialize_flows(tick)
    dog_model._initialize_flows(tick)
    for c in human_model.components: c.step(tick)
    for c in dog_model.components:   c.step(tick)
    cross_host.step(tick)           # dog↔human via shared water
```

---

## Three spec discrepancies to fix

The existing `guinea_worm_chad_interventions.py` has these values wrong relative to your spec. Apply the following patch:

```python
# After ENDEMIC_PATCHES — add district-level ABATE baseline:
ABATE_BASE_COVERAGE = np.array(
    [0.80, 0.75, 0.65, 0.55, 0.45, 0.35, 0.30, 0.20],
    dtype=np.float32,
)  # 8 districts, range 20-80% (Chari-Baguirmi highest)

def get_abate_coverage(year):
    return np.minimum(ABATE_BASE_COVERAGE + 0.05 * year, 0.90).astype(np.float32)

def get_intervention_coverages(year):
    return {
        "abate":       get_abate_coverage(year),      # was: scalar 0.30+0.10*year
        "filter":      min(0.60 + 0.05 * year, 0.90), # was: 0.40 (spec says ~60%)
        "containment": min(0.70 + 0.03 * year, 1.0),  # unchanged
        "tethering":   min(0.50 + 0.05 * year, 1.0),  # was: 0.40 (spec says 50%)
    }
```

In the annual report print lines, replace `{cov['abate']:.0%}` with `{cov['abate'].mean():.0%}` since ABATE is now a per-district array. Likewise in `plot_results`, use `c["abate"].mean()` when building the coverage timeline.

---

## Why SEIS doesn't break any of this

`SEISRecovery` runs immediately after `SEIR.Recovered` each tick, moving every newly-recovered agent back to susceptible:

```python
class SEISRecovery:
    def step(self, tick):
        recovered_mask = (people.state[:count] == SEIR.State.RECOVERED.value)
        recovered_idx  = np.nonzero(recovered_mask)[0]
        if len(recovered_idx) == 0:
            return
        people.state[recovered_idx] = SEIR.State.SUSCEPTIBLE.value
        r_by_node = np.bincount(people.nodeid[recovered_idx],
                                minlength=nodes.count).astype(nodes.R.dtype)
        nodes.R[tick + 1] -= r_by_node
        np.maximum(nodes.R[tick + 1], 0, out=nodes.R[tick + 1])
        nodes.S[tick + 1] += r_by_node
```

`CaseContainment` and `DogTethering` both reset their `contained`/`tethered` flags for any agent no longer in INFECTIOUS state, so reinfection after SEIS cycling correctly triggers a fresh surveillance decision.
