Here's the complete component setup and how to assemble the model with all four interventions.

---

## Components Added to `scripts/guinea_worm_components.py`

### Six new classes, organized in two groups:

**Intervention state components** (maintain per-patch or per-agent state):

| Component | Level | Effect | Model |
|---|---|---|---|
| `ABATELarvicide` | Patch | `abate_factor = 1 - coverage*0.80` during dry season | Both |
| `WaterFilterDistribution` | Agent | `has_filter` flag, 95% risk reduction | Human only |
| `CaseContainment` | Agent | `contained` flag on 70% of cases, 90% FOI reduction | Human only |
| `DogTethering` | Agent | `tethered` flag on 50% of dogs, 85% FOI reduction | Dog only |

**Intervention-aware transmission** (replace built-in transmission):

| Component | Replaces | Role |
|---|---|---|
| `IntervenedWaterTransmission` | `SEIR.Transmission` | Within-species FOI with all 4 modifiers |
| `IntervenedDualHostTransmission` | `DualHostWaterTransmission` | Cross-species FOI with all 4 modifiers |

### How interventions modify the force of infection

```
I_eff = sum over infectious agents:
    1.0                      if not contained/tethered
    (1 - efficacy)           if contained (0.10) or tethered (0.15)

prev = I_eff / N
coupled_prev = gravity-weighted mixing across patches

FOI = beta * season(t) * coupled_prev * abate_factor

P(infect agent j) = (1 - exp(-FOI)) * (1 - has_filter * 0.95)
```

### Complete model assembly with interventions

```python
from guinea_worm_components import (
    SEISRecovery, PatchImportation, build_chad_dry_season,
    # Interventions
    ABATELarvicide, WaterFilterDistribution,
    CaseContainment, DogTethering,
    # Intervention-aware transmission (replaces SEIR.Transmission
    # and DualHostWaterTransmission)
    IntervenedWaterTransmission, IntervenedDualHostTransmission,
)

# --- Per-district ABATE coverage (same water sources for both species) ---
abate_coverage = np.array([
    0.60,  # Chari-Baguirmi (high endemicity)
    0.50,  # Moyen-Chari
    0.40,  # Salamat
    0.30,  # Mandoul
    0.30,  # Logone_Oriental
    0.20,  # Tandjile
    0.20,  # Mayo-Kebbi_Est
    0.80,  # Lac (high dog burden)
], dtype=np.float32)

# --- Seasonality (raw 1D array for custom transmission) ---
_, season_365 = build_chad_dry_season(nticks, nnodes)
season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]

# ==== HUMAN MODEL COMPONENTS ====
human_model.components = [
    # 1. SEIS state machine
    SEIR.Susceptible(human_model),
    SEIR.Exposed(human_model, expdurdist, infdurdist),
    SEIR.Infectious(human_model, infdurdist),
    SEIR.Recovered(human_model),
    SEISRecovery(human_model),            # R→S (no lasting immunity)

    # 2. Interventions (set state BEFORE transmission reads it)
    CaseContainment(human_model,
                    detection_rate=0.70,
                    containment_efficacy=0.90),
    ABATELarvicide(human_model,
                   coverage_by_patch=abate_coverage,
                   efficacy=0.80,
                   season_start_day=0,    # Jan 1
                   season_end_day=150),   # ~May 30
    WaterFilterDistribution(human_model,
                            adoption_rate=0.60,
                            filter_efficacy=0.95),

    # 3. Vital dynamics
    BirthsByCBR(human_model, birthrates=human_birthrates,
                pyramid=human_pyramid[0]),
    MortalityByCDR(human_model, mortalityrates=human_deathrates),

    # 4. Within-species transmission (reads intervention state)
    IntervenedWaterTransmission(
        human_model, expdurdist,
        seasonality_values=season_tiled,
        containment_efficacy=0.90,
        filter_efficacy=0.95,
    ),
]

# ==== DOG MODEL COMPONENTS ====
dog_model.components = [
    SEIR.Susceptible(dog_model),
    SEIR.Exposed(dog_model, expdurdist, infdurdist),
    SEIR.Infectious(dog_model, infdurdist),
    SEIR.Recovered(dog_model),
    SEISRecovery(dog_model),

    DogTethering(dog_model,
                 village_coverage=0.50,
                 tether_efficacy=0.85),
    ABATELarvicide(dog_model,
                   coverage_by_patch=abate_coverage,  # same water sources
                   efficacy=0.80,
                   season_start_day=0,
                   season_end_day=150),

    BirthsByCBR(dog_model, birthrates=dog_birthrates,
                pyramid=dog_pyramid[0]),
    MortalityByCDR(dog_model, mortalityrates=dog_deathrates),

    IntervenedWaterTransmission(
        dog_model, expdurdist,
        seasonality_values=season_tiled,
        tether_efficacy=0.85,
    ),
]

# ==== CROSS-SPECIES COUPLING WITH INTERVENTIONS ====
cross_host = IntervenedDualHostTransmission(
    human_model, dog_model, expdurdist,
    dog_to_human=DOG_TO_HUMAN_COUPLING,
    human_to_dog=HUMAN_TO_DOG_COUPLING,
    seasonality_values=season_tiled,
    containment_efficacy=0.90,
    tether_efficacy=0.85,
    filter_efficacy=0.95,
)

# ==== IMPORTATION (unchanged) ====
human_importation = PatchImportation(...)
dog_importation = PatchImportation(...)

# ==== TICK LOOP (unchanged structure) ====
for tick in range(nticks):
    human_model._initialize_flows(tick)
    dog_model._initialize_flows(tick)

    for component in human_model.components:
        component.step(tick)
    for component in dog_model.components:
        component.step(tick)

    cross_host.step(tick)          # cross-species with interventions
    human_importation.step(tick)
    dog_importation.step(tick)
```

### Key design decisions

**Component ordering within each model:**
1. SEIS state machine (Susceptible → Exposed → Infectious → Recovered → SEISRecovery)
2. Intervention state updates (CaseContainment/DogTethering, ABATE, Filters) — must run BEFORE transmission so the FOI calculation sees current intervention state
3. Vital dynamics (births, deaths)
4. Transmission (reads intervention state to compute modified FOI)

**SEIS compatibility:** `CaseContainment` and `DogTethering` use a tri-state flag (0=pending, 1=contained/tethered, 2=escaped). The flag resets to 0 when an agent leaves the INFECTIOUS state, so after the SEIS cycle (I→R→S) a reinfected agent gets a fresh containment decision.

**Both transmission pathways are covered:** `IntervenedWaterTransmission` handles within-species gravity-coupled FOI, and `IntervenedDualHostTransmission` handles cross-species FOI. Both apply all four intervention modifiers using `hasattr` checks, so the same classes work transparently for human and dog models.
