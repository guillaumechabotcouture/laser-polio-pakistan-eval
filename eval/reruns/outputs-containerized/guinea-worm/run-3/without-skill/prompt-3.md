All files are in place. Here is a complete description of everything that was built.

---

## File structure

```
guinea_worm_model/
├── __init__.py
├── parameters.py               GWParameters dataclass
├── model.py                    GuineaWormModel (LaserFrame + component registry)
├── run_simulation.py           6-scenario comparison runner
└── components/
    ├── __init__.py
    ├── seis_dynamics.py        SEISHumanDynamics, SEISDogDynamics
    ├── water_sources.py        WaterSourceDynamics (copepod reservoir)
    ├── abate_larvicide.py      AbateLarvicide
    ├── water_filter.py         WaterFilterIntervention
    ├── case_containment.py     CaseContainment
    └── dog_management.py       DogTethering
```

---

## SEIS loop (no immunity)

```
S ──(λ·dt)──► E ──(σ·dt)──► I ──(γ·dt)──► S
                                            ▲
                                   recovery returns to S
```
- **σ** = 1/365 d⁻¹ (12-month incubation)
- **γ** = 1/30 d⁻¹ (30-day worm emergence)
- Transitions are binomial draws: `p = 1 − exp(−rate·dt)` so counts stay integer and non-negative. Applies to both humans and dogs.

---

## How the four interventions plug in

### Component ordering each tick

```
[AbateLarvicide]          → writes abate_copepod_multiplier
[WaterFilterIntervention] → writes filter_foi_multiplier
[CaseContainment]         → writes contained_frac, contained_water_reduction
[DogTethering]            → writes tethered_frac, tethered_water_reduction
        ↓
[SEISHumanDynamics]       → reads filter_foi_multiplier for λ_h
[SEISDogDynamics]         → reads infected_copepod_frac for λ_d
        ↓
[WaterSourceDynamics]     → reads all four multipliers to compute
                            effective contamination → updates infected_copepod_frac
```

### Combined effect on force of infection

| Pathway | Multiplier formula | ~Impact at defaults |
|---|---|---|
| **Filters** on `λ_h` | `1 − coverage × efficacy` = `1 − 0.60×0.95` | −57 % |
| **ABATE** on copepod pool | `1 − coverage × efficacy` = `1 − 0.50×0.80` (mean) | −40 % |
| **Case containment** on human→water | `(1−c_f) + c_f×(1−c_r)` = `0.30+0.07` | −63 % |
| **Dog tethering** on dog→water | `(1−t_f) + t_f×(1−t_e)` = `0.50+0.075` | −42 % |

Multipliers from all four sources compound through the shared `infected_copepod_frac` state variable, so the full-stack scenario naturally captures interaction effects rather than additive reductions.

### Key design choices
- **`contained_frac` initialised to 0.0**, not to the parameter value — so baseline runs (no `CaseContainment` registered) correctly have zero containment effect.
- **`hasattr` guards** in `SEISHumanDynamics` and `WaterSourceDynamics` allow any intervention to be omitted without code changes; absent multipliers fall back to 1.0 (no effect).
- **ABATE is seasonal** — the multiplier resets to 1.0 outside day 60–300; the timer is also reset so the first in-season tick immediately triggers application.
- **Dog tethering is conditional** on `I_dogs > 0` per district, matching the field reality that tethering is only activated where infected dogs have been observed.
