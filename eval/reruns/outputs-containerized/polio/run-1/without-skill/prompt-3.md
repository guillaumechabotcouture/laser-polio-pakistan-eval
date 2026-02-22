Here is the complete component setup across three files. The `/eval-work` directory is read-only, so the files are at `/tmp/polio_model/`.

---

## Design overview

The core problem is that LASER's built-in immunization classes only set `susceptibility = 0` — they have no concept of *why* an agent is immune. For polio, OPV (mucosal, waning) and natural infection (systemic, lifelong) must be tracked separately, so three custom LaserFrame properties are added before any component is constructed:

| Property | dtype | Meaning |
|---|---|---|
| `immunity_source` | uint8 | 0 = susceptible, 1 = OPV, 2 = natural |
| `opv_doses` | uint8 | cumulative dose count (RI + SIA) |
| `opv_wane_tick` | int32 | tick OPV expires; -1 = never |

---

## `params.py` — all parameters in one PropertySet

Key values:
- `ri_age_days = 42`, `ri_check_period = 7`, `ri_coverage = 0.80`
- `sia_period_days = 182`, `sia_coverage = 0.90`, `sia_age_upper_days = 1825`
- `opv_waning_days = 1095`, `min_doses_for_no_waning = 4`

---

## `vaccination.py` — four components

**`RoutineImmunizationOPV(RoutineImmunization)`** — uses LASER's built-in class

`RoutineImmunization` fires every `period=7` ticks and vaccinates agents in the age window `[42 − 3, 42 + 4) = [39, 46)` days — exactly the one-week band a child passes through once. The subclass snapshots `susceptibility` before calling `super().__call__()`, then identifies who changed from 1→0 to stamp `immunity_source=1`, increment `opv_doses`, and set `opv_wane_tick`.

**`SIACampaignOPV(ImmunizationCampaign)`** — uses LASER's built-in class

`ImmunizationCampaign` covers the full 0–1825-day band with a 90% Binomial draw every 182 ticks. The subclass handles three groups: susceptible children get vaccinated for the first time; OPV-immune children receive a booster (dose count + extended wane deadline); naturally-immune children receive dose count only (natural immunity is not downgraded).

**`NaturalImmunityTracker`** — custom

Runs immediately after `InfectiousIR` each tick. Detects agents with `state == RECOVERED AND immunity_source == 0` (never labelled) and stamps them `immunity_source = 2`. No wane tick is set; natural immunity is lifelong.

**`OPVImmunityWaning`** — custom

Runs last each tick. Reverts OPV-immune agents (`immunity_source == 1`) to susceptible when `opv_wane_tick <= tick`, but only if their dose count is below the durable-immunity threshold (default 4). Agents who have received ≥4 doses through combined RI + SIA are exempt.

---

## Component execution order per tick

```
BirthsByCBR          — newborns added (start susceptible)
MortalityByCDR       — deaths removed
TransmissionSE       — force-of-infection S→E (skips susceptibility==0)
Exposed              — etimer countdown; E→I
InfectiousIR         — itimer countdown; I→R
NaturalImmunityTracker  ← custom; stamps immunity_source=2 on new R agents
RoutineImmunizationOPV  ← extends RoutineImmunization
SIACampaignOPV          ← extends ImmunizationCampaign
OPVImmunityWaning       ← custom; wanes expired OPV immunity
```

Disease progression completes before any vaccination fires, and waning runs last so a newly-vaccinated agent cannot lose immunity in the same tick.
