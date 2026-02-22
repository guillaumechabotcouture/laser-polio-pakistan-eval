The model is at `/tmp/polio_model.py` (the `/eval-work` filesystem is read-only). It ran to completion. Here is the complete annotated listing:

---

## Design overview

### Why no built-in LASER immunization classes

`laser_core` v0.4.0 is a **low-level framework** — it provides three primitives (`LaserFrame`, `PropertySet`, `SortedQueue`) but no immunization-specific base classes. Every component is a custom callable that follows LASER's convention:

```python
class SomeComponent:
    def __init__(self, model): ...
    def __call__(self, model, tick): ...
```

`SortedQueue` *is* the right LASER tool for scheduling routine immunization events (it is a Numba-accelerated min-heap used extensively in LASER disease models for timed events).

---

### Disease states

```python
S     = 0   # Susceptible – no immunity
E     = 1   # Exposed / incubating
I     = 2   # Infectious (excreting poliovirus)
R_OPV = 3   # OPV immune  – mucosal, tracked separately
R_NAT = 4   # Natural immune – systemic + mucosal, stronger
```

### LaserFrame properties per agent

| Property | dtype | Purpose |
|---|---|---|
| `node` | int32 | district assignment |
| `age_days` | int32 | current age in days |
| `inf_status` | uint8 | S/E/I/R_OPV/R_NAT (0–4) |
| `inf_timer` | int16 | days until next stage transition |
| `ri_due` | int32 | tick when RI is due — **SortedQueue sorts on this** |
| `alive` | bool | flag; avoids `squash()` which would invalidate queue indices |

---

### Component pipeline

| Order | Component | Key design |
|---|---|---|
| 1 | `AgeingComponent` | masked `+= 1` over alive slice |
| 2 | `DeathsComponent` | Binomial draw; sets `alive=False`, zero-out `inf_status` |
| 3 | `BirthsComponent` | `pop.add()` + sets `ri_due` + **`ri_queue.push(idx)`** |
| 4 | `TransmissionComponent` | Per-node FOI; R_OPV ×0.15, R_NAT ×0.02 susceptibility |
| 5 | `ProgressionComponent` | E→I→**R_NAT** (natural always wins over prior OPV) |
| 6 | `RoutineImmunizationComponent` | **`SortedQueue` pop loop**: only today's RI cohort |
| 7 | `SIAComponent` | Every 180 days, ages 0–5, 90 % Bernoulli |
| 8 | `ReporterComponent` | Writes all 5 compartments to `ts` array properties |

---

### Key vaccination details

**RoutineImmunizationComponent — SortedQueue scheduling**

```python
# At birth (BirthsComponent):
pop.ri_due[start:end] = tick + params.ri_age_days   # set value first
for idx in range(start, end):
    model.ri_queue.push(idx)                         # heap indexes into ri_due

# At 6 weeks (RoutineImmunizationComponent):
while len(queue) > 0 and queue.peekv() <= tick:
    idx = int(queue.popi())
    if not pop.alive[idx]: continue                  # child died before dose
    if rng.random() >= 0.80: continue                # 20% miss the dose
    if pop.inf_status[idx] == S:
        pop.inf_status[idx] = R_OPV
```

This is O(due_today · log N) rather than scanning all 300 k agents every day.

**SIAComponent — periodic campaign**

```python
# Fires every 180 days on children 0-5 years
eligible   = alive & (age <= 5 * 365)
vaccinated = eligible & (rng.random(n) < 0.90)
status[vaccinated & (status == S)] = R_OPV
# R_NAT is never downgraded to R_OPV
```

**OPV vs natural immunity asymmetry**

```python
# TransmissionComponent, per node:
new_e_opv = R_OPV_mask & (draw < foi * 0.15)   # 15% rel. susceptibility
new_e_nat = R_NAT_mask & (draw < foi * 0.02)   # 2%  rel. susceptibility

# ProgressionComponent:
# I -> R_NAT regardless of prior vaccination history
#      (poliovirus replication induces stronger systemic immunity than OPV)
```

**Sample output from the 730-day run:**

```
=== Day 730 state ===
  Total alive        :   312,855
  S  - Susceptible   :    78,915
  R_OPV - OPV immune :    25,618   ← tracked separately
  R_NAT - Nat. immune:   208,322   ← tracked separately

=== Cumulative vaccination ===
  RI doses  (6-wk OPV):    12,902
  SIA doses (campaigns):   138,100
  SIA campaigns run    :         4
```
