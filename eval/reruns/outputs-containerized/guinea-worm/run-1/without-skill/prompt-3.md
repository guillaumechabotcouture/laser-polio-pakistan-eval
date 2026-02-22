All 7 files are written to `/home/evaluser/guinea_worm/` (the working directory `/eval-work` is read-only in this environment). Here is a complete walkthrough of every component.

---

## Package layout

```
guinea_worm/
├── __init__.py
├── params.py                   ← parameters + Chad district registry
├── model.py                    ← GuineaWormModel (extends laser.generic.Model)
└── components/
    ├── __init__.py
    ├── seis.py                 ← SEIS core dynamics (3 components)
    └── interventions.py        ← 4 prevention interventions
run_simulation.py               ← example entry point
```

---

## SEIS vs SEIR — the key LASER difference

Guinea worm confers no lasting immunity, so recovered individuals return directly to S. The `laser.generic` framework provides `InfectiousIR` (I→R) for SEIR, but not SEIS. These components implement the I→S path from scratch and pass `states={"S","E","I"}` (no R) to `Model.__init__`, so `_initialize_flows` only copies the three compartments forward each tick.

---

## Component execution order (matters)

| # | Component | What it does |
|---|-----------|-------------|
| 1 | `AbateLarvicide` | Sets `nodes.abate_active` before FOI is computed |
| 2 | `DogTethering` | Sets `nodes.tether_active` before FOI is computed |
| 3 | `ClothPipeFilter` | No-op step; `has_filter` assigned to agents in `__init__` |
| 4 | `SEISTransmission` | Reads all four modifiers, draws S→E exposures |
| 5 | `ExposedProgression` | Decrements `etimer`; transitions E→I, sets `itimer=infectious_days` |
| 6 | `CaseContainment` | Identifies new I cases (`itimer==infectious_days`), sets `is_contained` |
| 7 | `InfectiousRecovery` | Decrements `itimer`; transitions I→S, clears `is_contained` |

---

## Intervention mechanics

### 1 — ABATE larvicide (`AbateLarvicide`)
Sets `nodes.abate_active` to `True` throughout the dry-season window (Julian days 152–304). `SEISTransmission` then applies:
```
abate_factor[j] = 1 - abate_active[j] × abate_coverage[j] × 0.80
```
For Mandoul (80% coverage): factor = `1 - 0.80×0.80 = 0.36` → 64% FOI reduction.

### 2 — Water filters (`ClothPipeFilter`)
Assigns `people.has_filter` once at `__init__` via a Bernoulli(0.60) draw. During each transmission step, filter users have their FOI scaled:
```
filter_factor_i = (1 - 0.95) = 0.05   # 95% ingestion reduction
P_expose_i = 1 - exp(-foi_base[j] × filter_factor_i)
```

### 3 — Case containment (`CaseContainment`)
Runs immediately after `ExposedProgression`. Detects new I cases (`itimer == infectious_days`) and marks 70% of them `is_contained = True`. `SEISTransmission` then counts their contribution as:
```
effective_I_H[j] = n_uncontained[j] + n_contained[j] × (1 - 0.90)
```
Contained cases still contribute 10% contamination (brief water contact before detection).

### 4 — Dog tethering (`DogTethering`)
Assigns tethering capacity to 50% of nodes at `__init__`. Each tick, tethering activates only in assigned nodes where `dog_infectious_count > 0`:
```
tether_factor[j] = (1 - tether_active[j]) + tether_active[j] × (1 - 0.85)
effective_I_D[j] = dog_infectious_count[j] × tether_factor[j]
```

---

## Properties set up on `model.people` and `model.nodes`

| Frame | Property | Type | Purpose |
|-------|----------|------|---------|
| `people` | `state` | `int8` | S=0, E=1, I=2 |
| `people` | `node_id` | `int32` | patch assignment |
| `people` | `etimer` | `int32` | incubation countdown |
| `people` | `itimer` | `int32` | worm-emergence countdown |
| `people` | `has_filter` | `bool` | owns/uses water filter |
| `people` | `is_contained` | `bool` | water access restricted |
| `nodes` | `S/E/I` | `int32 (nticks+1, n)` | compartment time series |
| `nodes` | `abate_active` | `bool` | ABATE in effect this tick |
| `nodes` | `abate_coverage` | `float32` | fraction of sources treated |
| `nodes` | `tether_active` | `bool` | dog tethering in effect |
| `nodes` | `dog_infectious_count` | `int32` | infectious dogs per patch |
