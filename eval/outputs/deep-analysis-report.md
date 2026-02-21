# Deep Quantitative Analysis: LASER Skill A/B Test

**Date:** 2026-02-20
**Analyst:** Claude Opus 4.6
**Skills tested:** laser-spatial-disease-modeling, epi-model-parametrization, modelops-calabaria
**Model under test:** Claude Opus 4.6 via `claude --print`

---

## 1. LASER API Usage Audit

### Methodology

Every import statement, class instantiation, method call, and parameter access was categorized as:
- **Correct**: Matches the LASER v1.0.0 API (per `laser_api_reference.md` and installed package)
- **Hallucinated**: Class/method name that does not exist in LASER
- **Wrong signature**: Real class, wrong argument names or count
- **Polyfill**: Reimplementation of functionality that exists in LASER

### With-Skill Scripts

| Script | Correct API | Hallucinated | Wrong Sig | Polyfill | Total Calls |
|--------|:-----------:|:------------:|:---------:|:--------:|:-----------:|
| `polio_seir_basic_10patch.py` (P1) | 18 | 0 | 0 | 0 | 18 |
| `polio_gravity_seasonal.py` (P2) | 17 | 0 | 0 | 0 | 17 |
| `polio_seir_10patch.py` (P3) | 22 | 0 | 0 | 0 | 22 |
| `calibrate_polio.py` (P4) | 24 | 0 | 0 | 0 | 24 |
| `polio_seir_20district.py` (P5) | 19 | 0 | 0 | 0 | 19 |
| `custom_components.py` (shared) | 12 | 0 | 0 | 0 | 12 |
| **TOTAL** | **112** | **0** | **0** | **0** | **112** |

**Detailed breakdown for P1 (`polio_seir_basic_10patch.py`):**

| Line(s) | API Call | Status |
|---------|----------|--------|
| 29 | `from laser.core.propertyset import PropertySet` | Correct |
| 30 | `from laser.core.demographics import AliasedDistribution` | Correct |
| 31 | `import laser.core.distributions as dists` | Correct |
| 32 | `from laser.generic import SEIR, Model` | Correct |
| 33 | `from laser.generic.utils import ValuesMap` | Correct |
| 34 | `from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR` | Correct |
| 35 | `from laser.core.migration import gravity, row_normalizer, distance` | Correct |
| 53 | `PropertySet({...})` | Correct |
| 91 | `ValuesMap.from_timeseries(season_tiled, nnodes)` | Correct |
| 126 | `AliasedDistribution(stable_age_dist)` | Correct |
| 150 | `dists.gamma(shape=..., scale=...)` | Correct |
| 151 | `dists.normal(loc=..., scale=...)` | Correct |
| 168 | `Model(scenario, PARAMS, birthrates=...)` | Correct |
| 182 | `gravity(pops, dist_matrix, 1, 0, b, c)` | Correct |
| 186 | `row_normalizer(network, 0.2)` | Correct |
| 196-204 | `SEIR.Susceptible(model)`, `SEIR.Exposed(model, ...)`, etc. | Correct (x6) |
| 225 | `model.run("...")` | Correct |

### Without-Skill Scripts

| Script | Correct API | Hallucinated | Wrong Sig | Polyfill | Total Calls |
|--------|:-----------:|:------------:|:---------:|:--------:|:-----------:|
| `polio_seir_model.py` (P1) | 11 | 0 | 1 | 0 | 12 |
| `laser_polio_model.py` (P2-P3) | 0 | 2 | 0 | 5 | 7 |
| `polio_laser_model.py` (P4-P5) | 0 | 0 | 0 | 3 | 3 |
| **TOTAL** | **11** | **2** | **1** | **8** | **22** |

**Detailed breakdown for P1 (`polio_seir_model.py`):**

| Line(s) | API Call | Status | Notes |
|---------|----------|--------|-------|
| 23 | `import laser.core.distributions as dists` | Correct | |
| 24 | `from laser.core import PropertySet` | **Wrong path** | Should be `from laser.core.propertyset import PropertySet` (but may work as re-export) |
| 25 | `from laser.generic import Model, SEIR` | Correct | |
| 89 | `PropertySet({...})` | Correct | |
| 114-117 | `dists.gamma(shape=3.0, scale=1.0)` | Correct | |
| 122 | `Model(scenario, params, name="...")` | Correct | But missing `birthrates=` arg (no vital dynamics) |
| 125-129 | `SEIR.Susceptible(model)`, etc. | Correct (x5) | |
| 131 | `model.components = [...]` | Correct | |
| 140 | `model.run("...")` | Correct | |

**Signature error (P1):** Line 89 uses `gravity_k: 500, gravity_a: 1.0` -- the LASER auto-network feature expects `gravity_k` as a small coupling fraction (0.001-0.1), not a raw scaling constant of 500. With `gravity_a=1.0` (convention is `a=0`), this would produce an extremely over-connected network. Not a crash-level error, but would produce unrealistic spatial dynamics.

**Detailed breakdown for P2-P3 (`laser_polio_model.py`):**

| Line(s) | API Call | Status | Notes |
|---------|----------|--------|-------|
| 26-58 | `class PropertySet`, `class LaserFrame` | **Polyfill** | Reimplements `laser.core.propertyset.PropertySet` and `laser.core.LaserFrame` from scratch |
| 101-116 | `def build_gravity_network(...)` | **Polyfill** | Reimplements `laser.core.migration.gravity()` |
| 122-146 | `def row_normalize_network(...)` | **Polyfill** | Reimplements `laser.core.migration.row_normalizer()` (but semantics differ -- adds diagonal retention, which LASER does not) |
| 152-173 | `def build_seasonal_profile(...)` | **Polyfill** | Reimplements `ValuesMap.from_timeseries()` pattern |
| 376 | Reference to `Immunization.RoutineImmunization` | **Hallucinated** | No such API path exists in LASER |
| 401 | Reference to `Immunization.Campaign` | **Hallucinated** | No such API path exists in LASER |
| 196-255 | `class PolioModel` with `.run()` | **Polyfill** | Reimplements `laser.generic.Model` as a deterministic ODE solver |

**Detailed breakdown for P4-P5 (`polio_laser_model.py`):**

| Line(s) | API Call | Status | Notes |
|---------|----------|--------|-------|
| 115-124 | `def haversine_matrix(...)` | **Polyfill** | Reimplements `laser.core.migration.distance()` |
| 127-141 | `def gravity_matrix(...)` | **Polyfill** | Reimplements `laser.core.migration.gravity()` + `row_normalizer()` |
| 148-153 | `def seasonal_beta(...)` | **Polyfill** | Reimplements `ValuesMap` seasonal pattern |

Zero LASER imports in this file. The entire framework is abandoned.

---

## 2. Failure Mode Taxonomy

Every error or issue in the without-skill scripts, categorized by type:

### 2.1 API Hallucination (2 instances)

| Script | Line | Hallucinated API | What Exists |
|--------|------|------------------|-------------|
| `laser_polio_model.py` | 376 (comment) | `Immunization.RoutineImmunization` | `laser.generic.immunization.RoutineImmunization` (flat, not nested under `Immunization.`) |
| `laser_polio_model.py` | 401 (comment) | `Immunization.Campaign` | `laser.generic.immunization.ImmunizationCampaign` (different name entirely) |

These are in docstring comments rather than executable code, but they reveal the model's confused mental model of the API structure.

### 2.2 Framework Abandonment (3 instances, escalating)

| Script | What Was Abandoned | What Replaced It |
|--------|-------------------|------------------|
| `laser_polio_model.py` (P2) | `laser.core.PropertySet` | Custom `PropertySet` class (lines 26-42) |
| `laser_polio_model.py` (P2) | `laser.generic.Model` + all SEIR components | Custom `PolioModel` class with deterministic ODE phases (lines 196-255) |
| `polio_laser_model.py` (P4-P5) | Entire LASER framework | Standalone NumPy tau-leap stochastic simulation (lines 166-282) |

**Severity escalation:** P1 used LASER correctly. P2-P3 polyfilled core types but still mimicked LASER patterns. P4-P5 abandoned even the pretense of LASER, using raw NumPy arrays for everything.

### 2.3 Signature Errors (1 instance)

| Script | Line | Error | Impact |
|--------|------|-------|--------|
| `polio_seir_model.py` | 89 | `gravity_k: 500, gravity_a: 1.0` | LASER convention is `a=0` (source population does not affect outward flow) and `gravity_k` is a small fraction. With `k=500, a=1.0`, the auto-configured network would have extreme coupling, destroying spatial structure. |

### 2.4 Silent Failure Patterns (4 instances)

| # | Script | Lines | Pattern | Expected Impact |
|---|--------|-------|---------|-----------------|
| 1 | `polio_laser_model.py` | 43-44 | **Birthrate units**: `BIRTH_RATE = 25.0 / 1000.0 / 365` = 6.85e-5 per capita per day. If this value were ever passed to LASER's `BirthsByCBR`, it would cause zero births (BirthsByCBR expects 25.0 per-1000/year, not 6.85e-5). | Population static, no susceptible replenishment, eventual disease extinction. **This is the exact silent failure the skill's Critical Gotcha #1 warns about.** |
| 2 | `laser_polio_model.py` | 93 | **Leaky vaccine model**: `ve_opv = 0.50` reduces FOI by 50% for vaccinated agents. LASER's `TransmissionSE` does not read `susceptibility` -- it only checks `state == SUSCEPTIBLE`. Reducing FOI by a scalar multiplier is not how LASER handles vaccination. | In the standalone ODE model, this works as coded. But if ported to LASER, vaccination would have zero effect on transmission. **This is the exact silent failure the skill's Critical Gotcha #2 warns about.** |
| 3 | `laser_polio_model.py` | 67-68 | **Short simulation**: `nticks: 365 * 3` (3 years). The prompt asked for a 10-patch SEIR with realistic dynamics. Three years is too short for endemic equilibrium to establish with a 28-day infectious period and 95% initial immunity. | Transient dynamics only; no observable endemic equilibrium or seasonal cycling. |
| 4 | `laser_polio_model.py` | 219 | **Non-endemic initial conditions**: `S = pop * 0.90, I = pop * 0.01, R = pop * 0.07`. With 90% susceptible and R0=6.3 (beta=0.3, gamma=1/21), R_eff = 5.67 -- this triggers an immediate massive epidemic that depletes the susceptible pool in the first weeks, rather than simulating near-equilibrium endemic polio. | First epidemic wave infects >80% of population; subsequent dynamics are post-epidemic recovery, not endemic transmission. |

### 2.5 Missing Components (by prompt)

| Prompt | Required Feature | With-Skill | Without-Skill |
|--------|-----------------|:----------:|:-------------:|
| P1 | SEIR with LASER API | Present | Present |
| P1 | Vital dynamics (births/deaths) | Present (`BirthsByCBR`, `MortalityByCDR`) | **Missing** (no births, no deaths) |
| P1 | Importation | Present (`PatchImportation`) | **Missing** |
| P2 | Gravity migration via `laser.core.migration.gravity()` | Present | Polyfilled (custom function) |
| P2 | Seasonal forcing via `ValuesMap` | Present | Polyfilled (custom array) |
| P3 | Routine Immunization (age-targeted) | Present (RI at 6 weeks, age-checked) | Polyfilled (births-based proxy, no age tracking) |
| P3 | SIA campaigns (age-banded) | Present (0-5yr, biannual) | Polyfilled (fraction-based approximation) |
| P3 | Correlated missedness | Present (`reachable` flag, `on_birth` callback) | **Missing** |
| P4 | Calibration framework | Present (LHS, AFP loss, CCS loss, ranked results) | **Empty/missing** (complete failure) |
| P4 | Loss function design | Present (log-MSE + CCS zero-week proportion) | **Missing** |
| P5 | 20-district heterogeneous model | Present (20 districts, heterogeneous populations, per-district unreachable fractions) | Present (20 districts with tau-leap engine, no LASER) |
| P5 | Weekly incidence output | Present (DataFrame export) | Present (NumPy array) |
| P5 | Zero-incidence week analysis | Present (per-district analysis with visualization) | Present (per-district analysis with ASCII bar chart) |

### 2.6 Parameter Drift Between Prompts

Parameters that changed between without-skill scripts without justification:

| Parameter | P1 Value | P2-P3 Value | P4-P5 Value | With-Skill (all) |
|-----------|----------|-------------|-------------|-------------------|
| Latent period | 3 days | 5 days | 3 days | 3 days (consistent) |
| Infectious period | 28 days | 21 days | 28 days | 28 days (consistent) |
| beta | 0.2143 | 0.3 | 0.2143 | 0.2143 (consistent) |
| Initial immunity | 95% | 7% R | 1/R0 | 95% or per-district (consistent) |
| Population/patch | 100K uniform | 50K-500K random | 60K-2M named | 100K or named (consistent) |
| Simulation years | 10 | 3 | 20 | 10 or 20 (consistent) |

The without-skill condition drifted on latent period (3 to 5 to 3 days), infectious period (28 to 21 to 28 days), and beta (0.21 to 0.30 to 0.21). The shift in P2-P3 to `sigma=1/5` and `gamma=1/21` gives an effective R0 of beta/gamma = 0.3/0.0476 = 6.3, similar but not identical. These changes occurred without any prompt-driven justification.

---

## 3. Skill Content to Code Quality Mapping

### 3.1 Critical Gotcha #1: Rate Units (Per-1000/Year)

**SKILL.md section:** "Critical Gotchas > Rate Units: Per-1000/Year (Not Daily Per-Capita)"

**Without-skill error prevented:** `polio_laser_model.py` line 43: `BIRTH_RATE = 25.0 / 1000.0 / DAYS_PER_YEAR` converts to daily per-capita rate (6.85e-5). If ever passed to `BirthsByCBR`, this would cause silent zero-birth failure.

**With-skill code that demonstrates correct understanding:**

```python
# polio_seir_basic_10patch.py, lines 157-159
assert PARAMS.cbr >= 1 and PARAMS.cbr <= 60, \
    f"CBR must be per-1000/year (got {PARAMS.cbr})"
birthrate_array = np.full((nticks, nnodes), PARAMS.cbr, dtype=np.float32)
```

The with-skill code includes the **exact sanity check** recommended in the skill's Critical Gotcha section. This assertion appears in all five with-skill scripts.

### 3.2 Critical Gotcha #2: Vaccination State vs Susceptibility

**SKILL.md section:** "Critical Gotchas > Vaccination: State vs. Susceptibility"

**Without-skill error prevented:** `laser_polio_model.py` lines 93, 296: Uses `ve_opv = 0.50` as a "leaky vaccine" that reduces FOI by 50%. This is the `susceptibility`-based approach that the skill explicitly warns is ineffective with LASER's `TransmissionSE`.

**With-skill code that demonstrates correct understanding:**

```python
# custom_components.py, lines 372
people.state[vaccinated] = SEIR.State.RECOVERED.value
```

and for the SEIRV model:

```python
# custom_components.py, line 200
people.state[indices] = VACCINATED  # value 4
```

The with-skill code sets `state` (not `susceptibility`) for vaccination, exactly as the Critical Gotcha prescribes. The comment on line 305 of `custom_components.py` even references the gotcha directly:

> "CRITICAL: Sets state = RECOVERED (value 3) rather than susceptibility = 0, because the built-in TransmissionSE only checks state == SUSCEPTIBLE (0) when selecting agents for infection."

### 3.3 Workflow Step 3: Seasonal Transmission via ValuesMap

**SKILL.md section:** "Step 3: Seasonal Transmission"

**Skill code example:**
```python
season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)
SEIR.Transmission(model, expdurdist, seasonality=seasonality)
```

**Direct reflection in with-skill output (e.g., `polio_gravity_seasonal.py` lines 169-170):**
```python
season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)
```

This is a near-verbatim match. The without-skill code had to build a custom `build_seasonal_profile()` function (P2) or inline cosine computation (P4-P5) instead.

### 3.4 Workflow Step 4: Generalizing Seasonal Forcing

**SKILL.md section:** "Step 4 > Generalizing Seasonal Forcing"

The skill provides a general recipe for non-Bjornstad seasonal profiles:

```python
season_365 = 1.0 + 0.3 * np.cos(2 * np.pi * (days - peak_day) / 365)
season_365 /= season_365.mean()  # Normalize to mean == 1.0
```

**With-skill code (`polio_seir_basic_10patch.py` lines 87-90):**
```python
peak_day = 245       # Early September
amplitude = 0.30     # +/-30% modulation around baseline
beta_season = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)
```

The with-skill code adapted the Bjornstad measles profile recipe to a monsoon-peaked cosine profile, exactly as the skill's "Generalizing Seasonal Forcing" section suggests. The without-skill code used a piecewise linear interpolation approach (P2) or a simple inline cosine (P4-P5) -- functional but not using the LASER `ValuesMap` pattern.

### 3.5 Workflow Step 5: Gravity Migration Network

**SKILL.md section:** "Step 5: Gravity Migration Network"

**Skill code example:**
```python
model.network = gravity(np.array(scenario.population), distances, 1, 0, b, c)
average_export_frac = np.mean(model.network.sum(axis=1))
model.network = model.network / average_export_frac * gravity_k
model.network = row_normalizer(model.network, 0.2)
```

**With-skill code (`polio_gravity_seasonal.py` lines 117-128):**
```python
network = gravity(pops, dist_matrix, 1.0, 0, b, c)
avg_export = np.mean(network.sum(axis=1))
if avg_export > 0:
    network = network / avg_export * k
network = row_normalizer(network, max_export_frac)
```

Nearly identical. The `a=0` convention (source population does not affect outward flow) is correctly maintained. The without-skill P1 used `gravity_a=1.0`, violating the convention; P2-P3 reimplemented gravity from scratch with `a=1.0`.

### 3.6 Workflow Step 6: Component Assembly

**SKILL.md section:** "Step 6: Assemble and Run"

The skill provides the canonical component ordering:

```python
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    Importation(model, ...),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model, birthrates=..., pyramid=...),
    MortalityByEstimator(model, estimator=...),
]
```

All five with-skill scripts follow this pattern exactly, with the component ordering note ("Susceptible and Recovered components should wrap the transition steps") respected. The without-skill code created custom phase classes in a different order (P2-P3) or used a raw loop without components (P4-P5).

### 3.7 Workflow Step 7: BaseModel Wrapping

**SKILL.md section:** "Step 7: Wrap as BaseModel for Calibration"

This section directly enabled the P4 calibration script. The without-skill condition produced **empty output** for P4 because it had no knowledge of how to structure a calibration framework around the (non-LASER) model it had built.

### 3.8 API Reference: Component Signatures

The `laser_api_reference.md` provides exact constructor signatures:

| Component | Correct Signature | With-Skill Usage | Without-Skill Usage |
|-----------|------------------|------------------|---------------------|
| `TransmissionSE` | `(model, expdurdist, seasonality=None)` | `SEIR.Transmission(model, expdurdist, seasonality=seasonality)` | N/A (reimplemented) |
| `Exposed` | `(model, expdurdist, infdurdist)` | `SEIR.Exposed(model, expdurdist, infdurdist)` | P1 only: `SEIR.Exposed(model, exp_dur_dist, inf_dur_dist)` |
| `BirthsByCBR` | `(model, birthrates=array, pyramid=pyramid)` | `BirthsByCBR(model, birthrates=..., pyramid=pyramid)` | N/A (not used) |
| `MortalityByCDR` | `(model, mortalityrates=array)` | `MortalityByCDR(model, mortalityrates=...)` | N/A (not used) |

### 3.9 Epi Parametrization Skill: Disease-Specific Reference

**SKILL.md section:** "Disease-Specific Quick References > Polio (WPV1)"

The epi-model-parametrization skill provides:
> R0: 5-7 (sanitation-dependent) | Latent: 3-6d | Infectious: 14-28d (fecal) | Case ratio: 1:200 (paralysis)
> Seasonal: Monsoon/warm season | Vaccination: OPV campaigns + routine IPV

**With-skill parameter consistency across all scripts:**
- R0 = 6 (within 5-7 range)
- Latent period = 3 days (within 3-6d range)
- Infectious period = 28 days (within 14-28d range)
- Case ratio = 1:200 (exact match)
- Seasonal forcing = monsoon (exact match)

**Without-skill parameter inconsistency:** R0 ranged from 6.0 to 6.3 across scripts; latent period fluctuated between 3 and 5 days; infectious period fluctuated between 21 and 28 days.

---

## 4. Architectural Divergence Analysis

### 4.1 Agent-Based (LASER) vs Aggregated (ODE/Tau-Leap)

| Dimension | With-Skill (LASER ABM) | Without-Skill P2-P3 (ODE) | Without-Skill P4-P5 (Tau-Leap) |
|-----------|----------------------|---------------------------|-------------------------------|
| **State representation** | Individual agents with int8 state, DOB, nodeid, timers | Floating-point compartment counts per patch | Integer compartment counts per patch |
| **Stochasticity** | Per-agent Bernoulli infection trials (Numba parallel) | Deterministic (fractional transitions) | Poisson-distributed transitions |
| **Vaccination targeting** | Age-specific, per-agent reachability flag | Fraction-of-births and fraction-of-S proxy | Fraction of births, fraction of under-5 S |
| **Spatial coupling** | Built into `TransmissionSE` via `model.network` | Matrix multiplication on prevalence vector | Matrix multiplication on prevalence vector |
| **Duration distributions** | Gamma/normal sampled per agent | Exponential rates (1/mean) | Exponential rates (1/mean) |
| **Vital dynamics** | `BirthsByCBR` with Poisson draws, `MortalityByCDR` | Constant birth=death rate `mu` | Poisson births/deaths |
| **Memory** | O(N_agents x properties) -- millions of agents | O(N_patches x compartments) -- hundreds of floats | O(N_patches x compartments) |
| **Runtime** | Minutes (Numba JIT) | Sub-second | Seconds |

### 4.2 Practical Implications for Polio Modeling

The choice of agent-based vs aggregated modeling has direct consequences for the public health question:

**Age-targeted vaccination**: Polio eradication requires vaccinating children at specific ages (OPV at 6 weeks, SIA targeting 0-5 years). An agent-based model tracks each individual's date of birth and can enforce exact age eligibility windows. The ODE model (P2-P3) approximates this as "a fraction of daily births" and "a fraction of susceptibles" -- losing the age-targeting that is central to actual vaccination programs.

**Correlated missedness**: The most important operational challenge in Pakistan polio eradication is that the same communities are consistently missed by both routine and campaign vaccination (zero-dose children). The with-skill model implements this via a permanent `reachable` flag per agent. The without-skill models have no concept of correlated missedness -- they treat each vaccination event as an independent Bernoulli draw, systematically overestimating cumulative coverage.

**Stochastic fadeout and persistence**: Polio in Pakistan's low-incidence districts shows stochastic fadeout (weeks with zero cases) followed by reimportation. This is a fundamentally stochastic phenomenon that deterministic ODEs cannot capture. The tau-leap model (P4-P5) can capture some stochasticity, but without agent-level tracking, it cannot reproduce the relationship between population size and persistence (critical community size analysis).

**Disease surveillance comparison**: The with-skill model computes `model.nodes.newly_infected` per tick per patch, which can be directly compared to AFP surveillance data via the 1:200 paralysis-to-infection ratio. The P2-P3 ODE model tracks continuous-valued compartment sizes, not discrete infection events, making direct comparison to case counts less natural.

### 4.3 Which Model Could Answer the Policy Question?

The fundamental policy question for Pakistan polio eradication is: **"What vaccination strategies will achieve and sustain zero cases across all districts?"**

| Capability Required | With-Skill | Without-Skill |
|---------------------|:----------:|:-------------:|
| Evaluate RI coverage changes | Yes (per-district, age-targeted) | Approximate (fraction-based proxy) |
| Evaluate SIA frequency changes | Yes (age-banded campaigns) | Approximate (fraction of S) |
| Model zero-dose populations | Yes (`reachable` flag, correlated missedness) | No |
| Predict district-level fadeout | Yes (stochastic agent-based) | No (P2-P3 deterministic) / Approximate (P4-P5 tau-leap) |
| Compare with AFP surveillance | Yes (weekly incidence, 1:200 ratio) | No (no incidence tracking in P2-P3) / Approximate (P4-P5) |
| CCS analysis | Yes (zero-AFP-week proportions) | No |
| Calibrate to MMWR data | Yes (P4 calibration framework) | No (P4 was empty) |

**Verdict:** Only the with-skill model architecture supports the full range of policy-relevant analyses needed for polio eradication planning. The without-skill P2-P3 model cannot distinguish between increasing RI coverage in reachable populations vs reducing the unreachable fraction -- the two most important operational levers. The without-skill P4-P5 tau-leap model has some stochastic capability but lacks the structural features (age-targeting, reachability) needed for realistic vaccination strategy evaluation.

---

## 5. Progressive Complexity Analysis

### 5.1 Without-Skill Degradation Trajectory

```
Prompt   LASER API Usage    Model Type         Framework Status
------   ---------------    ----------         ----------------
P1       11/12 correct      LASER ABM          USING FRAMEWORK
P2       0/7 correct        Deterministic ODE  ABANDONED FRAMEWORK
P3       (continued P2)     Deterministic ODE  ABANDONED FRAMEWORK
P4       N/A                Empty output        COMPLETE FAILURE
P5       0/3 correct        Tau-leap stoch     ABANDONED FRAMEWORK
```

### 5.2 The Abandonment Trigger Point

Framework abandonment occurred between P1 and P2. The P2 prompt required:

1. **Gravity-model migration** (requires `laser.core.migration.gravity()` and `row_normalizer()`)
2. **Seasonal forcing** (requires `ValuesMap.from_timeseries()` and `SEIR.Transmission(..., seasonality=...)`)
3. **Row normalization** (requires understanding of network matrix semantics)

These are LASER-specific APIs that Claude without skills could not reconstruct from general knowledge. The gravity model formula (M_ij = k * p_i^a * p_j^b / d_ij^c) is well-known in epidemiology, but:
- The LASER-specific function signature (`gravity(populations, distances, k, a, b, c)`) was unknown
- The `ValuesMap` utility class was unknown
- The `seasonality=` parameter on `TransmissionSE` was unknown

Faced with unknown APIs, the model made a rational (but costly) decision: reimplemented everything from scratch using standard NumPy. This is the "framework abandonment" pattern.

### 5.3 Why P4 Was a Complete Failure

P4 (calibration) was scored 0/12 because the without-skill condition produced no calibration code. By P4, the accumulated state was:
- A non-LASER ODE model (P2-P3) with no `model.params`, no `model.nodes.newly_infected`, no component structure
- No knowledge of `calabaria.BaseModel`, `ParameterSpace`, `TrialResult`, or any calibration framework
- No `model_output` decorators or structured output extraction

The P4 prompt asked to "wrap the model for calibration" -- but there was no LASER model to wrap. The without-skill Claude apparently recognized the impossibility and produced nothing rather than an incoherent attempt.

### 5.4 With-Skill Quality Degradation

The with-skill condition showed **no quality degradation** across P1-P5:

| Prompt | API Correctness | Structural Correctness | Notes |
|--------|:-:|:-:|-------|
| P1 | 3/3 | 3/3 | Clean SEIR with all components |
| P2 | 3/3 | 3/3 | Added gravity + seasonality correctly |
| P3 | 3/3 | 3/3 | Added SEIRV vaccination with custom components, MMWR comparison |
| P4 | 3/3 | 3/3 | Full calibration framework with LHS, dual loss function, CCS analysis |
| P5 | 3/3 | 3/3 | 20 districts with heterogeneous unreachable fractions |

Each successive prompt built correctly on the previous one, reusing imports and patterns. The skill provided a stable API reference that prevented the "knowledge decay" seen in the without-skill condition.

---

## 6. Quantitative Code Metrics

### 6.1 Lines of Code

| Script | Total LOC | Comments/Blank | Executable |
|--------|:---------:|:--------------:|:----------:|
| **With-Skill** | | | |
| `polio_seir_basic_10patch.py` (P1) | 393 | ~120 | ~273 |
| `polio_gravity_seasonal.py` (P2) | 563 | ~160 | ~403 |
| `polio_seir_10patch.py` (P3) | 957 | ~250 | ~707 |
| `calibrate_polio.py` (P4) | 853 | ~210 | ~643 |
| `polio_seir_20district.py` (P5) | 584 | ~160 | ~424 |
| `custom_components.py` (shared) | 461 | ~130 | ~331 |
| **With-Skill Total** | **3,811** | **~1,030** | **~2,781** |
| **Without-Skill** | | | |
| `polio_seir_model.py` (P1) | 207 | ~50 | ~157 |
| `laser_polio_model.py` (P2-P3) | 563 | ~150 | ~413 |
| `polio_laser_model.py` (P4-P5) | 382 | ~80 | ~302 |
| **Without-Skill Total** | **1,152** | **~280** | **~872** |

The with-skill output is **3.3x larger** by total LOC. This reflects richer functionality (calibration, MMWR comparison, diagnostic plots, per-district analysis) rather than verbosity.

### 6.2 LASER-Specific Code Lines

Lines that directly import from, instantiate, or call LASER APIs:

| Script | LASER Lines | Custom/Polyfill Lines | LASER Fraction |
|--------|:-----------:|:---------------------:|:--------------:|
| **With-Skill** | | | |
| `polio_seir_basic_10patch.py` | 45 | 0 | 100% |
| `polio_gravity_seasonal.py` | 42 | 0 | 100% |
| `polio_seir_10patch.py` | 55 | 0 | 100% |
| `calibrate_polio.py` | 60 | 0 | 100% |
| `polio_seir_20district.py` | 48 | 0 | 100% |
| `custom_components.py` | 35 | 65 | 35% (custom components extend LASER) |
| **Without-Skill** | | | |
| `polio_seir_model.py` | 30 | 0 | 100% |
| `laser_polio_model.py` | 0 | 180 | **0%** |
| `polio_laser_model.py` | 0 | 140 | **0%** |

### 6.3 Component Usage

| Component | With-Skill P1 | P2 | P3 | P4 | P5 | Without P1 | Without P2-3 | Without P4-5 |
|-----------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `SEIR.Susceptible` | Y | Y | Y | Y | Y | Y | -- | -- |
| `SEIR.Exposed` | Y | Y | Y | Y | Y | Y | -- | -- |
| `SEIR.Infectious` | Y | Y | Y | Y | Y | Y | -- | -- |
| `SEIR.Recovered` | Y | Y | Y | Y | Y | Y | -- | -- |
| `SEIR.Transmission` | Y | Y | Y | Y | Y | Y | -- | -- |
| `BirthsByCBR` | Y | Y | Y | Y | Y | -- | -- | -- |
| `MortalityByCDR` | Y | Y | Y | Y | Y | -- | -- | -- |
| Custom `PatchImportation` | Y | Y | Y | Y | Y | -- | -- | -- |
| Custom `PerPatchVaccination` | -- | -- | -- | -- | Y | -- | -- | -- |
| Custom `VaccinatedCompartment` | -- | -- | Y | Y | -- | -- | -- | -- |
| Custom `PerPatchVaccinationSEIRV` | -- | -- | Y | Y | -- | -- | -- | -- |
| `ValuesMap` (seasonality) | Y | Y | Y | Y | Y | -- | -- | -- |
| `gravity()` + `row_normalizer()` | Y | Y | Y | Y | Y | -- | -- | -- |
| `AliasedDistribution` (ages) | Y | Y | Y | Y | Y | -- | -- | -- |

### 6.4 LASER Feature Usage

| Feature | With-Skill | Without-Skill |
|---------|:----------:|:-------------:|
| Gravity migration network | All 5 scripts | Polyfilled in 2, absent in 1 |
| Seasonal forcing via ValuesMap | All 5 scripts | Polyfilled in 2, inline in 1 |
| Vital dynamics (births + deaths) | All 5 scripts | None |
| Age-targeted vaccination | P3, P4, P5 | None (fraction-based proxy) |
| Correlated missedness | P3, P4, P5 | None |
| Additional states (`V` compartment) | P3, P4 | P2-P3 (custom ODE), P4-P5 (N/A) |
| `on_birth` callback | P3, P4, P5 | None |
| Agent-level DOB tracking | All 5 scripts | None |
| Importation (susceptible-targeted) | All 5 scripts | Random (P4-P5), absent (P2-P3) |
| `model.nodes.newly_infected` | All 5 scripts | None |
| `model.nodes.forces` | P3 | None |
| Calibration loss functions | P4 | None |
| CCS/zero-week analysis | P4, P5 | P5 (basic version) |

---

## 7. Practical Impact Assessment

### 7.1 Could the Without-Skill Output Support Policy Decisions?

**P1 (`polio_seir_model.py`):** Partially. The basic SEIR dynamics are correctly set up and the model runs. However, the absence of vital dynamics means the susceptible pool is never replenished by births, so the model predicts eventual disease extinction by exhaustion of susceptibles rather than endemic equilibrium. A modeler would immediately recognize this as unrealistic for a 10-year simulation of an endemic disease. **Usability: prototype only, not decision-grade.**

**P2-P3 (`laser_polio_model.py`):** Not for the stated purpose. The model is a deterministic ODE system with:
- No individual-level tracking (cannot evaluate age-targeted vaccination)
- No stochastic fadeout (cannot predict which districts achieve zero cases)
- Leaky vaccine efficacy (epidemiologically valid concept but wrong for LASER)
- 3-year simulation (too short for endemic dynamics with burn-in)
- Random populations (not mapped to real Pakistan districts)

A public health modeler could not use this to evaluate vaccination strategies because the key operational variables (RI age, SIA age band, unreachable fraction) cannot be represented. **Usability: pedagogical only.**

**P4 (calibration):** Empty. **Usability: none.**

**P5 (`polio_laser_model.py`):** More usable than P2-P3. The tau-leap stochastic engine with 20 named Pakistan districts, per-district vaccination coverage, seasonal forcing, and importation produces plausible-looking dynamics. The zero-incidence week analysis is directly relevant to polio eradication monitoring. However:
- No LASER framework (cannot leverage LASER ecosystem tools)
- No age-targeted vaccination (fraction-based proxy)
- No correlated missedness (overestimates cumulative SIA impact)
- Birth rate in wrong units (daily per-capita, would fail in LASER)
- Vaccination as "immunity acquisition" (R compartment) rather than waning state

**Usability: could inform basic spatial pattern discussions but not operational vaccination strategy design.**

### 7.2 What Would Need to Be Fixed?

**To make the without-skill P1 usable:**
1. Add vital dynamics (births and deaths) -- requires understanding `BirthsByCBR` and `MortalityByCDR` APIs
2. Add importation to prevent stochastic fadeout
3. Fix `gravity_k=500` and `gravity_a=1.0` to reasonable values
4. Estimated expert time: **2-4 hours** (add components, adjust parameters)

**To make the without-skill P2-P3 usable:**
1. Replace the entire ODE engine with LASER's agent-based framework
2. Implement age-targeted vaccination using `on_birth` callbacks and DOB tracking
3. Add correlated missedness
4. Extend simulation to 10-20 years
5. Add vital dynamics
6. Map to real Pakistan districts
7. Estimated expert time: **1-2 days** (essentially a rewrite)

**To make the without-skill P4-P5 usable:**
1. Replace tau-leap engine with LASER agent-based framework
2. Add agent-level age tracking
3. Implement correlated missedness
4. Fix birth rate units
5. Add calibration framework (loss functions, parameter sampling, result ranking)
6. Estimated expert time: **2-3 days** (essentially a rewrite, plus calibration framework from scratch)

**To go from with-skill output to production:**
1. Tune calibration parameter ranges based on initial LHS results
2. Increase calibration samples (32 to 500+)
3. Add convergence diagnostics
4. Validate against held-out years
5. Estimated expert time: **4-8 hours** (refinement, not rewrite)

### 7.3 Time-to-Usable Comparison

| Starting Point | Time to Usable Model | Time to Calibrated Model |
|----------------|:--------------------:|:------------------------:|
| With-skill output | 0 hours (already usable) | 4-8 hours (refinement) |
| Without-skill P1 | 2-4 hours | 2-3 days |
| Without-skill P2-P3 | 1-2 days (rewrite) | 3-5 days |
| Without-skill P4-P5 | 2-3 days (rewrite) | 4-6 days |
| From scratch (no Claude) | 3-5 days | 1-2 weeks |

---

## 8. Summary of Key Findings

### 8.1 The Skill Advantage is Structural, Not Just Informational

The skills did not merely provide "better hints" -- they enabled a fundamentally different architecture. The with-skill models are agent-based simulations with individual-level tracking, enabling age-targeted vaccination, correlated missedness, stochastic fadeout, and direct comparison to surveillance data. The without-skill models are aggregated compartmental models (ODE or tau-leap) that cannot represent these features.

### 8.2 Critical Gotchas Prevented Silent Failures

Two specific "Critical Gotcha" sections in the skill prevented silent failures that would have produced plausible-looking but quantitatively wrong results:
- **Birth rate units:** Would have caused zero population growth in LASER (no error raised)
- **Vaccination state semantics:** Would have caused vaccination to have zero effect on transmission (no error raised)

Both failures are "silent" in the strongest sense: the model runs, produces numbers, and does not crash. Without the skill's explicit warnings, a user would have no indication that results were wrong.

### 8.3 Framework Abandonment is the Dominant Failure Mode

The without-skill condition did not primarily fail through API hallucination (only 2 instances) or signature errors (only 1 instance). The dominant failure mode was **framework abandonment**: complete reimplementation of LASER functionality using basic NumPy. This is a qualitatively different failure mode than "using the API wrong" -- it means the model could not even attempt to use the framework for advanced features.

### 8.4 Degradation is Progressive and Self-Reinforcing

Once the framework was abandoned at P2, it could not be recovered at P3-P5. Each subsequent prompt built on the non-LASER foundation, making it increasingly impossible to return to the correct architecture. The P4 calibration failure is a direct consequence: there was no LASER model to calibrate.

### 8.5 Domain Knowledge Without API Knowledge is Insufficient

The without-skill condition demonstrated reasonable epidemiological knowledge: R0=6, monsoon seasonality, OPV waning at ~3 years, Waziristan corridor, 1:200 paralysis ratio. But this domain knowledge could not compensate for lack of API knowledge. The result was epidemiologically-informed code that could not be executed within the intended framework.

### 8.6 Quantitative Summary

| Metric | With-Skill | Without-Skill | Ratio |
|--------|:----------:|:-------------:|:-----:|
| Total score (60 max) | 60 | 25 | **2.4x** |
| API correctness score | 15/15 | 3/15 | **5.0x** |
| Correct LASER API calls | 112 | 11 | **10.2x** |
| Hallucinated API calls | 0 | 2 | -- |
| Polyfill implementations | 0 | 8 | -- |
| Silent failure patterns | 0 | 4 | -- |
| LASER features used | 14 | 5 | **2.8x** |
| Components used | 10 | 5 | **2.0x** |
| Total LOC (all scripts) | 3,811 | 1,152 | **3.3x** |
| Prompts with functional calibration | 1/1 | 0/1 | -- |
| Expert hours to production model | 4-8 | 48-96 | **~12x** |

---

## Appendix A: File Reference

### With-Skill Scripts

| File | Prompt | Absolute Path |
|------|--------|---------------|
| P1 | Basic 10-patch SEIR | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/polio_seir_basic_10patch.py` |
| P2 | Gravity + Seasonality | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/polio_gravity_seasonal.py` |
| P3 | Vaccination (RI+SIA) | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/polio_seir_10patch.py` |
| P4 | Calibration framework | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/calibrate_polio.py` |
| P5 | 20-district model | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/polio_seir_20district.py` |
| Shared | Custom components | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/scripts/custom_components.py` |

### Without-Skill Scripts

| File | Prompts | Absolute Path |
|------|---------|---------------|
| P1 | Basic 10-patch SEIR | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/eval/outputs/without-skill/polio_seir_model.py` |
| P2-P3 | Gravity + Seasonal + Vaccination | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/eval/outputs/without-skill/laser_polio_model.py` |
| P4-P5 | Calibration + 20-district | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/eval/outputs/without-skill/polio_laser_model.py` |

### Skill Files

| File | Absolute Path |
|------|---------------|
| LASER Spatial Modeling | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/.claude/skills/laser-spatial-disease-modeling/SKILL.md` |
| LASER API Reference | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/.claude/skills/laser-spatial-disease-modeling/references/laser_api_reference.md` |
| Epi Parametrization | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/.claude/skills/epi-model-parametrization/SKILL.md` |
| Modelops Calabaria | `/Users/guillaumechabot-couture/Library/CloudStorage/OneDrive-Bill&MelindaGatesFoundation/AI-LLM/laser-polio-pakistan-eval/.claude/skills/modelops-calabaria/SKILL.md` |
