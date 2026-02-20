---
name: laser-spatial-disease-modeling
description: This skill should be used for building spatial disease transmission
  models using the LASER (Light Agent Spatial modeling for ERadication) framework
  (v1.0.0+). This skill is appropriate when the user asks to model disease spread
  across geographic populations, simulate SEIR-type diseases with spatial coupling,
  set up gravity-model migration networks, wrap LASER models as BaseModel for
  calabaria calibration, or reproduce spatial epidemiological phenomena like
  traveling waves and critical community size. Trigger phrases include "LASER model",
  "spatial disease model", "gravity model transmission", "SEIR spatial simulation",
  "wavelet phase analysis", "critical community size", "build SEIR model",
  "wrap model as BaseModel", or "LASER BaseModel".
---

# Spatial Disease Transmission Modeling with LASER

## Overview

The LASER framework (v1.0.0+, December 2025) enables spatially-explicit, agent-based disease transmission modeling across geographic populations. This skill provides a complete workflow for building, calibrating, and analyzing spatial SEIR models, based on the canonical England and Wales pre-vaccine measles dataset. The workflow demonstrates reproduction of key spatiotemporal phenomena including traveling waves of infection and the relationship between disease extinction/reimportation and population size (critical community size).

**Packages:** `pip install laser-generic` (installs both `laser-core` and `laser-generic`)

---

## Critical Gotchas

### Rate Units: Per-1000/Year (Not Daily Per-Capita)

`BirthsByCBR`, `MortalityByCDR`, and `ConstantPopVitalDynamics` expect rates in **per-1000 population per year**. They divide by 1000 internally. If you pass pre-converted daily per-capita rates (e.g., 0.00008 instead of 30):

- `calc_capacity()` sees near-zero growth → capacity ≈ initial population
- `LaserFrame.add()` has no free slots → **no births occur**
- The model still runs and produces plausible-looking results from the initial susceptible pool
- **This is a silent failure** — no error is raised

**Sanity check** (add before model construction):
```python
assert np.all(birthrates >= 1) and np.all(birthrates <= 60), \
    f"Birthrates must be per-1000/year (typical range 10-50), got {birthrates.min():.4f}-{birthrates.max():.4f}"
```

### Vaccination: State vs. Susceptibility

`ImmunizationCampaign` and `RoutineImmunization` set `susceptibility = 0`, but **all Transmission kernels** (`TransmissionSE`, `TransmissionSI`, `TransmissionSIx`) only check `state == SUSCEPTIBLE` (int8 == 0). They do **not** check `susceptibility`.

**Result:** Built-in vaccination via these classes has **zero effect on transmission**.

**Correct approach:** Vaccination must set `state = State.RECOVERED.value` (3) and update `nodes.S` / `nodes.R` counts. Use either:
- `RoutineImmunizationEx` (built-in, sets state correctly)
- The `VaccinationCampaign` class in `scripts/custom_components.py` (supports correlated missedness)

> **Do NOT use** `ImmunizationCampaign` or `RoutineImmunization` if you need vaccination to actually reduce transmission.

---

## Workflow

### Step 1: Environment Setup and Data Loading

```python
import numpy as np
import pandas as pd
import numba as nb
from pathlib import Path
from laser.core.propertyset import PropertySet
import laser.core.distributions as dists
from laser.core.demographics import AliasedDistribution, KaplanMeierEstimator
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByEstimator
from laser.core.migration import gravity, row_normalizer
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
import pywt
```

Additional imports available in v1.0.0:
```python
# Alternative migration models
from laser.core.migration import competing_destinations, stouffer, radiation, distance
# Built-in importation (see Step 3 note)
from laser.generic.importation import Infect_Random_Agents, Infect_Agents_In_Patch
# Direct component imports (alternative to SEIR.* aliases)
from laser.generic.components import Susceptible, TransmissionSE, Exposed, InfectiousIR, Recovered
```

Data requirements:
- **Population data**: Per-patch population counts, geographic coordinates (lat/lon)
- **Birth/death data**: Crude birth rates per patch (or population-averaged)
- **Case data** (for calibration): Historical incidence time series per patch
- **Distance matrix**: Pairwise geodesic distances between patches (km). Can be computed with `distance()` from `laser.core.migration`.

---

### Step 2: Build the Geographic Scenario

The scenario is a GeoDataFrame with one row per spatial patch, containing population, geometry, and initial compartment counts.

```python
scenario = gpd.GeoDataFrame(cells, crs="EPSG:4326")
# Required columns: nodeid, name, population, geometry
# Initial conditions for endemic equilibrium:
scenario["E"] = 0
scenario["I"] = 3  # Seed infectious in every patch
scenario["R"] = np.round(0.95 * scenario.population).astype(np.uint32)
scenario["S"] = scenario.population - scenario["E"] - scenario["I"] - scenario["R"]
```

---

### Step 3: Seasonal Transmission

#### Primary approach: Built-in seasonality (v1.0.0)

The `SEIR.Transmission` (= `TransmissionSE`) class accepts a `seasonality` parameter directly. No custom component is needed for seasonal forcing:

```python
from laser.generic.utils import ValuesMap

# Build the 365-day Bjornstad seasonal profile (see Step 4)
# Then tile it across all simulation ticks:
nticks = 40 * 365
season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, len(scenario))

# Pass directly to the built-in Transmission component:
SEIR.Transmission(model, expdurdist, seasonality=seasonality)
```

> **Note:** Variable names `expdurdist` and `infdurdist` match the LASER API parameter names for exposed-duration and infectious-duration distributions respectively.

The built-in transmission already handles spatial coupling via `model.network` — it computes `transfer = ft[:, None] * model.network` and adjusts force of infection by incoming/outgoing pressure.

#### Advanced: Custom SeasonalTransmission

For cases requiring non-standard seasonal logic (e.g., `tick % 365` cycling, per-node seasonal profiles, or modified spatial coupling), a custom `SeasonalTransmission` component is provided in `scripts/custom_components.py`. See that file for the implementation.

#### Importation

**Built-in option:** LASER v1.0.0 provides `Infect_Random_Agents` and `Infect_Agents_In_Patch` in `laser.generic.importation`. These read parameters (`importation_period`, `importation_count`, `importation_start`, `importation_end`) from `model.params`.

**Custom option (recommended for epidemiological accuracy):** The `Importation` class in `scripts/custom_components.py` specifically targets **susceptible** agents only, which is more epidemiologically precise than the built-in classes that infect random agents regardless of state.

```python
# Built-in (simpler, infects random agents):
from laser.generic.importation import Infect_Random_Agents

# Custom (targets susceptibles only):
from custom_components import Importation
```

---

### Step 4: Configure Parameters and Seasonal Forcing

The seasonal forcing profile from Bjornstad et al. (2002) captures 26 biweekly periods of school-term-driven transmission intensity:

```python
log_betas = np.array([
    0.155, 0.571, 0.46, 0.34, 0.30, 0.34, 0.24, 0.15,
    0.31, 0.40, 0.323, 0.238, 0.202, 0.203, 0.074,
    -0.095, -0.218, -0.031, 0.433, 0.531, 0.479, 0.397,
    0.444, 0.411, 0.291, 0.509
])
beta_season = np.repeat(log_betas, int(np.floor(365 / len(log_betas))))
beta_season = np.append(beta_season, beta_season[-1])
beta_season = np.exp(beta_season - np.mean(beta_season))
```

#### Generalizing Seasonal Forcing

The Bjornstad profile above is one example of **school-term-driven forcing**. Other transmission contexts require different seasonal shapes:

- **Climate-driven forcing** (warm/wet season peak): Cosine or piecewise profiles centered on peak transmission month
- **Behavioral forcing** (school terms, holiday travel): Step functions or smoothed term profiles
- **No seasonality**: `ValuesMap.from_scalar(1.0, nticks, nnodes)`

**General recipe:** Build a 365-day multiplier array, normalize to `mean == 1.0`, tile across nticks, and wrap in a `ValuesMap`:

```python
# Example: peaked-season profile (peak at day 200, amplitude 0.3)
days = np.arange(365)
peak_day = 200
season_365 = 1.0 + 0.3 * np.cos(2 * np.pi * (days - peak_day) / 365)
season_365 /= season_365.mean()  # Normalize to mean == 1.0

season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)
```

---

### Step 5: Gravity Migration Network

Spatial coupling follows a gravity law: $M_{i,j} = k \cdot p_j^b / d_{ij}^c$ (with source exponent `a=0` fixed).

> **Auto-configuration (v1.0.0):** If `gravity_k`, `gravity_a`, `gravity_b`, `gravity_c` are in params, `Model.__init__()` will compute the distance matrix from scenario centroids and set up `model.network` automatically. Use the manual approach below when you need custom normalization.

**Manual setup (for custom normalization):**

```python
model.network = gravity(
    np.array(scenario.population), distances,
    1, 0, model.params.gravity_b, model.params.gravity_c
)
# Normalize so k represents average export fraction directly
average_export_frac = np.mean(model.network.sum(axis=1))
model.network = model.network / average_export_frac * model.params.gravity_k
model.network = row_normalizer(model.network, 0.2)  # Cap at 20% export per node
```

**Alternative migration models (v1.0.0):** Besides gravity, the framework also offers `competing_destinations()`, `stouffer()`, and `radiation()` models in `laser.core.migration`.

---

### Step 6: Assemble and Run

> **Component ordering matters:** Components execute in list order each tick. Susceptible and Recovered components should wrap the transition steps to preserve the `S + E + I + R = N` population invariant.

```python
parameters = PropertySet({
    "prng_seed": 4, "nticks": 40 * 365,
    "exp_shape": 40, "exp_scale": 0.25,   # Exposed ~10 days (gamma)
    "inf_mean": 8, "inf_sigma": 2,         # Infectious ~8 days (normal)
    "beta": 3.5,
    "cbr": average_cbr,
    "gravity_k": 0.01, "gravity_b": 0.5, "gravity_c": 1.5,
    "capacity_safety_factor": 3.0,
})

expdurdist = dists.gamma(shape=parameters.exp_shape, scale=parameters.exp_scale)
infdurdist = dists.normal(loc=parameters.inf_mean, scale=parameters.inf_sigma)

# Build seasonality ValuesMap (primary approach)
nticks = parameters.nticks
nnodes = len(scenario)
season_tiled = np.tile(beta_season, nticks // 365 + 1)[:nticks]
seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

model = Model(scenario, parameters, birthrates=birthrate_map.values)
# NOTE: birthrate_map.values must be per-1000/year (typical range 10-50).
# Same for deathrates if using MortalityByCDR. See Critical Gotchas.
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    Importation(model, infdurdist, period=30, count=3, end_tick=10 * 365),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
    MortalityByEstimator(model, estimator=survival),
    # Vaccination options (do NOT use ImmunizationCampaign — see Critical Gotchas):
    # RoutineImmunizationEx(model, period=7, coverage=0.85, age=270),
    # VaccinationCampaign(model, period=180, coverage=coverage_array),
]
# Set up gravity network (Step 5) then:
model.run("Simulation")
```

---

### Step 7: Wrap as BaseModel for Calibration

After building and testing the LASER model (Steps 1-6), wrap it inside calabaria's `BaseModel` for structured calibration, scenario management, and cloud scaling.

**Why:** calabaria provides structured parameter spaces (with bounds and transforms), Optuna-based optimization (ask/tell loop), scenario management (`@model_scenario`), and optional cloud scaling via modelops.

**Install:** `pip install modelops-calabaria`

#### Bridge Pattern: LASER Model inside BaseModel

```python
from calabaria import BaseModel, model_output, model_scenario, ScenarioSpec
from calabaria.parameters import ParameterSpace, ParameterSpec
from calabaria.parameters import ConfigurationSpace, ConfigSpec

class MySpatialSEIR(BaseModel):
    PARAMS = ParameterSpace([
        ParameterSpec("beta", lower=2.0, upper=6.0, kind="float", doc="..."),
        ParameterSpec("gravity_k", lower=1e-4, upper=0.1, kind="float", doc="..."),
        # ... additional uncertain parameters
    ])
    CONFIG = ConfigurationSpace([
        ConfigSpec("nticks", default=7300, doc="Simulation duration (days)"),
        # ... additional fixed settings
    ])

    def __init__(self, scenario_gdf, distances, birthrates, pyramid):
        super().__init__()
        self.scenario_gdf = scenario_gdf
        # ... store pre-built data from Steps 1-5

    def build_sim(self, params, config):
        # Construct LASER Model + PropertySet + gravity network + seasonality
        # + components using params and config values
        return model  # LASER Model object

    def run_sim(self, state, seed):
        laser.core.random.seed(seed)
        state.run()

    @model_output("weekly_incidence")
    def weekly_incidence(self, state):
        # Extract post-burn-in weekly incidence → pl.DataFrame
        ...
```

A complete disease-agnostic template is in `scripts/laser_basemodel.py`. Customize the component list, parameter ranges, and output extractors for your disease.

**Next step:** Use the `modelops-calabaria` skill for calibration workflow (Sobol sweeps, Optuna optimization, scenario analysis).

---

### Step 8: Quick Verification Run

Before calibration, verify the wrapped model produces sensible output with a single test run:

```python
model = MySpatialSEIR(scenario_gdf, distances, birthrates, pyramid)

# Run with plausible parameter values
outputs = model.simulate(
    {"beta": 3.5, "gravity_k": 0.01, "gravity_b": 0.5,
     "gravity_c": 1.5, "seasonal_amplitude": 1.0},
    seed=42,
)

# Check outputs
print(outputs["weekly_incidence"].head())
print(f"Total cases: {outputs['weekly_incidence']['cases'].sum()}")

# Verify non-zero incidence and spatial variation
by_patch = outputs["weekly_incidence"].group_by("patch").agg(
    pl.col("cases").sum()
)
print(by_patch)
```

If total cases are zero or all patches are identical, revisit Steps 3-6 (seasonal forcing, gravity network, initial conditions).

---

## Key Concepts

- **Critical Community Size (CCS)**: The minimum population for sustained endemic transmission without stochastic fadeout. Disease-specific (e.g., ~300K-500K for measles).
- **Traveling Waves**: Epidemics propagate from large cities to smaller populations, producing distance-dependent phase lags in wavelet analysis.
- **Gravity Model**: Spatial coupling scales with destination population and inversely with distance. The `a=0` convention means source population does not affect outward flow.
- **Seasonal Forcing**: Transmission oscillations driven by climate, behavior, or school terms. See Step 4 for general recipe.
- **BaseModel Bridge**: Wrapping a LASER Model inside calabaria's `BaseModel` enables structured calibration via Optuna, scenario management, and cloud scaling.

---

## Bundled Resources

- **`scripts/custom_components.py`** - `Importation` (susceptible-targeted seeding), `VaccinationCampaign` (correct state-based vaccination with correlated missedness), and `SeasonalTransmission` (advanced customization example) component classes
- **`scripts/calibration_metrics.py`** - CCS logistic fitting, wavelet phase similarity scoring, combined ranking, and `compute_calibration_loss()` bridge for calabaria `TrialResult`
- **`scripts/laser_basemodel.py`** - Disease-agnostic `SpatialSEIRModel(BaseModel)` template for wrapping LASER models for calabaria calibration
- **`references/laser_api_reference.md`** - Complete LASER v1.0.0 API documentation (Model, LaserFrame, PropertySet, all component variants, vital dynamics, migration models, distributions)
- **`references/wavelet_analysis.md`** - Wavelet transform functions, phase difference computation, and traveling wave detection workflow

---

## Troubleshooting

1. **No births occurring / population static**: Birthrates are likely in wrong units. Must be **per-1000/year** (typical range 10–50). Check with: `assert np.all(birthrates >= 1) and np.all(birthrates <= 60)`. Also verify `capacity_safety_factor` ≥ 2 — if capacity equals initial population, `LaserFrame.add()` silently returns no new slots.
2. **`ValueError` from `LaserFrame.add()`**: Population growth exceeded pre-allocated capacity. Increase `capacity_safety_factor` (try 3–4 for high-growth populations or long simulations).
3. **Vaccination has no effect on transmission**: `ImmunizationCampaign` / `RoutineImmunization` only set `susceptibility = 0`, which Transmission kernels ignore. Use `RoutineImmunizationEx` or `VaccinationCampaign` (custom) which set `state = RECOVERED`.
4. **All epidemics die out**: `beta` or `importation_count` may be too low, or initial `I` may be zero in some patches. Verify importation is active during burn-in period.
5. **No spatial structure**: `model.network` may be all zeros. Verify `gravity_k` is non-trivially small (try 0.001–0.1).
6. **Out of memory**: Reduce `capacity_safety_factor` or `nticks`. Each agent consumes memory for all state arrays.
7. **Wavelet NaN**: Time series with insufficient non-zero values should be padded using `pad_data()` from the wavelet reference.
8. **BaseModel `model_output` returns wrong type**: All `@model_output` methods must return `pl.DataFrame`. Use `polars`, not numpy arrays or pandas.

---

## References

- [LASER Documentation](https://laser.idmod.org/laser-generic/)
- [Grenfell et al. (2001) - Travelling waves in measles](https://www.nature.com/articles/414716a)
- [Bjornstad et al. (2002) - Estimating transmission rates](https://doi.org/10.1890/0012-9615(2002)072[0169:DOMEES]2.0.CO;2)
- [Conlan et al. (2010) - Waiting time distributions and measles persistence](https://pmc.ncbi.nlm.nih.gov/articles/PMC2842776/)
- [LASER GitHub - laser-core](https://github.com/laser-base/laser-core)
- [LASER GitHub - laser-generic](https://github.com/laser-base/laser-generic)
