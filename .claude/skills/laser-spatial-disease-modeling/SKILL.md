---
name: laser-spatial-disease-modeling
description: This skill should be used for building spatial disease transmission models using the LASER (Light Agent Spatial modeling for ERadication) framework (v1.0.0+). This skill is appropriate when the user asks to model disease spread across geographic populations, simulate measles or other SEIR-type diseases with spatial coupling, calibrate compartmental models against historical data, perform wavelet analysis of disease time series, set up gravity-model migration networks, or reproduce spatial epidemiological phenomena like traveling waves and critical community size. Trigger phrases include "LASER model", "spatial disease model", "measles England Wales", "gravity model transmission", "SEIR spatial simulation", "wavelet phase analysis", "critical community size", or "calibrate epidemic model".
---

# Spatial Disease Transmission Modeling with LASER

## Overview

The LASER framework (v1.0.0+, December 2025) enables spatially-explicit, agent-based disease transmission modeling across geographic populations. This skill provides a complete workflow for building, calibrating, and analyzing spatial SEIR models, based on the canonical England and Wales pre-vaccine measles dataset. The workflow demonstrates reproduction of key spatiotemporal phenomena including traveling waves of infection and the relationship between disease extinction/reimportation and population size (critical community size).

**Packages:** `pip install laser-generic` (installs both `laser-core` and `laser-generic`)

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
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    Importation(model, infdurdist, period=30, count=3, end_tick=10 * 365),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
    MortalityByEstimator(model, estimator=survival),
]
# Set up gravity network (Step 5) then:
model.run("Simulation")
```

---

### Step 7: Calibration

Random parameter sampling with two fitness metrics. Reusable metric functions are in `scripts/calibration_metrics.py`.

#### Parameter Ranges

| Parameter | Range | Description |
|-----------|-------|-------------|
| `beta` | U(3, 4.5) | Transmission rate (R0 ~ 24-36) |
| `amplitude` | U(0.5, 1.5) | Seasonal forcing multiplier |
| `k` | 10^U(-3.5, -1) | Gravity coupling constant (log-uniform) |
| `b` | U(0.25, 1.0) | Destination population exponent |
| `c` | U(1, 2) | Distance decay exponent |

**Metric 1 - CCS Similarity**: Logistic curve fit to proportion of zero-incidence weeks vs. log10(population). The `fit_mean_var()` and `similarity_metric()` functions in `scripts/calibration_metrics.py` handle this.

**Metric 2 - Wavelet Phase Similarity**: Cross-wavelet phase differences from a reference city detect traveling waves. Full implementation in `references/wavelet_analysis.md`.

**Combined Ranking**: Sum of CCS and phase similarity ranks; the lowest combined rank identifies the best-fit simulation.

---

### Step 8: Analyze Results

Weekly incidence aggregation and proportion of zero-incidence weeks (after burn-in) are the primary outputs:

```python
incidence = model.nodes.newly_infected
num_weeks = incidence.shape[0] // 7
weekly_incidence = incidence[:num_weeks*7, :].reshape(
    num_weeks, 7, incidence.shape[1]
).sum(axis=1)
prop_zero = np.mean(weekly_incidence[1040:, :] == 0, axis=0)  # After 20-year burn-in
```

---

## Key Concepts

- **Critical Community Size (CCS)**: The minimum population for sustained endemic transmission without stochastic fadeout. For measles, approximately 300,000-500,000.
- **Traveling Waves**: Epidemics propagate from large cities to smaller populations, producing distance-dependent phase lags in wavelet analysis.
- **Gravity Model**: Spatial coupling scales with destination population and inversely with distance. The `a=0` convention means source population does not affect outward flow.
- **Seasonal Forcing**: School-term-driven transmission oscillations captured as 26 biweekly log-beta values from Bjornstad et al. (2002).

---

## Bundled Resources

- **`scripts/custom_components.py`** - `Importation` (susceptible-targeted seeding) and `SeasonalTransmission` (advanced customization example) component classes
- **`scripts/calibration_metrics.py`** - CCS logistic fitting, wavelet phase similarity scoring, and combined ranking functions
- **`references/laser_api_reference.md`** - Complete LASER v1.0.0 API documentation (Model, LaserFrame, PropertySet, all component variants, vital dynamics, migration models, distributions)
- **`references/wavelet_analysis.md`** - Wavelet transform functions, phase difference computation, and traveling wave detection workflow

---

## Troubleshooting

1. **Out of memory**: The `capacity_safety_factor` parameter or `nticks` should be reduced. Each agent consumes memory for state arrays.
2. **All epidemics die out**: The `beta` or `importation_count` values may be too low, or initial `I` may be zero in some patches.
3. **No spatial structure**: The `model.network` matrix may be all zeros. The `gravity_k` value should be verified as non-trivially small.
4. **Wavelet NaN**: Time series with insufficient non-zero values should be padded using `pad_data()` from the wavelet reference.
5. **Poor calibration fits**: Parameter ranges may need widening, or the burn-in period (10+ years recommended) may be insufficient.

---

## References

- [LASER Documentation](https://laser.idmod.org/laser-generic/)
- [Grenfell et al. (2001) - Travelling waves in measles](https://www.nature.com/articles/414716a)
- [Bjornstad et al. (2002) - Estimating transmission rates](https://doi.org/10.1890/0012-9615(2002)072[0169:DOMEES]2.0.CO;2)
- [Conlan et al. (2010) - Waiting time distributions and measles persistence](https://pmc.ncbi.nlm.nih.gov/articles/PMC2842776/)
- [LASER GitHub - laser-core](https://github.com/laser-base/laser-core)
- [LASER GitHub - laser-generic](https://github.com/laser-base/laser-generic)
