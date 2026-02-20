# Pakistan Polio SEIR Model: Calibration Report

## Executive Summary

A 13-district agent-based SEIR model for WPV1 transmission in Pakistan was
built using the LASER framework (v1.0.0) and calibrated against GPEI/MMWR
surveillance data (2020-2025). The model predicts **~20 expected paralytic
cases per year** nationally, compared to **~36 observed** (ratio 0.57,
within 2x). Province-level patterns match epidemiological reality: KP and
Sindh show sustained endemic transmission, Punjab is disease-free, and
Balochistan's per-capita rate exceeds observed provincial averages.

The model was built iteratively over a single Claude Code session using the
`laser-spatial-disease-modeling` skill, which provided correct LASER API
signatures, component ordering conventions, and spatial coupling patterns.
Several critical bugs in the initial implementation were discovered and fixed
during calibration, including a fundamental error in how birth/death rate
units were passed to LASER's vital dynamics components.

---

## 1. Model Description

### Districts

13 districts representing Pakistan's key epidemiological zones:

| District | Province | Population | Unreachable | RI Coverage | SIA Coverage |
|----------|----------|:---:|:---:|:---:|:---:|
| Quetta | Balochistan | 1.1M | 35% | 55% | 60% |
| Peshawar | KP | 2.0M | 15% | 78% | 72% |
| Karachi | Sindh | 3.5M | 16% | 80% | 74% |
| Hyderabad | Sindh | 1.7M | 16% | 78% | 72% |
| Bannu | KP/South | 0.9M | 28% | 50% | 55% |
| N. Waziristan | KP/FATA | 0.55M | 40% | 30% | 40% |
| S. Waziristan | KP/FATA | 0.6M | 40% | 30% | 40% |
| DI Khan | KP/South | 1.5M | 30% | 45% | 50% |
| Lahore | Punjab | 2.5M | 2% | 95% | 90% |
| Islamabad | ICT | 0.5M | 3% | 92% | 87% |
| Multan | Punjab | 1.8M | 2% | 93% | 88% |
| Faisalabad | Punjab | 2.2M | 2% | 94% | 89% |
| Rawalpindi | Punjab | 1.8M | 2% | 93% | 88% |

Total modeled population: 20.65M (initial), growing to ~31.8M over 20 years.

### Disease Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| R0 | 6 | Standard WPV1 estimate |
| Beta (transmission rate) | 0.2143 (= R0 / infectious period) | Derived |
| Latent period | ~3 days (gamma, shape=3, scale=1) | Polio incubation |
| Infectious period | ~28 days (normal, mean=28, sd=3) | Fecal shedding duration |
| Crude birth rate | 29 per 1000/year | Pakistan DHS |
| Crude death rate | 7 per 1000/year | Pakistan vital statistics |
| Seasonal forcing | Cosine, peak Sep 2, amplitude +/-30% | Monsoon seasonality |

### Key Model Components

1. **SEIR compartmental dynamics** — Built-in LASER Susceptible, Exposed,
   Infectious, Recovered components with standard transitions.

2. **Gravity-model spatial coupling** — k=0.005, b=0.5, c=1.5, row-normalized
   at 20% max export. Computed from Haversine distances between district
   centroids.

3. **Correlated-missedness vaccination** — Custom `PerPatchVaccination`
   component with a per-agent `reachable` flag (int8). Unreachable agents
   (representing zero-dose children in hard-to-reach, refusal, or
   inaccessible communities) are permanently excluded from both RI and SIA.
   This models the critical real-world dynamic where the same children are
   missed by every campaign.

4. **Endemic corridor importation** — Custom `PatchImportation` seeds 3
   infections per endemic-corridor district every 30 days, targeting
   susceptible agents only. Endemic patches: Quetta, Peshawar, Karachi,
   Hyderabad, Bannu, N/S Waziristan, DI Khan.

5. **Vital dynamics** — `BirthsByCBR` (29/1000/yr) with Pakistan age pyramid
   and `MortalityByCDR` (7/1000/yr). Population grows at ~2.2%/year.

### Simulation Configuration

- Duration: 20 years (7,300 ticks)
- Burn-in: First 10 years discarded; analysis covers years 10-20
- Initial conditions: Per-district immunity fractions from MMWR, age
  distribution from Pakistan age pyramid (exponential decay, mean 65yr)

---

## 2. Comparison with MMWR/GPEI Surveillance Data

### Observed WPV1 Paralytic Cases (2020-2025)

Source: CDC MMWR Vol. 69-73, Pakistan NEOC (endpolio.com.pk), WHO IHR
Emergency Committee reports, GPEI fact sheets.

| Year | Total | Balochistan | KP | Sindh | Punjab | ICT |
|------|:---:|:---:|:---:|:---:|:---:|:---:|
| 2020 | 84 | 26 | 22 | 22 | 14 | 0 |
| 2021 | 1 | 1 | 0 | 0 | 0 | 0 |
| 2022 | 20 | 0 | 20 | 0 | 0 | 0 |
| 2023 | 6 | 0 | 4 | 2 | 0 | 0 |
| 2024 | 74 | 27 | 22 | 23 | 1 | 1 |
| 2025 | ~31 | 0 | ~20 | 9 | 1 | 0 |
| **Avg** | **36** | **9.0** | **14.7** | **9.3** | **2.7** | **0.2** |

### Model-to-Observed Translation

Model infections are converted to expected paralytic cases using the standard
WPV1 paralysis-to-infection ratio of **1:200** (CDC Yellow Book; WHO Polio
Fact Sheet). This ratio reflects WPV1's high paralytogenicity relative to
types 2 and 3 (~1:1000).

### Province-Level Comparison

| Province | Model infections/yr | Model paralytic/yr | Observed avg/yr | Ratio |
|----------|:---:|:---:|:---:|:---:|
| **KP** | **2,464** | **12.3** | **14.7** | **0.84** |
| **Sindh** | **1,384** | **6.9** | **9.3** | **0.74** |
| Balochistan | 229 | 1.1 | 9.0 | 0.13 |
| Punjab | 1 | ~0 | 2.7 | ~0 |
| ICT | 0 | 0 | 0.2 | — |
| **National** | **4,078** | **20.4** | **35.8** | **0.57** |

### Per-Capita Rate Comparison

Because our model represents high-risk subsets of each province (not full
provincial populations), per-capita rates in the model are expected to exceed
province-wide observed rates:

| Province | Model districts pop | Model rate (per M) | Province pop | Observed rate (per M) |
|----------|:---:|:---:|:---:|:---:|
| Balochistan | 1.1M | 1.04 | 14M | 0.64 |
| KP | 5.6M | 2.22 | 40M | 0.37 |
| Sindh | 5.2M | 1.33 | 50M | 0.19 |
| Punjab | 8.3M | ~0 | 120M | 0.02 |

Model per-capita rates are 1.6-7x higher than province-wide rates, which is
epidemiologically correct: the modeled districts are the known endemic
hotspots within each province.

### District-Level Detail

| District | Infections/yr | Expected paralytic/yr | IQR (paralytic) |
|----------|:---:|:---:|:---:|
| N. Waziristan | 989 | 4.9 | [1.1 - 5.6] |
| S. Waziristan | 905 | 4.5 | [1.0 - 5.4] |
| Karachi | 873 | 4.4 | [3.6 - 4.9] |
| Hyderabad | 511 | 2.6 | [2.1 - 2.8] |
| Peshawar | 320 | 1.6 | [1.1 - 2.0] |
| Quetta | 229 | 1.1 | [0.4 - 0.9] |
| DI Khan | 143 | 0.7 | [0.4 - 1.0] |
| Bannu | 108 | 0.5 | [0.3 - 0.8] |

S/N equilibrium values: N/S Waziristan at 0.18 (above 1/R0 = 0.167,
self-sustaining), Quetta at 0.166 (at endemic equilibrium), Karachi and
Hyderabad at 0.16 (near-equilibrium, sustained by importation/coupling).
Punjab districts at 0.03-0.04 (disease-free, R_eff ~ 0.2).

---

## 3. Assessment of Model Calibration

### What the model captures well

1. **Province-level transmission hierarchy** — KP > Sindh > Balochistan >
   Punjab, matching GPEI surveillance patterns. Punjab is correctly
   disease-free.

2. **Endemic equilibrium dynamics** — Karachi and N/S Waziristan show
   S/N oscillating around 1/R0 with periodic outbreak peaks, the hallmark
   of endemic SEIR dynamics.

3. **Monsoon seasonal forcing** — Clear infection peaks during Jul-Oct monsoon
   season, matching Pakistan's polio seasonality.

4. **Inter-annual variability** — The model produces year-to-year variation
   in case counts (IQR spans 2-5x), consistent with the extreme observed
   variability (1 case in 2021 vs 84 in 2020).

5. **National scale** — 20 expected paralytic cases/year vs 36 observed is
   within 2x, which is reasonable for a 13-district subset model.

### Where the model underperforms

1. **Balochistan absolute count** — Quetta produces only ~1 paralytic/year
   vs ~9 observed for all of Balochistan. This is primarily a coverage
   limitation: we model 1.1M of 14M people. Adding Killa Abdullah, Pishin,
   and Jaffarabad districts would likely close this gap.

2. **Punjab sporadic cases** — The model produces zero Punjab cases, while
   real data shows occasional cases (14 in 2020, 1 in 2024). This reflects
   spillover events from endemic areas that our gravity coupling (k=0.005)
   is too weak to generate.

3. **Campaign-driven dynamics** — Real inter-annual variability is driven
   largely by campaign effectiveness variation (mass SIA reach fluctuates
   with security access, funding, and political will). The model uses
   constant vaccination parameters and cannot capture these exogenous shocks.

---

## 4. Key Technical Findings During Development

### Critical Bug: Birth/Death Rate Unit Mismatch

LASER's `BirthsByCBR` and `MortalityByCDR` expect rates in **per-1000/year**
format (they divide by 1000 internally). The initial implementation passed
pre-converted daily per-capita rates (CBR/1000/365), which were interpreted
as near-zero rates. This caused:

- **Zero population growth** (capacity = initial population, no free slots for births)
- **No susceptible replenishment** from births
- **Artificially stable infection dynamics** from the initial susceptible pool
  depleting over time

The bug was not immediately obvious because the model still produced plausible-
looking disease dynamics from the large initial susceptible pool. It was only
detected when population tracking showed exactly 18M at every time step with
`agent capacity = active count = 18,000,000`.

**Impact:** This bug fundamentally altered the model's susceptible dynamics.
Without proper births, vaccination had no ongoing effect (no newborns to
vaccinate), and endemic equilibrium was driven by initial conditions rather
than demographic replenishment. Results from this phase were qualitatively
wrong despite appearing reasonable.

### Critical Bug: TransmissionSE State Check

LASER's built-in `TransmissionSE` kernel only checks `state == SUSCEPTIBLE`
(int8 value 0) when selecting agents for infection. It does **not** check
the `susceptibility` property. This means the built-in `ImmunizationCampaign`
and `RoutineImmunization` components (which set `susceptibility = 0`) have
**zero effect on transmission**.

Vaccination must set `state = RECOVERED` (value 3) to actually prevent
infection, and must also update `nodes.S[tick+1]` and `nodes.R[tick+1]`
for compartment bookkeeping.

### Structural Issue: Independent Bernoulli SIA Model

The initial vaccination model applied SIA as independent Bernoulli draws
per agent per campaign. With 2 campaigns/year over 5 years (10 draws),
even modest per-campaign coverage (e.g., 40%) produces near-complete
cumulative coverage: (1-0.40)^10 = 0.6%. This massively overestimates
real SIA effectiveness because:

- The same hard-to-reach children are missed by every campaign
- Vaccine refusal communities are systematically excluded
- Doses administered to already-immune children are wasted

The fix was to add a per-agent `reachable` flag, permanently set at birth
based on per-district unreachable fractions. Only reachable agents can
receive RI or SIA. This correctly models the structural immunity gap that
sustains endemic polio transmission in Pakistan.

### Minor Bug: PatchImportation Broadcasting

The initial importation component subtracted `n_infect` from ALL patches'
S counts rather than just the target patch. This was fixed by using
per-patch indexing: `nodes.S[tick+1, patch_id] -= n_infect`.

---

## 5. Assessment of LASER Skill Effectiveness

### What the skill provided

1. **Correct API signatures** — Import paths (`from laser.generic import
   SEIR, Model`), component constructors, distribution creation
   (`dists.gamma(shape=, scale=)`), gravity model arguments. Without the
   skill, these would likely be hallucinated or wrong.

2. **Component ordering convention** — The skill's explicit ordering
   (Susceptible → Exposed → Infectious → Recovered → Births → Deaths →
   Transmission) preserves the S+E+I+R=N population invariant. Getting
   this wrong produces compartment accounting errors.

3. **GeoDataFrame scenario pattern** — The required columns (nodeid, name,
   population, geometry, S, E, I, R) and the `Model(scenario, params)`
   constructor pattern were correctly specified.

4. **Gravity network setup** — The manual normalization pattern
   (`gravity()` → normalize by average export → `row_normalizer()`) was
   directly usable.

5. **ValuesMap for seasonality** — The
   `ValuesMap.from_timeseries(season_tiled, nnodes)` pattern for passing
   seasonal forcing to `SEIR.Transmission` was correct.

### Where the skill was insufficient or misleading

1. **Birth/death rate units** — The skill's example uses `birthrate_map.values`
   without specifying the expected unit format. The `calc_capacity` and
   `BirthsByCBR` functions both expect per-1000/year, but this was not
   documented in the skill. This led to the most severe bug in the
   implementation.

2. **Vaccination state semantics** — The skill does not mention that
   `TransmissionSE` ignores the `susceptibility` property and only checks
   `state == SUSCEPTIBLE`. This required reading the LASER source code to
   discover, and is critical for any model with vaccination.

3. **Measles bias** — The skill is built around the England & Wales measles
   example, using school-term seasonal forcing (Bjornstad 26 biweekly
   periods), CCS analysis, and wavelet phase similarity. These don't apply
   to polio. The domain adaptation (monsoon forcing, paralysis ratio,
   vaccination structure) had to be done entirely from epidemiological
   knowledge outside the skill.

4. **Agent capacity** — The skill mentions `capacity_safety_factor` but
   doesn't explain how it interacts with `calc_capacity` or what happens
   when capacity is exhausted (births are silently dropped).

### Overall Assessment

The skill was **essential for API correctness** (imports, signatures,
component patterns) and **structural correctness** (component ordering,
scenario format, network setup). Without it, significant time would have
been spent debugging hallucinated API calls.

The skill was **not helpful for domain adaptation** — polio epidemiology
required knowledge beyond the measles-centric skill content. The correlated
missedness model, endemic corridor importation, and paralysis ratio
translation were all domain-specific additions.

The skill had a **blind spot on vital dynamics units** that caused the most
impactful bug in the implementation. This suggests the skill should include
explicit documentation of expected rate formats for `BirthsByCBR`,
`MortalityByCDR`, and `calc_capacity`.

---

## 6. Limitations

1. **13-district subset** — The model covers 20.65M of Pakistan's ~240M
   population. Province-level absolute counts are not directly comparable
   to national surveillance data.

2. **Constant vaccination parameters** — Real campaign coverage fluctuates
   with security access, funding, and political conditions. The model
   cannot capture campaign-driven inter-annual variability.

3. **No waning immunity** — The SEIR model treats recovery as permanent
   immunity. Real polio immunity can wane, especially mucosal immunity from
   OPV, requiring booster doses.

4. **No age-structured transmission** — Polio transmission is concentrated
   in under-5s, but the model's force of infection applies equally to all
   ages. Age-structured contact patterns would improve realism.

5. **Gravity coupling may be too weak** — The model produces zero Punjab
   cases despite occasional real-world spillover. Increasing `gravity_k`
   or adding explicit transportation corridors could address this.

6. **Population growth dilution** — With 2.2%/year growth, the population
   denominator changes significantly over the 10-year analysis window. The
   summary statistics use initial population, slightly inflating per-capita
   rates.

---

## 7. Figures

- `polio_seir_diagnostics.png` — 6-panel overview: seasonal forcing,
  incidence rates, R_eff, S/N trajectories, active infections, post-burn-in
  weekly incidence with monsoon shading.

- `polio_dynamics_detail.png` — 6-panel deep dive: annual incidence heatmap,
  S/N+I dual-axis for Peshawar and Quetta, monthly incidence rate,
  cumulative incidence, force of infection for endemic districts.

- `polio_mmwr_comparison.png` — 4-panel MMWR comparison: annual paralytic
  cases (model IQR vs observed), province-level bar chart, per-capita rates,
  model time series by province.

---

## 8. Reproducibility

```bash
# Install dependencies
pip install laser-generic geopandas

# Run the model (requires ~3 minutes, ~4GB RAM)
/opt/anaconda3/bin/python3 scripts/polio_seir_10patch.py

# Outputs:
#   Console: summary statistics and MMWR comparison table
#   eval/outputs/polio_seir_diagnostics.png
#   eval/outputs/polio_dynamics_detail.png
#   eval/outputs/polio_mmwr_comparison.png
```

Key source files:
- `scripts/polio_seir_10patch.py` — Main simulation, plotting, and comparison
- `scripts/custom_components.py` — PerPatchVaccination (with correlated
  missedness) and PatchImportation components

---

## References

- CDC MMWR Vol. 69 No. 46 (2020): Progress Toward Poliomyelitis Eradication — Pakistan, Jan 2019-Sep 2020
- CDC MMWR Vol. 71 No. 42 (2022): Pakistan, Jan 2021-Jul 2022
- CDC MMWR Vol. 72 No. 33 (2023): Pakistan, Jan 2022-Jun 2023
- CDC MMWR Vol. 73 No. 36 (2024): Pakistan, Jan 2023-Jun 2024
- Pakistan NEOC (endpolio.com.pk): Year-end case tallies 2024-2025
- WHO 43rd Polio IHR Emergency Committee (Nov 2025)
- CDC Yellow Book: Poliomyelitis — paralysis-to-infection ratio 1:200 for WPV1
- LASER Framework: https://laser.idmod.org/laser-generic/
