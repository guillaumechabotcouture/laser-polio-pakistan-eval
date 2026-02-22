# Suite 1: Atomic Competency Map — Scoring Rubric

## Scoring Scale (per prompt)

Each micro-test is scored 0/1/2:

| Score | Meaning | Criteria |
|-------|---------|----------|
| **0** | Wrong/missing | The API call is absent, uses an invented API, or has a fundamentally wrong signature |
| **1** | Partial | Recognizes the right API but has wrong arguments, wrong units, or is missing a key detail |
| **2** | Correct | Valid LASER v1.0.0 API usage that would work if executed |

**Maximum score per condition:** 40 (20 prompts x 2 points each)

## Per-Prompt Check Items

### A01 — Model(scenario, params)
- [ ] Imports `Model` from `laser.generic` (not invented path)
- [ ] Imports `PropertySet` from `laser.generic`
- [ ] Creates a GeoDataFrame with population, S, E, I, R, geometry columns
- [ ] `Model(scenario, params)` — two positional args, GeoDataFrame first
- **Score 2:** All four checks pass
- **Score 1:** Right classes but wrong constructor args or missing scenario columns
- **Score 0:** Invented API or no Model instantiation

### A02 — GeoDataFrame scenario
- [ ] Uses `geopandas.GeoDataFrame` (not plain DataFrame)
- [ ] Includes all required columns: population, S, E, I, R, geometry
- [ ] S/E/I/R counts are integers derived from population (not fractions)
- [ ] Geometry uses `shapely.geometry.Point` objects
- **Score 2:** Valid GeoDataFrame with all columns and proper types
- **Score 1:** Has the columns but wrong types (fractions instead of counts, or missing geometry)
- **Score 0:** Plain DataFrame or missing critical columns

### A03 — PropertySet
- [ ] Imports `PropertySet` from `laser.generic`
- [ ] Uses keyword arguments: `PropertySet(nticks=..., beta=..., ...)`
- [ ] Shows attribute access (e.g., `params.beta` or `params.nticks`)
- [ ] Accepts arbitrary custom parameters (seasonal_amplitude)
- **Score 2:** Correct import, keyword construction, and access pattern
- **Score 1:** Right class but wrong construction pattern (e.g., dict-style)
- **Score 0:** Invented class or fundamentally wrong usage

### A04 — dists.gamma
- [ ] Imports `dists` from `laser.generic`
- [ ] Uses `dists.gamma(shape=5, scale=2)` — named parameters shape and scale
- [ ] Does NOT use scipy-style `(a=, loc=, scale=)` parameterization
- [ ] Mean = shape * scale = 10 (mathematically consistent)
- **Score 2:** Correct import path and parameter names
- **Score 1:** Right module but wrong parameter names (e.g., `a=5, b=2`)
- **Score 0:** Uses scipy directly or invents a distribution API

### A05 — dists.normal
- [ ] Imports `dists` from `laser.generic`
- [ ] Uses `dists.normal(loc=28, scale=5)` or equivalent
- [ ] Parameters: `loc` for mean, `scale` for std deviation
- **Score 2:** Correct import and parameter names
- **Score 1:** Right module but wrong parameter names
- **Score 0:** Uses scipy directly or invents API

### A06 — gravity network
- [ ] Imports `gravity` from `laser.generic`
- [ ] Calls `gravity(pops, dists, k, a, b, c)` — exactly 6 positional/keyword args
- [ ] Constructs a proper distance matrix (4x4, symmetric, zero diagonal)
- [ ] All 6 parameters present: populations, distances, k, a, b, c
- **Score 2:** Correct 6-argument signature with proper distance matrix
- **Score 1:** Right function but wrong number of args or missing distance matrix
- **Score 0:** Invents gravity function or uses completely wrong API

### A07 — row_normalizer
- [ ] Imports `row_normalizer` from `laser.generic`
- [ ] Calls `row_normalizer(network, max_frac)` — 2 args
- [ ] `max_frac` is 0.20 (not 20, not a percentage)
- [ ] Applied to the gravity network output
- **Score 2:** Correct import, 2-arg call, fraction (not percentage)
- **Score 1:** Right function but wrong arg format (e.g., passes 20 instead of 0.2)
- **Score 0:** Invents normalization function or manual implementation only

### A08 — ValuesMap.from_timeseries
- [ ] Imports `ValuesMap` from `laser.generic`
- [ ] Uses `ValuesMap.from_timeseries(array, num_patches)` class method
- [ ] Input array is 1D with 365 elements
- [ ] `num_patches` parameter is 4
- **Score 2:** Correct class method call with proper array and patch count
- **Score 1:** Right class but wrong construction method or wrong array shape
- **Score 0:** Invents ValuesMap API or no ValuesMap usage

### A09 — SEIR.Transmission
- [ ] Uses `SEIR.Transmission(model, ...)` or equivalent component addition
- [ ] Passes a gamma distribution for the infectious period duration
- [ ] Includes seasonality parameter (ValuesMap)
- [ ] Distribution created with `dists.gamma(shape=5, scale=5.6)`
- **Score 2:** Correct Transmission setup with distribution and seasonality
- **Score 1:** Right component but missing distribution or seasonality
- **Score 0:** Invents Transmission API or no component setup

### A10 — SEIR component ordering
- [ ] Adds all five components: Susceptible, Exposed, Infectious, Recovered, Transmission
- [ ] Order is S, E, I, R, then Transmission (state-trackers before dynamics)
- [ ] Explains that order matters (state accounting before transmission)
- [ ] Components added via `SEIR.Susceptible(model)`, etc.
- **Score 2:** All 5 components in correct order with valid explanation
- **Score 1:** All components present but wrong order or no explanation
- **Score 0:** Missing components or invented component names

### A11 — BirthsByCBR (units)
- [ ] Imports `BirthsByCBR` from `laser.generic`
- [ ] Passes CBR as 44 (per-1000/year), NOT 0.044 or 0.00012
- [ ] Explicitly states units are per-1000/year
- [ ] States that BirthsByCBR divides by 1000 internally
- **Score 2:** Correct import, value=44, and explicit unit explanation
- **Score 1:** Right class but wrong units (e.g., converts to daily per-capita first)
- **Score 0:** Invents birth API or fundamentally wrong usage

### A12 — MortalityByCDR (units)
- [ ] Imports `MortalityByCDR` from `laser.generic`
- [ ] Passes CDR as 13 (per-1000/year), NOT 0.013 or 0.0000356
- [ ] Explicitly states units are per-1000/year
- [ ] States that MortalityByCDR divides by 1000 internally
- **Score 2:** Correct import, value=13, and explicit unit explanation
- **Score 1:** Right class but wrong units (pre-converted)
- **Score 0:** Invents mortality API or fundamentally wrong usage

### A13 — calc_capacity (units)
- [ ] Imports `calc_capacity` from `laser.generic`
- [ ] Passes birthrate as 44 (per-1000/year), consistent with BirthsByCBR
- [ ] Passes initial population and simulation duration
- [ ] Explains that calc_capacity expects per-1000/year birthrate
- **Score 2:** Correct import, per-1000/year units, proper function call
- **Score 1:** Right function but wrong units or missing parameters
- **Score 0:** Invents capacity function or no calc_capacity usage

### A14 — Custom component pattern
- [ ] Class has `__init__(self, model)` method
- [ ] Class has `step(self, tick)` method
- [ ] `__init__` stores reference to model (`self.model = model`)
- [ ] `step` uses `tick % 90 == 0` for periodic importation
- [ ] Modifies model state to add infections (sets agent state to INFECTIOUS)
- **Score 2:** Both methods with correct signatures, working importation logic
- **Score 1:** Right method signatures but incomplete or broken logic
- **Score 0:** Wrong method signatures or no component class

### A15 — AliasedDistribution
- [ ] Imports `AliasedDistribution` from `laser.generic`
- [ ] Passes probability weights (e.g., `[0.30, 0.25, 0.20, 0.15, 0.10]`)
- [ ] Shows sampling (e.g., `dist.sample(1000)` or equivalent)
- [ ] Weights sum to 1.0
- **Score 2:** Correct import, proper weights, working sample call
- **Score 1:** Right class but wrong construction or sampling API
- **Score 0:** Invents distribution class or uses numpy/scipy instead

### A16 — distance (haversine)
- [ ] Imports `distance` from `laser.generic`
- [ ] Passes lat/lon coordinates (not lon/lat)
- [ ] Function returns distance in kilometers
- [ ] Uses correct coordinate values for Lahore and Karachi
- **Score 2:** Correct import and coordinate order
- **Score 1:** Right function but wrong coordinate order or wrong import
- **Score 0:** Implements haversine manually or uses geopy instead

### A17 — State enum values
- [ ] Shows SUSCEPTIBLE = 0
- [ ] Shows EXPOSED = 1
- [ ] Shows INFECTIOUS = 2
- [ ] Shows RECOVERED = 3
- [ ] Correct import path (from laser.generic or from SEIR-related module)
- **Score 2:** All four values correct with valid import
- **Score 1:** Some values correct but others wrong or guessed
- **Score 0:** Invents state system or all values wrong

### A18 — model.network manual setup
- [ ] Creates Model instance first
- [ ] Computes gravity network with `gravity()` function
- [ ] Assigns to `model.network` (correct attribute name)
- [ ] Network is a 2D numpy array matching patch count
- **Score 2:** Correct attribute assignment with gravity-computed network
- **Score 1:** Right approach but wrong attribute name or missing gravity step
- **Score 0:** No network assignment or completely invented API

### A19 — Vaccination gotcha
- [ ] Sets agent `state` to RECOVERED (value 3), NOT just susceptibility
- [ ] Explicitly mentions that TransmissionSE checks state, not susceptibility
- [ ] Selects susceptible agents (state == 0) before vaccination
- [ ] Samples 80% of susceptible agents
- **Score 2:** Correctly sets state=RECOVERED and explains the gotcha
- **Score 1:** Mentions state change but also/instead modifies susceptibility
- **Score 0:** Only modifies susceptibility field (the common mistake)

### A20 — Seasonal normalization
- [ ] Creates a 365-element numpy array
- [ ] Uses sinusoidal or piecewise seasonal pattern
- [ ] Peak is approximately 1.5x, trough approximately 0.5x
- [ ] Annual mean is approximately 1.0 (verified in code)
- **Score 2:** Correct array with verified mean ~1.0 and correct amplitude
- **Score 1:** Seasonal array but mean is not ~1.0 or amplitude is wrong
- **Score 0:** No seasonal array or fundamentally wrong approach

## Key Discriminators

These prompts are expected to show the **largest skill advantage** (most likely
to separate WITH vs WITHOUT conditions):

| Prompt | Why it discriminates |
|--------|---------------------|
| **A06** | `gravity()` has a non-obvious 6-arg signature |
| **A07** | `row_normalizer` is an obscure utility function |
| **A08** | `ValuesMap.from_timeseries` is a class method, not a constructor |
| **A11** | Per-1000/year units are a common mistake without documentation |
| **A12** | Same units trap as A11 |
| **A13** | Same units trap, compounded with calc_capacity |
| **A15** | `AliasedDistribution` is a niche class unlikely in training data |
| **A19** | The state-vs-susceptibility gotcha requires insider knowledge |

These prompts are expected to show the **smallest skill advantage** (Claude's
base knowledge may suffice):

| Prompt | Why it may not discriminate |
|--------|---------------------------|
| **A02** | GeoDataFrame construction is standard geopandas |
| **A03** | PropertySet is a simple keyword-arg container |
| **A17** | SEIR state values 0-3 are conventional |
| **A20** | Seasonal arrays are standard numpy |

## Score Sheet Template

### Condition A: WITH skill

| Prompt | Score (0-2) | Notes |
|--------|-------------|-------|
| A01 | | |
| A02 | | |
| A03 | | |
| A04 | | |
| A05 | | |
| A06 | | |
| A07 | | |
| A08 | | |
| A09 | | |
| A10 | | |
| A11 | | |
| A12 | | |
| A13 | | |
| A14 | | |
| A15 | | |
| A16 | | |
| A17 | | |
| A18 | | |
| A19 | | |
| A20 | | |
| **Total** | **/40** | |

### Condition B: WITHOUT skill

| Prompt | Score (0-2) | Notes |
|--------|-------------|-------|
| A01 | | |
| A02 | | |
| A03 | | |
| A04 | | |
| A05 | | |
| A06 | | |
| A07 | | |
| A08 | | |
| A09 | | |
| A10 | | |
| A11 | | |
| A12 | | |
| A13 | | |
| A14 | | |
| A15 | | |
| A16 | | |
| A17 | | |
| A18 | | |
| A19 | | |
| A20 | | |
| **Total** | **/40** | |

### Summary Comparison

| Category | Prompts | WITH | WITHOUT | Delta |
|----------|---------|------|---------|-------|
| Core objects | A01, A02, A03 | /6 | /6 | |
| Distributions | A04, A05, A15 | /6 | /6 | |
| Spatial coupling | A06, A07, A16, A18 | /8 | /8 | |
| SEIR components | A09, A10, A14 | /6 | /6 | |
| Demographics | A08, A11, A12, A13, A20 | /10 | /10 | |
| State semantics | A17, A19 | /4 | /4 | |
| **Overall** | **A01-A20** | **/40** | **/40** | |

## Expected Outcome Hypotheses

| Category | Expected skill advantage | Rationale |
|----------|------------------------|-----------|
| Core objects | Moderate (+1-2) | Model/PropertySet are discoverable; GeoDataFrame columns are not |
| Distributions | Moderate (+1-2) | dists.gamma param names differ from scipy; AliasedDistribution is obscure |
| Spatial coupling | **Strong (+3-4)** | gravity() 6-arg signature and row_normalizer are highly skill-dependent |
| SEIR components | Moderate (+1-2) | Component ordering is documented in skill but also somewhat conventional |
| Demographics | **Strong (+3-5)** | Per-1000/year unit trap is the #1 gotcha; ValuesMap API is non-obvious |
| State semantics | **Strong (+2-3)** | Vaccination gotcha requires insider knowledge of TransmissionSE internals |

**Overall prediction:** WITH skill scores 30-38/40; WITHOUT skill scores 12-22/40.
The largest gap should appear in Demographics and Spatial coupling categories.
