The complete model is written to `/tmp/eastland_seir.py`. Here's a full walkthrough of every design decision:

---

## File: `/tmp/eastland_seir.py`

### Structure & API choices

**Scenario GeoDataFrame** — LASER's `Model` requires a GeoDataFrame with columns `nodeid`, `population`, `geometry`, plus initial compartment columns `S`, `E`, `I`, `R`. Point geometries at the given coordinates (lon, lat) in EPSG:4326 are used; `Model.__init__` calls `get_centroids()` internally to extract lat/lon for the distance matrix.

**Parameters** — `PropertySet` is a dot-access dict from `laser.core`. `gravity_k/a/b/c` are set here for transparency (the Model reads them), and then the network is explicitly rebuilt below to apply the correct `row_normalizer` threshold.

**`calc_capacity`** — called automatically inside `Model.__init__` when `birthrates` is passed. It computes how many agent slots to pre-allocate to accommodate 10 years of Poisson births (CBR=30/1000/yr × 10 yr × population) with a safety factor. No explicit user call needed.

**Gravity network** — after `Model()` returns, `model.network` is overridden:
```python
dist_matrix   = distance(lats, lons, lats, lons)          # Haversine
raw_network   = gravity(populations, dist_matrix, k=0.01, a=1, b=1, c=1.5)
model.network = row_normalizer(raw_network, 0.15)
```

**Seasonal forcing** — `ValuesMap.from_timeseries(data_1yr, nnodes, nticks=NTICKS)` cycles a 365-day cosine series over 3650 ticks. The actual signature is `(data, nnodes, nticks=None)` — `nnodes` is the positional second argument. The seasonality ValuesMap is passed to `SEIR.Transmission`.

**`BirthsByCBR`** — `pyramid=None, track=False` disables age-structure tracking (no `AliasedDistribution` needed). The birth step still fires daily and Poisson-samples newborns into S.

**`MortalityByCDR`** — takes a `(nticks, nnodes)` array of CDR values. Converts internally: `p(death/tick) = 1 − (1 − CDR/1000)^(1/365)`.

**`PulseImportation`** — custom component following the `__init__(model) + step(tick)` contract:
```python
class PulseImportation:
    def __init__(self, model, patch=0, count=3, period=60): ...
    def step(self, tick):
        if tick == 0 or tick % self.period != 0: return
        # find susceptibles in patch, S→I, update S/I counts at tick+1
```
It seeds directly into **I** (not E) to represent an imported *infectious* case, sampling `itimer ~ gamma(5,2)`. The compartment arrays at `tick+1` are updated after `model._initialize_flows` has already seeded `tick+1 = tick`, so deltas accumulate correctly.

**Component order**:
```
Susceptible → Exposed → Infectious → Recovered
→ Transmission → BirthsByCBR → MortalityByCDR → PulseImportation
```
