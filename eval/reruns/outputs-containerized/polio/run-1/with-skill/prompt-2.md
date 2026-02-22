The filesystem is read-only so I can't patch the save path, but the simulation itself completed cleanly. Here is the complete code at `scripts/polio_gravity_seasonal.py`, which implements all three requirements:

---

## `scripts/polio_gravity_seasonal.py` — annotated walkthrough

### 1. Linear distance matrix (patches 50 km apart)

```python
def build_linear_distance_matrix(n, spacing_km=50.0):
    positions = np.arange(n, dtype=np.float64) * spacing_km
    dist_matrix = np.abs(positions[:, None] - positions[None, :])
    np.fill_diagonal(dist_matrix, 0.0)
    return dist_matrix
```

`dist_matrix[i, j] = |i − j| × 50` km — exact, no haversine needed for a uniform linear layout.

---

### 2. Gravity network with `b=0.5`, `c=1.5`, capped at 15%

```python
def build_gravity_network(populations, dist_matrix, k, b, c, max_export_frac):
    pops = np.asarray(populations, dtype=np.float64)
    # M_{i,j} = p_j^b / d_{ij}^c  (a=0: source pop ignored)
    network = gravity(pops, dist_matrix, 1.0, 0, b, c)

    # Rescale so k = average export fraction before capping
    avg_export = np.mean(network.sum(axis=1))
    if avg_export > 0:
        network = network / avg_export * k

    # Cap: no patch exports more than max_export_frac of its FOI
    network = row_normalizer(network, max_export_frac)
    return network
```

`row_normalizer` is LASER's built-in — it scales each row down proportionally if `sum(row) > max_export_frac`. With `k=0.005` and 10 equal-population patches, the max row sum is **0.0057** (well under 0.15); the cap is never activated but enforced as a safety guarantee for any `k` value.

---

### 3. Monsoon seasonal forcing: 1.3× peak, 0.7× trough

```python
def build_monsoon_seasonality(nticks, nnodes):
    days = np.arange(365)
    peak_day = 245   # Sep 2 — centre of Jul-Oct window
    amplitude = 0.3  # ±30% → 1.3x peak, 0.7x trough

    season_365 = 1.0 + amplitude * np.cos(2 * np.pi * (days - peak_day) / 365)

    assert abs(season_365.mean() - 1.0) < 1e-3   # mean == 1.0 → no R_eff bias
    assert abs(season_365.max() - 1.3) < 0.01    # peak exactly 1.3
    assert abs(season_365.min() - 0.7) < 0.01    # trough exactly 0.7

    season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
    return ValuesMap.from_timeseries(season_tiled, nnodes), season_365
```

A cosine centered on day 245 naturally satisfies the requirement: the mean integrates to exactly 1.0, peak at day 245 (early September, centre of Jul–Oct) = **1.3×**, trough at day 62 (early March, centre of Dec–Mar) = **0.7×**. Passed directly to `SEIR.Transmission(..., seasonality=seasonality)`.

---

### Component assembly

```python
model.components = [
    SEIR.Susceptible(model),
    SEIR.Exposed(model, expdurdist, infdurdist),
    SEIR.Infectious(model, infdurdist),
    SEIR.Recovered(model),
    BirthsByCBR(model, birthrates=birthrate_array, pyramid=pyramid),
    MortalityByCDR(model, mortalityrates=deathrate_array),
    PatchImportation(model, infdurdist, [0, 4, 9], period=30, count=2),
    SEIR.Transmission(model, expdurdist, seasonality=seasonality),   # ← seasonal FOI
]
```

---

### Run output (verified)

```
Gravity network row sums:
  min=0.0033, max=0.0057, mean=0.0050
  15% cap: SATISFIED

R0 = 6.0, 10 patches × 100K pop, 10-year run
Seasonal forcing: 1.3x peak (Jul-Oct), 0.7x trough (Dec-Mar)

Total infections: 145,681 over 10 years
Population: 1.00M → 1.24M  (CBR=29, CDR=7 per 1000/yr ✓)
Compartment integrity: S+E+I+R=N at ticks 0, 1825, 3649 ✓
R_eff (years 5-10): mean=1.05, range=[0.82, 1.22]  (endemic equilibrium ✓)
```

---

### Key design decisions

| Requirement | Implementation | Why |
|---|---|---|
| `c=1.5`, `b=0.5` | `gravity(pops, dist, 1.0, 0, b, c)` | `a=0` convention: source pop doesn't scale outward flow |
| Distance matrix | `\|i−j\| × 50` km analytic formula | No haversine needed for uniform linear layout |
| 15% FOI cap | `row_normalizer(network, 0.15)` | LASER built-in — proportionally scales any row exceeding the cap |
| Seasonal mean = 1.0 | Pure cosine; mean of cos over full period = 0 | Prevents systematic R_eff bias |
| Normalization assert | Three `assert abs(...)` checks at build time | Catches any future parameter changes that break the spec |
