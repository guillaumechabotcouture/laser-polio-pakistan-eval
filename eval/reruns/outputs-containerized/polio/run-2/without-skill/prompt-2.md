The file is at `/home/evaluser/polio_spatial_model.py`. Here is a walkthrough of every design decision.

---

## Structure

```
polio_spatial_model.py
├── seasonal_multiplier()          # pure function
├── build_gravity_network()        # pure function
├── SeasonalTransmissionComponent  # LASER component
├── ProgressionRecoveryComponent   # LASER component
├── ReporterComponent              # LASER component
└── PatchyPolioModel               # model class
    ├── LaserFrame  patches        # patch-level data store
    ├── ndarray     G              # (n,n) coupling matrix
    └── list        components     # ordered step pipeline
```

---

## Feature 1 — Gravity-model migration network

`build_gravity_network()` (`polio_spatial_model.py:108`)

```python
G_raw[i, j] = N_j**b / d_ij**c      # b=0.5, c=1.5
```

- Distances are computed analytically from `positions_km` — no geographic library needed.
- The diagonal is set to `+inf` before the power-law division (no self-export), then zeroed.
- All arithmetic is vectorised over the full `(n, n)` matrix.

---

## Feature 2 — Seasonal forcing

`seasonal_multiplier()` (`polio_spatial_model.py:60`)

A single cosine wave, centred on **1 September** (DOY 243, the midpoint of Jul–Oct):

```
multiplier(t) = 1.0 + 0.3 * cos(2π(t − 243)/365)
```

| Season window | DOY range | β multiplier |
|---|---|---|
| Jul–Oct peak | 182–304 | ≤ 1.3× |
| Dec–Mar trough | 335–90 | ≥ 0.7× |

Applied each tick: `beta_t = beta_base * seasonal_multiplier(doy)`.

---

## Feature 3 — Row-normalised coupling (≤ 15 % FOI export)

Back in `build_gravity_network()`:

```python
row_sums = G.sum(axis=1, keepdims=True)
G = np.where(row_sums > 0, G * (max_export / row_sums), G)
```

This guarantees every row sums to exactly `0.15`. Each patch retains the complementary `local_ret[i] = 1 − G[i,:].sum() = 0.85` of its own FOI.

The transmission step (`polio_spatial_model.py:179`) applies both parts:

```python
local_foi    = local_ret * prevalence          # 85% home patch
imported_foi = G.T @ prevalence               # 15% distributed by gravity
total_foi    = beta_t * (local_foi + imported_foi)
```

---

## LASER API used

| LASER call | Purpose |
|---|---|
| `LaserFrame(capacity=n, initial_count=n)` | Fixed patch container |
| `.add_scalar_property(name, dtype, default)` | Per-patch current state (S/E/I/R counts) |
| `.add_vector_property(name, length=T, dtype)` | Per-patch time-series output — shape `(T, n)` |
| `PropertySet({...})` | Typed parameter dictionary with attribute access |

---

## Running the model

```bash
pip install laser-core        # one-time, if not already installed
python polio_spatial_model.py
```

Default: 10 patches, 500 000 people, 3-year simulation. Change `num_patches`, `total_population`, `timesteps`, or any epidemiological/gravity parameter in the `PropertySet` block at the bottom of the file.
