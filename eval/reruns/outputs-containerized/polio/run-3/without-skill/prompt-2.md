The complete model is at `/home/evaluser/polio_spatial_model.py` (note: `/eval-work` is read-only). It runs without errors and all mathematical properties verify correctly.

---

## Architecture

The code follows the LASER design pattern exactly — `laser-core` is not installed so the idioms are implemented directly:

```
PatchModel
├── Pre-allocated float64 arrays: S, E, I, R, N  (shape: n_patches)
├── coupling matrix W                              (shape: n×n)
├── params dict                                    (shared state)
├── step_fns list                                  (ordered components)
└── run(duration_days)                             → record → [fn(m)] → t++
```

---

## Key design decisions

### 1. Gravity coupling (`build_gravity_coupling`)

```
g[i,j]  =  N_j^0.5 / d_ij^1.5       (off-diagonal raw weights)
W[i,j]  =  0.15 × g[i,j] / Σ_k g[i,k]
```

Row-normalising to 1 and then scaling by 0.15 guarantees every row sums to exactly 0.15. The diagonal stays 0. The complementary local weight (0.85) is applied inside the FOI calculation so the total always sums to 1:

```
λ_i = β_eff × [ 0.85 · I_i/N_i  +  (W @ prev)[i] ]
                 └─ local ────┘   └─ spatial (≤15%) ┘
```

### 2. Seasonal forcing (`make_seasonal_step`)

```
mult(doy) = 1.0 + 0.3 × cos( 2π(doy − 213) / 365 )
```

- doy 213 (1 Aug): `cos(0) = 1` → mult = **1.3** (monsoon peak) ✓  
- doy 31 (31 Jan): `cos(π) ≈ −1` → mult = **0.7** (winter trough) ✓  
- Writes `model.params["beta_eff"]` each tick; the SEIR step reads it.

### 3. SEIR step (`make_seir_step`)

Forward-Euler with daily Δt. `new_exposed` is clipped to `S` to prevent negative compartments. Vital dynamics (μ) are included with Pakistan's birth rate (~35/1000/year) to replenish susceptibles.

### Parameters used

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| R0 | 6 | Endemic, low-coverage scenario |
| σ | 1/10 | 10-day latent period |
| γ | 1/14 | 14-day infectious period |
| μ | 35/1000/365 | Pakistan birth ≈ death rate |
| Seed | 20 cases, patch 1 | Single-district outbreak start |
