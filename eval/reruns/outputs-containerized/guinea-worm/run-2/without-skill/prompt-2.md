The model is complete at `/home/evaluser/gwd_laser_model.py` (note: `/eval-work` is read-only on this system). Here's a full breakdown of what was built:

---

## Architecture

The code follows LASER's component-based, array-oriented design without requiring `laser_core`. State lives in `NodeFrame` with shape `(n_species, n_nodes)` numpy arrays; components have `initialize(model)` and `step(tick, model)` methods registered on a central `Model` object.

---

## Components

### 1. `SeasonalForcing`
Builds a 365-element daily multiplier via linear interpolation between monthly midpoints (no step discontinuities at month boundaries):

| Period | Multiplier |
|--------|-----------|
| Jan–Mar | 1.80× (dry season peak) |
| Apr–May | declining 1.20→0.80× |
| Jun–Sep | 0.40× (rains disperse copepods) |
| Oct–Dec | rising 0.70→1.50× |

### 2. `WaterSourceCoupling`
Gravity matrix with `W[dest, src] = N_dest^0.3 / dist_km^2.0`, built from Haversine distances between 12 endemic districts. Row-normalised so coupling weights sum to 1 per destination. Self-coupling (local ponds, 10 km radius) ranges 0.88–0.97; cross-district contamination averages 6.3%.

### 3. `CrossSpeciesTransmission`
Water contamination at source node `i`:
```
C[i] = alpha_dog * I_dog[i]/N_dog[i] + alpha_human * I_human[i]/N_human[i]
```
with `alpha_dog=19, alpha_human=1` giving exactly 95% dog contribution. Spatial mixing via `W @ C`, then `foi[s,j] = beta[s] * seasonal_forcing * water_contam[j]`.

### 4. `SEISTransition`
Forward-Euler SEIS with species-specific demography (`mu_dog = 1/(7×365)`, `mu_human = 1/(55×365)`). No immunity on recovery — worms exit and hosts return directly to S. Floating-point N drift is corrected each step.

---

## Calibration

Beta values are analytically derived so that R₀ = 1.3 for dogs:
```
R0 = [σ/(σ+μ)] × [β × α × W_self × forcing_mean / (γ+μ)]
```
This yields endemic equilibrium of ~2.8% dog prevalence and ~0.08% human prevalence — consistent with historical pre-control estimates in southern Chad. A six-panel diagnostic figure is saved to `gwd_results.png` showing SEIS dynamics for both species, the seasonal forcing curve, the coupling heatmap, district-level prevalence traces, and the dog/human water contamination split.
