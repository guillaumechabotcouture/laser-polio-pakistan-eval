The model ran successfully. Here's the complete script and results:

**Script:** `eastland_respiratory_model.py`

**Results summary:**
- **Total infections:** 468,608 / 530,000 (~88.4% attack rate) — consistent with R0≈5 theory (final size for R0=5 is ~96% in a fully susceptible population; with 90% initial susceptibility, ~88% is correct)
- **Peak:** ~14,042 new infections on day 37
- **All 4 patches** show nearly identical attack rates (~88.2–88.5%), with the gravity network providing spatial coupling
- **Epidemic burns through** the susceptible pool by ~day 150, then susceptible fraction stabilizes near 1/R0 = 0.2

**Key code decisions:**
- `beta = 0.5` (R0 / mean_infectious = 5 / 10)
- `SEIR.Exposed(model, latent_dist, inf_dist)` — latent gamma(4,1) for E→I transition, infectious gamma(5,2) for I timer assignment
- `SEIR.Transmission(model, inf_dist)` — infectious period gamma(5,2) as `expdurdist` parameter per your specification
- Gravity network auto-computed from `gravity_k=0.02, gravity_b=0.5, gravity_c=1.5`
- No vital dynamics, vaccination, or seasonality (clean 1-year epidemic)
