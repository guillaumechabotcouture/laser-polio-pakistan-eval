The model runs correctly and all three features are working as expected. Here's a summary of what was written in `scripts/polio_gravity_seasonal.py`:

**1. Gravity-model migration network** (`build_gravity_network`, lines 100-142)
- Uses `laser.core.migration.gravity()` with source exponent `a=0`, destination exponent `b=0.5`, distance decay `c=1.5`
- Nearest neighbors (50km) have the strongest coupling (~0.0017); the farthest pair (450km) has ~0.00006
- The banded structure visible in panel (B) reflects the linear patch arrangement

**2. Monsoon seasonal forcing** (`build_monsoon_seasonality`, lines 148-185)
- Cosine profile centered on day 245 (Sep 2): `1.0 + 0.3 * cos(2π(t-245)/365)`
- Peak = 1.3x during Jul-Oct (monsoon), trough = 0.7x during Dec-Mar (dry season)
- Mean is exactly 1.0, so baseline R0 is preserved
- Panel (C) shows incidence peaks aligning with monsoon shading

**3. Row-normalized network** (`row_normalizer` with 0.15 cap, line 138)
- All row sums are well below 15% (max = 0.57%), so the cap is satisfied
- The gravity weights are small enough that the cap doesn't bind in this configuration

**Key outputs:**
- ~15,000 infections/year across 1M total population
- R_eff oscillates between 0.78–1.25 (mean ~1.04), consistent with endemic equilibrium
- Interior patches (4–6) have slightly higher incidence due to more neighbors in the gravity network
- Population grows from 1M to 1.24M over 10 years (CBR 29 - CDR 7 = 22/1000 net growth)
