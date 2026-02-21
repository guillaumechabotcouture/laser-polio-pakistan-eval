The simulation completed successfully with epidemiologically sensible results. Here's the complete script at `scripts/polio_seir_basic_10patch.py` and a summary of what it does:

**Model Configuration:**
- **10 patches**, 100K population each (1M total), placed at realistic Pakistan coordinates
- **R0 = 6** (beta = 6/28 = 0.2143/day), latent period = 3 days (gamma), infectious period = 28 days (normal)
- **Initial conditions**: 95% recovered, 5 infectious, ~5K susceptible per patch
- **Vital dynamics**: CBR = 29/1000/yr, CDR = 7/1000/yr (Pakistan-like), giving ~2.2%/yr growth
- **Gravity-model** spatial coupling between patches
- **Monsoon seasonal forcing** (Jul-Oct peak, +/-30%)
- **Periodic importation** in 3 patches to prevent stochastic fadeout

**Key Results (10-year run):**
- ~15,000 infections/year across all patches (15,000 per million)
- Population grows from 1.0M to 1.24M (vital dynamics working correctly)
- Susceptible fraction stabilizes around S* = 1/R0 = 0.167 (endemic equilibrium)
- R_eff oscillates around 1.0 (mean 1.06) in the second half — classic endemic behavior
- Clear epidemic waves visible in years 6-10 once susceptibles rebuild sufficiently, with monsoon-season peaks

The diagnostic plot shows: (A) the monsoon seasonal forcing profile, (B) weekly incidence with growing epidemic waves, (C) susceptible fraction rising from 5% toward the 1/R0 = 16.7% equilibrium, and (D) R_eff crossing 1.0 and triggering epidemics.
