The model ran successfully in ~6 seconds. Here's a summary of the results:

**Verification - all checks pass:**
- Population: 530,000 → 643,820 (21.5% growth over 10 years, confirming CBR=30 > CDR=10 in correct units)
- All compartments non-negative
- All 4 patches received infections (spatial coupling working)
- Total infections: 476,021 over 10 years

**Seasonal forcing:**
- Winter peak (days 0-90): 1.300x
- Summer trough (days 150-240): 0.700x  
- Mean: 1.0000

**Key model components in `eastland_respiratory_model.py`:**

| Component | Details |
|-----------|---------|
| Patches | 4 (100k, 200k, 150k, 80k) at specified coordinates |
| Disease | R0=5, latent gamma(4,1), infectious gamma(5,2) |
| Init | 90% S, 0% E, 1% I, 9% R |
| Vital dynamics | BirthsByCBR(30/1000/yr), MortalityByCDR(10/1000/yr) |
| Gravity network | k=0.01, a=1, b=1, c=1.5, row_normalizer(0.15) |
| Seasonality | `ValuesMap.from_timeseries()` with 365-day cosine profile |
| Duration | 10 years (3,650 ticks) |

The dynamics show a large initial epidemic (starting from 90% susceptible with R0=5), which rapidly depletes the susceptible pool to the herd immunity threshold (~1/R0 = 0.2). After the epidemic, births gradually replenish susceptibles. The plot is saved at `eval/outputs/eastland_respiratory_seir.png`.
