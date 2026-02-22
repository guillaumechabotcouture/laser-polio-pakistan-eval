The model is working correctly. The plots show the expected behavior for an R0~5 respiratory disease:

- **(A)** All 4 patches experience synchronized epidemics peaking around day 20-30, with Newmarket (200k, largest) having the highest peak
- **(B)** A sharp epidemic wave peaking at ~23,000 new infections/day, resolving by day ~100
- **(C)** Susceptible fraction drops from 90% to ~2% across all patches, crossing below the herd immunity threshold (1/R0 = 0.2)
- **(D)** Clean SEIR compartment dynamics — S depletes rapidly, brief E and I waves, then nearly everyone is in R

**Key results:**
- **88.4% attack rate** across all patches (consistent with R0=5 and 90% initially susceptible)
- **Peak on day 24** with ~23,000 daily infections
- **468,648 total infections** out of 530,000 population
- All compartments remain non-negative and sum to population throughout

The complete code is at `eastland_respiratory_model.py` with all specifications:
- 4 patches (100k, 200k, 150k, 80k) with gravity-model spatial coupling
- R0~5 (`beta = 0.5 = R0/D_infectious`)
- Latent period: gamma(shape=4, scale=1) = 4 day mean
- Infectious period: gamma(shape=5, scale=2) = 10 day mean
- Initial: 90% S, 0% E, 1% I, 9% R
- Components: Susceptible, Exposed, Infectious, Recovered, Transmission (in correct execution order)
- 365-day simulation
