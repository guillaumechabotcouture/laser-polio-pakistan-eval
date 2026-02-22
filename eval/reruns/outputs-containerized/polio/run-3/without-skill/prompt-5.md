The model ran successfully. Here's a summary of what was built and the results:

---

## Model: `/home/evaluser/pakistan_polio_model.py`

### Architecture (LASER-style components)

| Component | Class | Role |
|---|---|---|
| State storage | `Patches` | Float64 S/E/I/R/N arrays per district — analogous to LASER's `LaserFrame` |
| Transmission | `TransmissionComponent` | Stochastic binomial SEIR + gravity-weighted spatial FoI |
| Demography | `DemographyComponent` | Poisson births, binomial deaths, OPV RI at birth |
| Interventions | `SIAComponent` | Mass campaigns every 182 days |
| Seeding | `ImportationComponent` | 3 infections every 60d for years 0–5 |

### Key design choices

**Spatial coupling** — Force of infection in patch *i*:
```
λ_i(t) = β(t) × [(1−ε)·(I_i/N_i)  +  ε·Σ_j W_ij·(I_j/N_j)]
```
where `W` is a row-normalised gravity matrix (`sqrt(P_i)·sqrt(P_j) / d_ij²`) and `ε = 0.015`.

**Seasonality** — `β(t) = β₀ × [1 + 0.30·cos(2π(DOY−228)/365)]` — peak mid-August, trough mid-February.

**Transition probabilities** — Exact exponential waiting-time expressions (`p = 1−exp(−rate)`) rather than Euler approximations, matching LASER's discrete-time stochastic convention.

### Simulation results (10-year analysis window, post burn-in)

- The combined 80% RI + biannual SIA + importation-seeding regime drives the system to near-elimination very rapidly (infectious counts drop from ~22,000 at year 0 to ~12 by year 5)
- **Overall zero-incidence proportion: 0.998** across all 521 district-weeks
- **22 total cases** over the 10-year analysis window, driven entirely by residual importation pressure and stochastic re-ignition
- Population grows from ~19.7M to ~26.2M over 20 years (Pakistan ~1.9%/yr net growth)

### Outputs
- `district_summary.csv` — 20×8 per-district metrics
- `weekly_incidence.csv` — 521 weeks × 20 districts matrix (2025-12-29 to 2035-12-17)
