The model runs successfully. All verifications pass:

- **Population growth**: 530k → 644k (21.6% over 10 years, consistent with net 20/1000/yr growth)
- **All compartments non-negative**: OK
- **470k total infections** across 10 years, spread across all 4 patches
- **Spatial coupling active**: all 4 patches have infections

One observation: at tick 1825 (midpoint), E=0 and I=0 — the initial epidemic burned through the susceptible pool and the disease went extinct. This is expected for a high-R0 disease (R0=5) in small populations without importation: the initial epidemic depletes susceptibles, then birth-replenished susceptibles accumulate but there's no reintroduction. The first-wave dynamics and demographic components are all working correctly.

Here's a summary of the key design decisions following the LASER skill:

| Element | Implementation |
|---------|---------------|
| **Beta** | `R0 / inf_period = 5 / 10 = 0.5` |
| **Latent period** | `dists.gamma(shape=4, scale=1)` → mean 4 days |
| **Infectious period** | `dists.gamma(shape=5, scale=2)` → mean 10 days |
| **Birth rates** | `30.0` per-1000/year (LASER convention, NOT daily per-capita) |
| **Death rates** | `10.0` per-1000/year via `MortalityByCDR` |
| **Capacity** | `calc_capacity` called internally by `Model()` via `birthrates=` parameter |
| **Transmission** | `SEIR.Transmission(model, expdurdist)` — takes the *exposed* duration distribution (assigns etimer to newly infected) |
| **Component order** | Susceptible → Exposed → Infectious → Recovered → Transmission → Births → Deaths |
| **Gravity network** | Auto-computed from `gravity_k/a/b/c` params + scenario centroids |
