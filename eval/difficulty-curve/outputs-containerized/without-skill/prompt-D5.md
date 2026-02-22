The complete script has been written to `/home/evaluser/eastland_seir.py`. Here is a walkthrough of every design decision:

---

**Imports**
```
laser.core.PropertySet             – model parameter container
laser.core.migration.distance      – Haversine pairwise-distance matrix
laser.core.migration.gravity       – gravity flow matrix builder
laser.core.migration.row_normalizer – makes rows ∈ [0, 1]
laser.core.utils.calc_capacity     – projects required agent-array capacity
laser.core.demographics.AliasedDistribution – age pyramid sampler
laser.core.distributions           – Numba-compiled duration distributions
laser.generic.Model                – main simulation driver
laser.generic.SEIR.*               – Susceptible, Exposed, Infectious,
                                     Recovered, Transmission
laser.generic.vitaldynamics.*      – BirthsByCBR, MortalityByCDR
```

**Key parameter decisions**

| Item | Value | Derivation |
|---|---|---|
| `beta` | 0.5 day⁻¹ | R0 / mean_infectious = 5 / (5×2) |
| `nticks` | 3 650 | 10 years × 365 days |
| `birthrates` / `deathrates` | (3650 × 4) arrays | CBR=30, CDR=10 per 1 000/yr |
| `calc_capacity` | called explicitly | sizes agent arrays for birth growth over 10 yr |
| `distance(LATS, LONS, LATS, LONS)` | 4×4 km matrix | Haversine, lat first |
| `gravity(POPS, dist, k=0.01, a=1, b=1, c=1.5)` | raw flow matrix | then `row_normalizer` clamps row-sums |
| `model.network = norm_network` | overrides auto-built network | explicit per spec |
| `pyramid` | uniform 80-bucket `AliasedDistribution` | synthetic; used only to back-assign initial dates-of-birth |
| `exp_dist` | `gamma(shape=4, scale=1)` | latent period ~4 days |
| `inf_dist` | `gamma(shape=5, scale=2)` | infectious period ~10 days |

**Component order rationale**

`MortalityByCDR` runs before `BirthsByCBR` because `BirthsByCBR.step()` sizes the birth cohort against the post-mortality population at `tick+1` to avoid double-counting.

**Output**
After `model.run()`, final SEIR compartment counts are read from `model.nodes.S/E/I/R[NTICKS, patch_index]` and printed for each of the 4 patches.
