The model runs successfully. Here's a summary of the results:

**Model architecture** — follows the LASER framework patterns:
- `LaserFrame`-style state arrays for compartments (S/E/I per species per district)
- Component pipeline: importation → force of infection → stochastic SEIS transitions
- Gravity-model spatial mixing matrix for inter-district water network coupling
- Binomial draws for all stochastic transitions (exact discrete-time: `p = 1 - exp(-rate)`)

**Key results over 10 years:**

| Metric | Value |
|---|---|
| Total human cases | 3,400 |
| Total dog infections | 864 |
| Human cases year 1 → 10 | 1,215 → 1 |
| Dog infections year 1 → 10 | 266 → 0 |
| Elimination achieved? | Not within 10 years (1 human case in year 10) |

**Dynamics** — Cases decline steeply as interventions scale up. Dog infections reach zero by year 9, but 1 residual human case in year 10 prevents formal elimination (requires 3 consecutive zero-case years). The model is approaching elimination and would likely achieve it within 1-2 additional years.

The code is at `guinea_worm_model.py` — a single self-contained file runnable with `python3 guinea_worm_model.py` using only numpy.
