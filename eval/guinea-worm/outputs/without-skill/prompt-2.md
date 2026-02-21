All tasks complete. Here's a summary of what was built:

## Guinea Worm SEIS Transmission Model — Complete

**33 tests passing, 6 modules, fully functional simulation**

### Architecture (mirrors LASER framework)

| Module | Purpose |
|--------|---------|
| `spatial.py` | Haversine distance + gravity-model coupling (`k*pop_i^a * pop_j^b / dist^c`) |
| `seasonal.py` | Chad dry/rainy season forcing (1.8x peak Jan-Mar, 0.4x Jun-Sep) with cosine interpolation |
| `population.py` | SEIS compartment container per species per district |
| `components.py` | Force-of-infection (cross-species), exposure, incubation, recovery step functions |
| `model.py` | Simulation runner with step-function registration |
| `guinea_worm_model.py` | Main assembly: 8 Chad districts, dogs + humans, 2-year simulation |

### Three requested features implemented:

1. **Water-source mediated transmission** — Gravity model spatial coupling with distance exponent c=2.0 and destination population exponent b=0.3. Row-normalized coupling matrix ensures each district's total external water access is bounded.

2. **Dry-season transmission forcing** — Monthly multipliers for Chad's climate: peak 1.8x (Jan-Mar stagnant pools), declining 1.2x (Apr-May), low 0.4x (Jun-Sep flowing water), baseline 1.0x (Oct-Dec). Cosine-interpolated daily for smooth transitions.

3. **Cross-host transmission** — Dogs (95% weight) and humans (5% weight) contaminate shared water bodies. Dogs have higher beta (0.08 vs 0.03) and higher per-capita prevalence, matching Chad's epidemiology where ~95% of guinea worm cases originate from dogs.

### SEIS dynamics confirmed
- S→E→I→S cycle with no lasting immunity
- Population conservation (S+E+I=N) verified at every timestep
- No negative compartments
- Spatial spread from seeded districts to all 8 districts
- Dogs show higher per-capita prevalence than humans
