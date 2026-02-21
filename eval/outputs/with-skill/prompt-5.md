The model runs successfully and produces realistic endemic dynamics. Here's a summary of the complete script at `scripts/polio_seir_20district.py`:

## Model Specifications (all met)

| Spec | Implementation |
|---|---|
| 20 districts, 50k-2M | Chaman (50k) to Karachi (2M), total 17.85M |
| SEIR, R0=6 | beta=6/28=0.2143, gamma(3,1) latent, normal(28,3) infectious |
| Gravity coupling | `gravity()` + `row_normalizer()`, k=0.005, b=0.5, c=1.5 |
| Monsoon forcing | Cosine peak at day 245 (Sep), ±30% amplitude |
| RI at 80% | OPV at 6 weeks (ri_age=42), 80% coverage among reachable |
| SIA every 6 months | 80% coverage, ages 0-5 years |
| Importation | 3 infections/60 days (1 each in Quetta, Peshawar, N_Waziristan) for first 5 years |
| 20-year sim, 10yr burn-in | nticks=7300, analysis period years 10-20 |

## Key Results (post burn-in, years 10-20)

**Endemic corridor districts** sustain transmission:
- DI Khan: 7,224 infections/yr, **0% zero-incidence weeks**
- N Waziristan: 5,923/yr, 0.4% zero weeks
- Quetta: 6,509/yr, 4.6% zero weeks

**Well-vaccinated districts** see rare spillovers:
- Islamabad: 3/yr, **96.2% zero weeks**
- Jacobabad: 0/yr, 99.4% zero weeks
- Punjab districts: 60-70% zero weeks

**Outputs**:
- `eval/outputs/weekly_incidence_20district.csv` — 10,420 rows (521 weeks x 20 districts)
- `eval/outputs/polio_20district_diagnostics.png` — 6-panel diagnostic figure
- Console: zero-incidence week analysis per district

The diagnostics show monsoon-season peaks (panel E), clear spatial heterogeneity in the heatmap (panel D), and R_eff oscillating around 1.0 in endemic districts (panel C).
