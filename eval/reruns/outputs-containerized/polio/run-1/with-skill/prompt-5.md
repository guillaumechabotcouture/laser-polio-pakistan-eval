The model ran successfully. Here's a summary of what was built and the key results:

---

## Model: `/tmp/pakistan_polio_model.py`

**Code structure** (~340 lines, 9 functions):

| Section | What it does |
|---|---|
| `DISTRICTS` | 20 Pakistan districts, pop 50k–2M, with per-district `unreachable_frac` |
| `build_monsoon_seasonality()` | Cosine profile peak day 245 (early Sep), ±30% amplitude, mean=1.0 |
| `build_scenario()` | GeoDataFrame seeded at endemic equilibrium (R/N≈83%) |
| `initialize_ages()` | Assigns pyramid-drawn ages to initial population (prevents cohort artifact) |
| `run_model()` | Builds gravity network, assembles 9 components, runs 7300-tick simulation |
| `compute_weekly_incidence()` | Aggregates `newly_infected` to weekly per-district post-burn-in |
| `analyze_zero_incidence()` | Tabulates weeks with zero cases per district |
| `print_summary()` | Annual incidence, S/N at end, compartment integrity checks |
| `plot_diagnostics()` | 6-panel figure (seasonal forcing, incidence, R_eff, heatmap, weekly, zero-weeks) |

**Key epidemiological results (post-burn-in, years 10–20):**

| District | Pop | Total cases | % Zero weeks |
|---|---|---|---|
| DI_Khan | 250k | 10,318 | **56%** (most endemic) |
| N_Waziristan | 100k | 10,541 | **70%** |
| Quetta | 200k | 8,416 | **66%** |
| S_Waziristan | 100k | 7,058 | **69%** |
| Karachi | 2M | 304 | **82%** |
| Lahore | 700k | 75 | **90%** |
| Islamabad | 110k | 8 | **98%** |

**What drives the pattern:** Districts with `unreachable_frac > 1/R0 ≈ 17%` (the herd-immunity gap) can sustain local chains even with 80% nominal coverage. The tribal/border districts (40–55% unreachable) form the endemic corridor. Punjab and ICT (2–5% unreachable) approach elimination.

**Three caveats on memory/scale:**
1. Total district population was scaled to ~6.2M (keeping the 50k–2M range) because the system has only ~1 GB available RAM. The epidemiological structure is preserved.
2. Importation stops at year 5 as specified; you can see from panel E that transmission then depends purely on endogenous dynamics in high-unreachable districts.
3. Output files: `/tmp/pakistan_polio_weekly_incidence.csv` (521 weeks × 20 districts), `/tmp/pakistan_polio_zero_weeks.csv`, `/tmp/pakistan_polio_20district.png`.
