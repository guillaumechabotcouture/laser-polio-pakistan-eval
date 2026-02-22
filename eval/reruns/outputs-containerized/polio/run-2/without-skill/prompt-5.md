The model runs successfully. Here is a complete summary of what was built and the key design decisions.

---

## Pakistan Spatial Polio Model — Complete

**File:** `/tmp/pakistan_polio_model.py`

### Architecture (LASER-style)
The code mirrors the LASER framework conceptually — patch-level columnar state arrays (`S`, `E`, `I`, `R`), discrete process functions called each time step, and a separation between model mechanics and output aggregation. Since `laser_core` is not installed, the implementation uses numpy/pandas directly.

### Key epidemiological parameters

| Parameter | Value | Notes |
|---|---|---|
| R0 | 6.0 | Baseline |
| Latent period | 3 days (σ = 0.333/d) | E→I |
| Infectious period | 28 days (γ = 0.0357/d) | I→R |
| Seasonal forcing | ±30%, peak DOY 258 | Mid-Sep monsoon peak |
| OPV birth dose | 80% admin × 65% efficacy → **52% effective** | Per-dose seroconversion realistic for Pakistan |
| SIA (accessible) | 70% × 65% → **45.5% effective** | 12 districts |
| SIA (hard-to-reach) | 15% × 65% → **9.8% effective** | 8 KP/Balochistan districts |
| Waning immunity | 3%/yr (R→S) | OPV mucosal immunity loss |
| Gravity coupling | α=1, decay=300 km | Normalised to 0.5% mixing |

### Why OPV efficacy was added (design rationale)
With 100% per-dose efficacy, 80% OPV + SIA every 6 months gives R_eff < 1 and instant elimination — mathematically correct but not useful for planning. In Pakistan, per-dose OPV seroconversion is ~65% due to enteric virus interference, heat-chain gaps, and maternal antibody interference. This raises R_eff to 2.6 under routine OPV alone, sustaining endemic transmission.

### Results (post-burn-in, years 10–20)

```
Hard-to-reach districts    → endemic (0%–25% zero-incidence weeks)
Accessible districts       → eliminated (100% zero-incidence weeks)
National zero-incidence wks → 0.0%   (always active in Peshawar/Swat)
National peak weekly cases  → 10.3
National mean weekly cases  → 3.7
```

This matches the historical Pakistan pattern: persistent reservoir in KP/FATA border areas, elimination achieved in urban Punjab/Sindh.

### Outputs
```
/tmp/pakistan_polio_model.png           4-panel figure
/tmp/weekly_incidence_by_district.csv   520 weeks × 20 districts
/tmp/weekly_incidence_per100k.csv       same, per 100k
/tmp/district_summary.csv              per-district summary stats
/tmp/model_metadata.csv                all parameters + scalar results
```
