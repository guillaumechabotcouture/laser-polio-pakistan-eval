All checks pass. The model is at `/tmp/guinea_worm_chad.py`. Here's a complete summary of what was built and what the results show:

---

## Model architecture

**LASER component pipeline** (5 components called each day):

| Component | What it does |
|---|---|
| `_component_importation` | Adds 1 infected traveller (→ E) every 90 days for 3 years |
| `_component_water` | Euler update: `dW/dt = ω_h·i_h + ω_d·i_d − δ_W·W` |
| `_component_spatial_mix` | Gravity-model diffusion of W between districts |
| `_component_foi` | Computes λ_h, λ_d with ABATE / filter / seasonality applied |
| `_component_transitions` | Binomial draws for S→E, E→I, I→S (stochastic, conserves N) |

**R0 calibration** — for water-source-mediated SEIS, R0 = β·ω / (γ·δ_W), which is population-size-independent when W is a per-capita concentration index. Verified: R0_H = 1.5000, R0_D = 2.5000.

**Dog reservoir dominance** — verified at 99.6% of water contamination from dogs on day 0 (case containment + tethering reduce dog shedding but dogs still dominate by population ratio; initial seeding is 97.6% dogs by count).

---

## 10-year results (seed = 42)

```
Year     Human    Dog    Total
─────────────────────────────
Year  1    847    801    1648   ← seeds + year-1 emergences
Year  2    531    147     678   ← declining fast
Year  3    211     61     272
Year  4     69     29      98
Year  5     46     10      56
Year  6     18      3      21
Year  7      9      2      11
Year  8      3      1       4
Year  9      2      1       3
Year 10      0      0       0   ← zero cases
```

**Elimination status:** Not achieved within the 10-year window — Year 10 is zero, but the criterion requires 3 consecutive zero years (Years 9 still has 3 cases). Elimination is imminent (Year 11–12 would likely confirm).

The 5-panel figure at `/tmp/guinea_worm_chad_model.png` shows incidence trends, district-level breakdown, intervention ramp-up curves, and 10-year cumulative burden.
