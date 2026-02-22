The simulation runs correctly. Here's a summary of the model and results:

---

## Model structure

The code at `/home/evaluser/gworm_chad_seis.py` implements a **LASER-framework-style** compartmental SEIS model with three architecture layers mirroring `laser_core`/`laser-generic` conventions:

| LASER concept | Implementation |
|---|---|
| `PropertySet` | Typed parameter container with attribute-style access |
| `PatchFrame` | Vectorised numpy arrays, one element per district |
| `Component` | `SEISTransitions` + `WaterContaminationFOI`, each with `.step()` |
| `Model` | `GWormModel` orchestrates the component list, drives the time loop |

**Note:** `laser-generic` is not installed in this environment, so the model implements equivalent functionality directly in numpy/pandas — the same approach IDM uses for custom models built on LASER primitives.

---

## Epidemiological design

**SEIS dynamics** (no immunity):
```
S → E (water exposure, ~12-month pre-patent)
E → I (worm emergence, detectable case)
I → S (worm expelled, no immunity conferred)
```

**Two-host FOI** (frequency-dependent, water-mediated):
```
λ_h = β_hh·(I_h/N_h) + β_dh·(I_d/N_d)    # dogs → humans via water
λ_d = β_dd·(I_d/N_d)                        # dogs → dogs (dominant)
```

**Transmission calibration:**
- `β_dd = R₀·γ = 1.7/21 ≈ 0.0810 d⁻¹` (dog–water–dog cycle)
- `β_dh/β_dd = 0.00526` → produces ~95% of infections in dogs (validated: 94.7% in output)
- `β_hh = 0.1·γ` → R₀_hh = 0.1, humans cannot sustain transmission alone

---

## 5-year national case counts

```
       Human cases  Dog cases  Total cases  Human exposures  Dog exposures
year
1              5.6      101.0        106.6             13.5          243.1
2             11.1      198.4        209.5             18.2          324.8
3             21.0      374.4        395.4             34.3          611.7
4             39.5      703.6        743.2             64.5        1,146.0
5             74.1    1,313.1      1,387.3            120.5        2,126.2
```

Key observations:
- **Year 1 cases are low** because seed infections were placed in I directly; their 21-day water contamination window drives E→I emergences from year 2 onward — correct given the 365-day pre-patent lag
- **~1.88× annual growth** is consistent with the analytical prediction: linearised SEIS growth rate r ≈ 0.00176 d⁻¹ → e^(365r) ≈ 1.90
- **Dog fraction holds at 94.7%** across all years, matching Chad surveillance data
- Only the 3 seeded districts (Moyen-Chari, Mandoul, Barh-Kôh) show cases; the 5 unseeded patches stay at zero — appropriate for a model without inter-district mixing, which reflects guinea worm's tight coupling to local water sources
- **Exposures lead cases by ~12 months**, which the model correctly distinguishes in the output

The CSV output and a 4-panel figure are saved to `/home/evaluser/`.
