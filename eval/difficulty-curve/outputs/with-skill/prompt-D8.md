The simulation ran successfully. Here's the summary:

**Model configuration** (all specs met):
- 4 patches: Eastport (100k), Northton (200k), Centerburg (150k), Westham (80k)
- Disease: R0=5, beta=0.5, latent gamma(4,1), infectious gamma(5,2)
- Initial: 90% S, 0% E, 1% I, 9% R
- Gravity network (k=0.01, a=1, b=1, c=1.5) with `row_normalizer(0.15)`
- Seasonal forcing: winter peak 1.3x, summer trough 0.7x via `ValuesMap`
- Custom `PulseImportation`: 3 infections into patch 0 every 60 ticks
- `BirthsByCBR` (CBR=30) and `MortalityByCDR` (CDR=10), `calc_capacity` for 10 years

**Verification results:**
- **Population growing**: 530k → 643k (+21.4%) — confirms CBR/CDR units correct
- **All compartments non-negative**: OK
- **All 4 patches have infections**: spatial coupling working
- **Total infections**: 480,602 over 10 years (~48k/yr)
- **Per-patch incidence**: ~90k per million/yr (proportional across patches)

**Diagnostics** (6-panel plot):
- **(A)** Seasonal profile: 1.3 winter peak, 0.7 summer trough
- **(B)** Initial explosive epidemic (90% susceptible + R0=5), then settles to endemic equilibrium with recurring outbreaks
- **(C)** Susceptible fraction drops from 0.9 to ~0.2 (1/R0), confirming correct endemic equilibrium
- **(D)** Population grows steadily (CBR > CDR)
- **(E)** R_eff drops from ~4.5 to ~1.0, oscillating around 1 in equilibrium
- **(F)** Heatmap shows initial epidemic flash then sustained low-level transmission

The complete script is at `eastland_respiratory_model.py`.
