The simulation ran successfully. Here's a summary of the model and results:

**Model architecture — SEIS (no lasting immunity):**
- `Susceptible` → `TransmissionSE` (S→E) → `Exposed` (E→I) → `InfectiousIS` (I→S)
- Pre-patent period: gamma(shape=40, scale=9) ≈ 360 days mean
- Infectious period: gamma(shape=7, scale=3) ≈ 21 days mean
- Dry-season seasonality (±50% sinusoidal, peak March)
- Gravity-based spatial coupling with steep distance decay (c=3)

**Two separate LASER models** — one for humans (2M, beta=0.055, R0≈1.16) and one for dogs (200K, beta=0.09, R0≈1.89), reflecting dogs as the dominant reservoir.

**Results — Annual case counts (worm emergences):**

| Year | Human | Dog | Total | Dog share |
|------|-------|-----|-------|-----------|
| 1 | 1 | 62 | 63 | 98.4% |
| 2 | 3 | 172 | 175 | 98.3% |
| 3 | 2 | 353 | 355 | 99.4% |
| 4 | 3 | 665 | 668 | 99.6% |
| 5 | 2 | 1,405 | 1,407 | 99.9% |

The results are epidemiologically consistent:
- **Humans**: barely sustained transmission (R0≈1.16), ~1-3 cases/year — consistent with near-elimination
- **Dogs**: exponential growth from the 50-dog seed (R0≈1.89) — shows what happens without interventions, confirming dogs as the reservoir driving persistence
- **Dog share ~98-99%** matches Chad's real-world pattern (~95% of detected infections in dogs)
- **Spatial spread**: cases start in the 3 seeded patches (Moundou, Sarh, Am Timan) and gradually spill into neighboring districts
