The model runs successfully. Here is the complete code (written to `/home/evaluser/eastland_respiratory_model.py`):

---

**Key design decisions and outputs:**

**Beta derivation:**
- `D_inf = shape × scale = 5 × 2 = 10 days`
- `beta = R0 / D_inf = 5 / 10 = 0.5 day⁻¹`

**calc_capacity:**  Called explicitly before `Model()` for visibility. `Model.__init__` also calls it internally when `birthrates` is passed. With CBR=30, CDR=10 → ~2%/yr net growth → ~22% over 10 years. Safety factor 2.0 gives capacity of ~177k–354k per patch vs. initial 80k–200k.

**Gravity network (after `row_normalizer(0.15)`):**
```
Eastport   → Northton 11.2%, Centerburg 3.0%, Westham 0.9%  (row_sum = 15%)
Northton   → Eastport  5.4%, Centerburg 8.1%, Westham 1.5%  (row_sum = 15%)
Centerburg → Eastport  1.7%, Northton   9.5%, Westham 3.8%  (row_sum = 15%)
Westham    → Eastport  1.2%, Northton   4.4%, Centerburg 9.4% (row_sum = 15%)
```
All patches hit exactly the 15% cap — expected, since with k=0.01 and a=1 the raw gravity flows are large enough that `row_normalizer` clips all of them.

**Simulation results (10 years):**
| Check | Result |
|---|---|
| Population | 530k → 644k (+21.5%), PASS |
| Compartments ≥ 0 | PASS |
| Total infections | 470,569 |
| Patches infected | 4/4 |
| Annual rate | ~8,800–8,900 per 100k/yr |

**Component order:**
```
Susceptible → Exposed → Infectious → Recovered → Transmission → BirthsByCBR → MortalityByCDR
```
`Susceptible` and `Recovered` bracket the transition steps to preserve `S+E+I+R = N` each tick. `Transmission` (= `TransmissionSE`) reads `model.network` to compute spatially-coupled force of infection.
