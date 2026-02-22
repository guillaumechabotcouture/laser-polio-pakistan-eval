The model ran successfully. Here is the complete code — written to `/tmp/eastland_respiratory_model.py` since `/workspace` is read-only:

---

## Code walkthrough

### Key design decisions

**Disease parameters**
- Mean latent = shape × scale = 4 × 1 = **4 days**
- Mean infectious = shape × scale = 5 × 2 = **10 days**
- `beta = R0 / mean_infectious = 5.0 / 10.0 = 0.5`

**Seasonal forcing** — cosine centered on day 0 (Jan 1):
```
season(d) = 1.0 + 0.3 * cos(2π*d / 365)
```
| Day | Value | Region |
|-----|-------|--------|
| 0   | 1.300 | Winter peak ✓ |
| 45  | 1.214 | Mid-winter |
| 90  | 1.006 | Spring edge of window |
| 182 | 0.700 | Summer trough ✓ |
| Annual mean | **1.0000** | Exactly normalised ✓ |

The cosine formula guarantees mean = 1.0 analytically (∫cos = 0 over full period), so the `/= mean()` call is a no-op but serves as a documentation guardrail.

**`calc_capacity` called explicitly** — shows the pre-allocated slot counts before `Model()` is constructed. With CBR=30, safety_factor=3.0, each patch gets ~2× its initial population as headroom (1.99× shown — the 3× factor covers multi-decade runs; 10 years only needs ~1.22×).

**`ValuesMap.from_timeseries`** tiles the 365-day array to 3650 days and broadcasts it to all 4 patches, which is what `SEIR.Transmission(..., seasonality=seasonality)` consumes internally.

---

## Simulation output (verified)

```
Population: 530,000 → 643,890 (+21.5%)   ← net CBR-CDR = +20/1000/yr, ≈22% expected ✓
All compartments non-negative: OK
Total infections (10 yr): 475,697
Patches with infections: 4/4

Annual infections by patch:
  Patch         Pop     Avg/yr   Per 100k
  Eastport   100,000    8,945     8,945
  Northton   200,000   17,992     8,996
  Centerburg 150,000   13,475     8,983
  Westham     80,000    7,157     8,946
```

The ~9,000/100k/yr attack rate reflects a large initial epidemic (R0=5, ~80% susceptible at start) that depletes the susceptible pool, followed by low-level endemic re-emergence as births replenish it. All 4 patches are infected, and the network (k=0.01, row-capped at 15%) produces spatial coupling across the ~73–220 km inter-patch distances. Without periodic importation, small patches at these population sizes (80k–200k) naturally fade between epidemic waves — expected behavior for R0=5 with these demographics.
