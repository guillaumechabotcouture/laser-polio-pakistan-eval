Since the workspace filesystem is read-only the code was presented inline above. Here's a precise summary of what's there and what's different from the two existing files.

---

## What's in the existing files vs. what the framework above provides

| | `calibrate_polio.py` | `polio_calibration.py` (old) | **Framework above** |
|---|---|---|---|
| AFP zero-week probability | `(199/200)^W` ✓ | raw `infections==0` ✗ | `(199/200)^W` ✓ |
| Logistic CCS in **scoring** | plotting only ✗ | shape penalty only ✗ | curve distance ✓ |
| Parameter sampler | LHS ✓ | LHS ✓ | LHS ✓ |
| Seasonal profile | `max(1+A*cos, 0)` | `1+A*cos` | `1+(amp-1)*cos`, mean=1 guaranteed |
| Borda-count ranker | ✗ | ✗ | ✓ |

---

## The three non-obvious choices

**`(199/200)^W` instead of binomial sampling** — a district with 50 infections/week has P(AFP=0) ≈ 0.78. Calling that a "zero week" for simulation purposes is correct; calling it a non-zero week (by thresholding on raw infections) would systematically inflate simulated incidence and make the model appear to need a lower beta to match observed zero-week proportions.

**Logistic curve distance as primary score** — direct SSE on per-district `prop_zero` punishes raw deviation uniformly. The CCS logistic tests whether the *shape* of the fadeout gradient (large city = continuous, small district = intermittent) is reproduced, which is the epidemiologically meaningful target.

**`amplitude - 1.0` as the cosine coefficient** — with the user's range U(1.0, 1.5), `amplitude=1.0` gives a coefficient of 0 (flat profile), `amplitude=1.5` gives a peak of 1.5× and trough of 0.5×. The mean is exactly 1.0 by construction (cosine integrates to zero), so no normalization call is needed and the unit check is a simple assertion.
