The model now produces realistic endemic dynamics. Here's a summary of the results:

**Model structure:**
- Stochastic SEIR with tau-leap (Poisson/binomial draws) at daily resolution
- 20 districts with heterogeneous populations (60k–2M) and per-district routine coverage (25%–92%, pop-weighted avg 80%)
- Gravity-model spatial coupling (Haversine distances, row-normalized)
- Monsoon cosine forcing peaking late August (±30%)
- SIA pulse vaccination every 180 days targeting under-5 susceptibles
- 3 imported infections every 60 days for years 1–5

**Key results (post burn-in, years 11–20):**

| Pattern | Districts | Zero-incidence weeks |
|---|---|---|
| Endemic (persistent) | N. Waziristan, S. Waziristan, Tank | 0–14% |
| Frequent outbreaks | Bannu, D.I. Khan, Khyber, Lakki Marwat | 12–32% |
| Sporadic sparks | Quetta, Multan, Zhob, Mohmand | 33–55% |
| Near-elimination | Karachi, Lahore, Islamabad, Rawalpindi | 60–99% |

This mirrors real Pakistan polio epidemiology: persistent circulation in low-coverage tribal/conflict areas (Waziristan corridor) with sporadic seeding into higher-coverage urban areas via spatial coupling. The code is a single runnable file at `polio_laser_model.py`.
