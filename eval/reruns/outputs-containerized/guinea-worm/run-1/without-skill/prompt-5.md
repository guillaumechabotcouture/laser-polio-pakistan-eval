The model runs cleanly. Here's a summary of what was built and what the simulation found:

---

## Model: `/home/evaluser/guinea_worm_chad.py`

**Architecture** — LASER-inspired SEIS compartmental model implemented from scratch in pure numpy/pandas/matplotlib (LASER framework is not installed; scipy is also absent). Uses a classical 4th-order Runge–Kutta integrator at dt = 1 day.

### Key design decisions

| Component | Implementation |
|---|---|
| Compartments | S, E, I per species per district (48 state variables) |
| Transmission | Water-source mediated: `FOI = β × water_contamination × season × (1−ABATE) × (1−filter)` |
| Cross-species | `water_h = 0.05·prev_h + 0.95·prev_d`; dogs dominate reservoir |
| Seasonal forcing | `f(doy) = 1.1 + 0.7·cos(2π·doy/365)` — peaks at 1.8× (Jan 1), troughs at 0.4× (Jul 1) |
| Spatial coupling | Gravity model W[i,j] ∝ pop_i·pop_j/dist_ij², normalised so 15% of FOI crosses districts |
| Interventions | ABATE and filters reduce FOI; containment/tethering reduce source contamination |

### Simulation results (10 years)

- **Year 1 total burden**: ~3,974 cases (3,303 human + 671 dog)
- **Year 10 total burden**: ~4.8 cases — **99.9% reduction**
- **Elimination status**: Not confirmed within 10 years (requires < 0.5 cases for 3 consecutive years), but Year 10 burden is 4.8 — **elimination is imminent** just beyond the horizon

The dynamics show epidemiologically coherent behaviour: human effective R0 drops below 1 immediately in Year 1 (0.63) due to filters + ABATE; dogs remain above R0 = 1 until Year 5 (dog chain requires higher ABATE coverage to suppress), explaining why dogs sustain transmission longer. The three seeded districts (Chari-Baguirmi, Moyen-Chari, Salamat) contribute ~80% of burden throughout. The summary figure is saved to `/home/evaluser/guinea_worm_chad.png`.
