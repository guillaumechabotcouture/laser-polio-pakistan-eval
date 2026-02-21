The model runs correctly and produces a 6-panel diagnostic plot. Here's a summary of what was added:

---

### Changes to `laser_polio_model.py`

**New compartment: V (Vaccinated/OPV-immune)** — the model is now SEIRV. V is tracked separately from R because OPV provides mucosal immunity that wanes (~3 yr), while natural infection (R) confers stronger, permanent immunity.

**New parameters** (lines 87-92):
| Parameter | Value | Meaning |
|---|---|---|
| `ri_coverage` | 0.80 | RI coverage at 6 weeks (OPV) |
| `sia_coverage` | 0.90 | SIA coverage for children 0-5 |
| `sia_interval` | 182 | SIA every 6 months |
| `fraction_under5` | 5/70 | Demographic proxy for age targeting |
| `ve_opv` | 0.50 | Leaky vaccine efficacy (50% FOI reduction) |
| `waning_rate` | 1/(3\*365) | V→S waning (~3 year duration) |

**Three new LASER-patterned phase components:**

1. **`RoutineImmunizationPhase`** — diverts 80% of daily births (mu\*N) from S→V, approximating OPV at 6 weeks in the aggregate model. Analogous to LASER's `Immunization.RoutineImmunization`.

2. **`SIAPhase`** — every 182 days, moves 90% of estimated susceptible children 0-5 (fraction\_under5 \* S) from S→V. Analogous to LASER's `Immunization.Campaign`. The 6 SIA events are visible as orange vertical lines in panel (c).

3. **`WaningPhase`** — V→S at rate 1/(3\*365) per day. Only OPV immunity wanes; natural immunity (R) is permanent.

**Modified components:**

- **`TransmissionPhase`** — now handles both S→E (full FOI) and V→E (FOI \* 0.5, leaky vaccine). Breakthrough infections in V move through E→I→R, so re-infected vaccinees gain natural immunity.
- **`DemographyPhase`** — deaths now drawn from V as well.
- **`RecorderPhase`** — snapshots V alongside S/E/I/R.

**Phase execution order:** Transmission → Latent → Recovery → Demography → RI → SIA → Waning → Record

**Results:** The model ran for 3 years across 10 patches, conducted 6 SIA campaigns, and ended with ~64,400 in V (OPV-immune) vs ~2.3M in R (naturally immune), with ~65,800 remaining susceptible.
