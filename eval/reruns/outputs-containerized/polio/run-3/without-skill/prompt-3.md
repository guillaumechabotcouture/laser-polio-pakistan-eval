---

## Design rationale

### Compartments and why five are needed

| Compartment | Immunity source | Susceptibility multiplier | Tracked separately? |
|---|---|---|---|
| `S` | None | 1.0 | — |
| `E` | — | — | — |
| `I` | — | — | — |
| `R_nat` | Natural infection | 0 (full) | **Yes** |
| `V_opv` | OPV (mucosal IgA) | `1 − opv_effectiveness` = 0.15 | **Yes** |

`V_opv` and `R_nat` are distinct because an OPV recipient can still be colonised, shed virus, and — on breakthrough infection — transitions through `E → I → R_nat`, gaining the stronger natural immunity. Merging them would underestimate residual transmission risk in high-SIA settings.

### Component mapping to LASER built-ins

| Component | LASER analogue | Custom or adapted? |
|---|---|---|
| `InitializeEquilibriumStates` | `InitializeEquilibriumStatesProcess` | Adapted — solves for `S/V/R` split accounting for partial OPV susceptibility |
| `VitalDynamicsWithOPV` | `BaseVitalDynamicsProcess` | Adapted — replaces direct-to-R birth routing with a **42-slot circular birth queue** for exact 6-week timing |
| `PolioTransmission` | `InfectionProcess` (SEIR) | Adapted — extends to SVEIR; adds `V_opv → E` breakthrough pathway with `λᵛ = λ(1 − ε)` |
| `SIACampaigns` | `SIACalendarProcess` | **Closely mirrors** the laser-measles pattern; replace the internal schedule with a WHO Polars DataFrame to use it directly |
| `OutcomeTracker` | `StateTracker` | Custom — tracks all five compartments including the `V_opv`/`R_nat` split |

### Key implementation details

**6-week RI delay (birth queue):** Slot `tick % 42` is **read first** (vaccinating the 42-day-old cohort) then **overwritten** with today's births. This is the minimal correct implementation of a delayed differential equation in a daily-step model without individual agent tracking.

**OPV vs natural immunity in the FOI:** `PolioTransmission` applies the complementary-log-log transformation (`1 − exp(−λ)`) to both `S` and `V_opv`, but uses `λᵛ = λ(1−ε)` for the latter. This exactly captures the mucosal IgA barrier without any ad-hoc approximation.

**SIA age targeting:** The compartmental approximation `S_under5 = S × under5_frac` is documented with a clear upgrade path to PHASER's explicit age-cohort tracking for production calibration work.
