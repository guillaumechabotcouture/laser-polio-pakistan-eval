The model is complete and all verification checks pass. Here's a summary of the results:

## Model Output

The complete script is at `scripts/guinea_worm_chad_interventions.py`. It builds two coupled LASER agent-based models (human + dog) with all the specified features.

### Key Results

**Disease trajectory** — Infections decline sharply as interventions scale up:
- Year 1: 393 human cases, 617 dog infections
- Year 5: 111 human, 277 dog (interventions at 70-82%)
- Year 10: 1 human, 3 dog (near-elimination)

**Elimination** was not quite achieved within 10 years (need 3 consecutive zero-case years), but transmission dropped to near-zero by year 10. The model shows the characteristic pattern: dogs maintain transmission longer than humans due to higher R0 and lower intervention coverage.

**Spatial pattern** — Infections concentrate in the 3 endemic seeded districts (Chari-Baguirmi, Moyen-Chari, Salamat), with only minor spillover to Mayo-Kebbi Est and Mandoul via the gravity network.

### Architecture

- **Two parallel LASER Models** (2.2M human agents, 180K dog agents) across 8 districts
- **SEIS dynamics** via `SEISRecovery` component (R→S each tick, no lasting immunity)
- **IntervenedWaterTransmission** replaces standard SEIR.Transmission — gravity-coupled FOI with ABATE, filter, containment/tethering modifiers
- **IntervenedDualHostTransmission** handles cross-species coupling (dog→human 0.15, human→dog 0.05)
- **Time-varying interventions** updated yearly at tick boundaries
- All 10 verification checks pass (SEIS R≈0, non-negativity, infections > 0, population growth, network connectivity)
