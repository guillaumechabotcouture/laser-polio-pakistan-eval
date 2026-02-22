The model is written to `/home/evaluser/chad_guinea_worm_model.py` (workspace root is read-only in this eval environment).

---

Here's what the file contains and the key design decisions:

## Architecture

**Two parallel LASER Models** (human + dog) stepped manually each tick — `model.run()` cannot coordinate two simultaneous models, so the tick loop is explicit.

```
tick loop:
  1. _initialize_flows(tick)  ← propagate counts t → t+1
  2. human_model components   ← SEIR state + births/deaths + interventions + within-species FOI
  3. dog_model components     ← same
  4. cross_host.step(tick)    ← dog↔human water bridge
  5. importation              ← 1 traveler every 90 days, years 1–3
  annual update (every 365 ticks): ramp all intervention coverages
```

## Component breakdown

| Component | What it does |
|-----------|-------------|
| `SEISRecovery` | Intercepts all R agents each tick and moves them back to S — implements zero immunity |
| `ABATELarvicide` | Writes `nodes.abate_factor` each tick; active Jan–May, kills 80% of copepods per coverage unit |
| `WaterFilterDistribution` | Assigns `has_filter` per agent at init + birth; passive flag read by transmission |
| `CaseContainment` | Assigns per-infectious-human containment status (detected 70%→97%); reduces their water contamination by 90% |
| `DogTethering` | Same pattern for dogs (40%→65% coverage); reduces contamination by 85% |
| `IntervenedWaterTransmission` | Full within-species FOI: effective prevalence → gravity-coupled → seasonal → ABATE → Bernoulli → filter protection |
| `IntervenedDualHostTransmission` | Cross-species: dog effective prevalence → human infections (Path A, dominant); human → dog (Path B, weak) |

## Key epidemiological choices

- **beta = R0 / INF_MEAN**: `HUMAN_BETA = 1.5/21 ≈ 0.071`, `DOG_BETA = 2.5/21 ≈ 0.119` — per-day transmission rate in a fully susceptible population
- **Cross-species coupling** (`DOG_TO_HUMAN=0.15`, `HUMAN_TO_DOG=0.05`): dogs dominate water contamination; ~95% of Chad infections are canine
- **Gravity network** with high distance decay (`c=2.0`): water source access is highly local; caps at 15% export per district
- **Importation ends at tick 1095** (3 × 365) exactly per spec

Run with:
```bash
python /home/evaluser/chad_guinea_worm_model.py
```
