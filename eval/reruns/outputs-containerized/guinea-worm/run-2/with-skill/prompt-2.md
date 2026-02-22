The model is at `/tmp/guinea_worm_water_model.py` and all verification checks pass. Here's what was implemented and why each choice was made:

---

## What was built

### Two new transmission components

**`WaterSourceTransmission`** (within-species, replaces `SEIR.Transmission`)

Makes the copepod-mediated biology explicit. Each tick:
1. Computes per-patch water contamination as `I(t) / N(t)` (infectious hosts releasing larvae)
2. Applies gravity-coupled sharing: `coupled_prev[i] = (1 - row_sum[i]) × prev[i] + network[i] @ prev`
3. Multiplies by the dry-season copepod concentration multiplier
4. Runs Bernoulli trials on susceptible agents

The loop `for i in range(nnodes): coupled_prev[i] = local_frac * prev[i] + network[i] @ prev` is identical for both within- and cross-species paths — deliberate, since both flow through the same physical water bodies.

**`CrossSpeciesWaterTransmission`** (cross-species — key extension)

The critical difference from the existing `DualHostWaterTransmission`: the source species' contamination is **gravity-coupled before being applied to the target species**. In the base model, dogs in patch `i` only affect humans in patch `i`. Here, dogs in patch `i` contaminate water shared with neighboring patches, so humans there face elevated FOI — matching how seasonal ponds and rivers straddle district boundaries.

### Design decisions

| Decision | Rationale |
|---|---|
| Same gravity network for within- and cross-species | Both processes flow through the same physical water bodies |
| `WaterSourceTransmission` before `CrossSpeciesWaterTransmission` | Already-infected agents are excluded from cross-species trials via `people.state` check |
| `SEISRecovery` imported from existing components | No duplication; SEIS mechanic unchanged |
| `DOG_TO_HUMAN_COUPLING = 0.15`, `HUMAN_TO_DOG_COUPLING = 0.05` | Preserves calibrated asymmetry from the existing model; 3:1 ratio reflects dogs' dominant contamination role; reaching the empirical 50:1 case ratio requires fitting beta against Carter Center data |
| `season_tiled` (raw array) instead of `ValuesMap` | Simpler direct indexing; `ValuesMap` only needed by `SEIR.Transmission` |
