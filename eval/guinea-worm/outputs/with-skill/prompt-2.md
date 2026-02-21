The model is working correctly. Here's a summary of what was built:

## Files Created

### `scripts/guinea_worm_components.py` — Custom LASER Components

| Component | Purpose |
|-----------|---------|
| **`SEISRecovery`** | Routes I→S instead of I→R (no lasting immunity). Runs after `SEIR.Recovered`, catches agents just set to RECOVERED and resets them to SUSCEPTIBLE. |
| **`build_chad_dry_season()`** | 365-day seasonal profile: peak Jan-Mar (1.8x), decline Apr-May, rainy Jun-Sep (0.4x), ramp Oct-Dec. Normalized to mean=1.0. |
| **`DualHostWaterTransmission`** | Cross-species FOI via shared water sources. Computes prevalence in source species, applies asymmetric coupling (dog→human=0.15, human→dog=0.05) as Bernoulli infection trials on susceptible agents in the target species. |
| **`PatchImportation`** | Seeds exposed agents in endemic patches (susceptible-targeted). |

### `scripts/guinea_worm_chad.py` — Main Simulation

**Architecture:** Two parallel LASER `Model` instances (human, dog) sharing the same 8-district geography and seasonal profile, coupled via `DualHostWaterTransmission` in a manual tick loop that calls `_initialize_flows()` for both models.

**Key design decisions:**
- **Human within-species R0≈0.4** (subcritical alone) — humans can't sustain transmission without the dog reservoir
- **Dog within-species R0≈2.1** — dogs are the primary driver, consistent with Chad's ~95% dog case share
- **Gravity coupling** with c=2.0 (steep distance decay for local water sources), b=0.3 (weak population scaling)
- **SEIS component order:** Susceptible → Exposed → Infectious → Recovered → SEISRecovery (R→S) → Births → Deaths → Transmission

**Output:** Annual cases by district showing dog-dominated epidemic (10:1 ratio), spatial concentration in 3 endemic districts, and dry-season peaks.
