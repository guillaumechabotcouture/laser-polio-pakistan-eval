The workspace is mounted read-only so the files can't be edited directly — the changes are shown above and correspond exactly to what the upgrade requires.

**Summary of what changes and why:**

- **`additional_states=["V"]` on `Model()`** is the most important line. Without it, LASER doesn't propagate `V[t]→V[t+1]` each tick, so `VaccinatedCompartment` would compute `0 - waned` and go negative immediately.
- **`PerPatchVaccinationSEIRV` replaces `PerPatchVaccination`** because it sets `state = VACCINATED` (4) instead of `RECOVERED` (3), enabling the model to later transition `V → S` via the waning timer. The built-in `RoutineImmunizationEx` is not usable here because it also sets `state = RECOVERED`, conflating OPV with natural immunity.
- **SIA coverage changes from 80% → 90%** per the specification; RI stays at 80%.
- **`VaccinatedCompartment` placed before `BirthsByCBR`** so waning is processed before vaccination in the same tick — agents who wane become susceptible and could (correctly) be re-vaccinated by an SIA on the same day.
- **All population totals updated to `S+E+I+R+V`** in `print_summary` and `plot_diagnostics` — failing to include V would undercount population and overstate S/N (R_eff would appear higher than it is).
