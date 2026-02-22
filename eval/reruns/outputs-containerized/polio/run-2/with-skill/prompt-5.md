All 5 health checks pass. Here's a summary of what was built and the results.

---

## Model: `/tmp/pakistan_polio_seir.py`

Run with:
```bash
python3 /tmp/pakistan_polio_seir.py --outdir /tmp/polio_out
```

### Architecture

**Custom components (all inline, no external dependencies):**

| Component | Notes |
|---|---|
| `PerPatchVaccination` | RI (weekly at age 270d) + SIA (every 180d, 0–5yr). Uses `add_scalar_property` — **not** `add_property`, which doesn't exist in LASER's `LaserFrame` and would silently crash the existing workspace file. Sets `state=RECOVERED` (not `susceptibility=0`) so `TransmissionSE` correctly treats vaccinated agents as non-susceptible. |
| `PatchImportation` | Seeds susceptible agents only in Quetta / Peshawar / N_Waziristan (1 per patch × 3 patches = 3 infections per event), every 60 days, first 5 years. |

**Built-in LASER components:** `SEIR.Transmission` with `seasonality=ValuesMap.from_timeseries(...)`, `BirthsByCBR` (29‰/yr), `MortalityByCDR` (7‰/yr).

### Verified outputs

**Health checks:** All 5 passed — population grew 1.54× over 20 years, no negative compartments, network active (sum=0.10), 100% of weeks had incidence, susceptible fraction declined from 15.7% → 14.2% confirming vaccination works.

**Epidemiological pattern (post-burn-in):**

| District | Province | Total Inc (10yr) | % Zero Weeks |
|---|---|---|---|
| Quetta | Balochistan | 64,558 | 0% (endemic) |
| DI Khan | KP | 73,384 | 5% (endemic) |
| Bannu | KP | 41,983 | 15% |
| N_Waziristan | KP | 41,233 | 27% |
| Peshawar | KP | 3,775 | 20% |
| Lahore | Punjab | 260 | 73% (sporadic) |
| Islamabad | ICT | 26 | 96% (near-elimination) |
| Jacobabad | Sindh | 4 | 99% |

The gradient is epidemiologically correct: KP tribal / Balochistan border districts sustain endemic transmission (high `unreachable_frac` > 1/R₀ ≈ 0.17), while Punjab / ICT see near-elimination under 80% RI + SIA. The zero-incidence weeks CSV and a 6-panel diagnostic PNG are written alongside the weekly incidence CSV.
