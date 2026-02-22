The file `scripts/polio_seir_basic_10patch.py` already exists and implements exactly what you specified. Let me walk through the key design decisions so you can understand and run it.

**`scripts/polio_seir_basic_10patch.py`** — 394 lines, complete and runnable.

---

**Parameter derivation** (`polio_seir_basic_10patch.py:53-67`):

```python
"beta": 6.0 / 28.0,     # R0=6, D_inf=28d → beta = R0/D ≈ 0.2143
"exp_shape": 3,          # Latent period ~3 days (gamma, shape=3, scale=1)
"exp_scale": 1.0,
"inf_mean": 28,          # Infectious period ~28 days (normal)
"inf_sigma": 3,
```

The gamma latent distribution `Gamma(shape=3, scale=1)` has mean = 3 days. This is the correct parameterization for a 3-day latent period.

**Initial conditions** (`polio_seir_basic_10patch.py:111-115`):

```python
scenario["I"] = np.full(N_PATCHES, 5, dtype=np.uint32)      # 5 infectious
scenario["E"] = np.zeros(N_PATCHES, dtype=np.uint32)
scenario["R"] = np.full(N_PATCHES, int(0.95 * 100_000), dtype=np.uint32)  # 95% immune
scenario["S"] = (population - E - I - R)                    # remainder susceptible
```

Note: with 95% initial immunity and R0=6, the initial effective R = 6 × 0.05 = 0.3 < 1, so the disease would die out stochastically without importation. The model includes `PatchImportation` (seeding 3 endemic patches every 30 days) to sustain transmission while the susceptible pool rebuilds through births.

**Key architectural choices**:

| Concern | Choice | Why |
|---|---|---|
| Spatial coupling | Haversine gravity network, k=0.005 | Produces traveling waves |
| Seasonal forcing | Monsoon cosine profile (peak Sep, ±30%) | Pakistan WPV1 seasonality |
| Vital dynamics | CBR=29, CDR=7 /1000/yr | Pakistan demographics |
| Ages | Pyramid-initialized at t=0 | Prevents artificial RI pulse |
| Importation | `PatchImportation` (susceptible-targeted) | Prevents stochastic extinction |

**To run it:**

```bash
python scripts/polio_seir_basic_10patch.py
```

Output goes to `eval/outputs/polio_seir_basic_10patch.png` — a 4-panel diagnostic showing seasonal forcing, weekly incidence per patch, susceptible fraction (S/N vs 1/R0 threshold), and R_eff over time.

**One thing to be aware of**: the verification suite from `verification_checks.py` is not yet wired into this script. If you want automated health checks after the run, add this at the end of `run_model()`:

```python
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / ".claude/skills/laser-spatial-disease-modeling/scripts"))
from verification_checks import verify_model_health
verify_model_health(model)
```
