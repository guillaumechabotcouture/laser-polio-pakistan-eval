The script is written at `/home/evaluser/polio_seir_pakistan.py`. Here's a breakdown of every design decision:

---

## Model design

**Transmission rate**  
β = R₀ / T_infectious = 6 / 28 ≈ **0.2143 day⁻¹**. LASER uses a daily Bernoulli FOI: p = 1 − exp(−β · I/N), plus gravity-coupled cross-patch transmission.

**Duration distributions** (`laser.core.distributions`)  
| Parameter | Distribution | Mean | SD |
|-----------|-------------|------|----|
| Latent (σ) | `Gamma(shape=3, scale=1)` | 3 d | 1.73 d |
| Infectious (γ) | `Gamma(shape=4, scale=7)` | 28 d | 14 d |

Gamma enforces positivity and is standard for sojourn times. `dists.gamma(...)` returns a callable `f(tick, node) → float` consumed by the components.

**Initial conditions per patch (100,000 pop)**  
- R = 95,000 (95 % immune)  
- I = 5 (infectious seed)  
- E = 0  
- S = 4,995  

**Effective R at t=0**: Re = 6 × 0.04995 ≈ **0.30 < 1** — the model will show endemic persistence at very low incidence rather than explosive transmission, which is the expected behavior for a near-eradication polio scenario.

**Spatial structure**: 10 patches in a 2×5 grid of ~10 km square polygons placed over Pakistan's central Punjab/KPK region. Gravity coupling (`k=500, a=1, b=1, c=2`) allows cross-district diffusion; these defaults can be calibrated to Pakistan mobility data.

**Component pipeline order**:  
`Susceptible → Exposed → Infectious → Recovered → Transmission`

**Outputs** (`model.nodes.{S,E,I,R}`, shape `[nticks+1, n_patches]`):  
- Aggregate SEIR time-series  
- Peak infectious day, net attack rate, final immunity coverage  
- 4-panel figure saved to `polio_seir_simulation.png`
