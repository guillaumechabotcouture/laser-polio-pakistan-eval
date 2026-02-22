The model runs correctly. All three constraints are verified:
- Row sums exactly 1.0, diagonal exactly 0.85, off-diagonal capped at 0.15
- Seasonal β ranging from 0.700× (Feb 13, trough) to 1.300× (Aug 15, peak)
- Epidemic propagates from the seed patch with realistic delay across the 450 km line

Since `/eval-work` is a read-only filesystem, the file lives at `/tmp/laser_polio_extended.py`. Here is the complete code:

---

```python
"""laser_polio_extended.py
========================
Extended LASER-style spatial SIR polio model.

Three extensions over a single-patch baseline:

  1. Gravity-model migration network (c=1.5, b=0.5)
       W[i, j] = N_j**b / d_ij**c   (j ≠ i)
       Rows = source patch, columns = destination patch.

  2. Seasonal transmission forcing
       β(t) = β₀ · [1 + A · cos(2π(DOY − DOY_peak) / 365)]
       A=0.3 → [0.7×, 1.3×] range
       Peak at DOY 227 (Aug 15, centre of Jul–Oct).
       Trough at DOY 44 (Feb 13, inside Dec–Mar).

  3. Row-normalised FOI export cap (15%)
       Σ_{j≠i} M[i,j] ≤ 0.15,  M[i,i] = 1 − export ≥ 0.85

Patches are equally spaced along a 1-D line, 50 km apart.
Each extension is a LASER-protocol component: __call__(model, step).
Dependencies: Python ≥ 3.9, NumPy only.
"""

from __future__ import annotations
import datetime
from dataclasses import dataclass, field
import numpy as np

try:
    from laser_core.model import Model as _LaserModel          # noqa
    from laser_core.laserframe import LaserFrame as _LF        # noqa
    _LASER_CORE_AVAILABLE = True
except ModuleNotFoundError:
    _LASER_CORE_AVAILABLE = False


# ===========================================================================
# Parameters
# ===========================================================================

@dataclass
class Params:
    """Central parameter store (mirrors laser_core.Parameters style)."""

    # Spatial layout
    n_patches: int = 10
    patch_spacing_km: float = 50.0

    # Demography
    initial_population: int = 100_000
    birth_rate: float = 1.0 / (70.0 * 365)
    death_rate: float = 1.0 / (70.0 * 365)

    # Epidemiology
    beta0: float = 0.35
    gamma: float = 1.0 / 14.0
    initial_prevalence: float = 1e-3        # seed patch only

    # Extension 1: gravity model
    gravity_b: float = 0.5
    gravity_c: float = 1.5

    # Extension 2: seasonal forcing
    # β(t) = β₀·[1 + 0.3·cos(2π(DOY−227)/365)]
    # → 1.3× at DOY 227 (Aug 15)   → 0.7× at DOY 44 (Feb 13)
    season_amplitude: float = 0.3
    season_peak_doy: int = 227

    # Extension 3: export cap
    max_export_fraction: float = 0.15

    # Simulation
    n_days: int = 365 * 5
    dt: float = 1.0
    start_date: datetime.date = field(
        default_factory=lambda: datetime.date(2020, 1, 1)
    )
    seed_patch: int = 0


# ===========================================================================
# Utilities
# ===========================================================================

def _distance_matrix(n: int, spacing_km: float) -> np.ndarray:
    """
    (n×n) inter-patch distances for a 1-D line of equally spaced patches.
    Diagonal = np.inf so 1/d_ii → 0 in gravity weights without special-casing.
    """
    coords = np.arange(n, dtype=float) * spacing_km
    D = np.abs(coords[:, None] - coords[None, :])
    np.fill_diagonal(D, np.inf)
    return D


def _day_of_year(start: datetime.date, step: int) -> int:
    return (start + datetime.timedelta(days=int(step))).timetuple().tm_yday


# ===========================================================================
# Extension 1 & 3 — GravityCouplingComponent
# ===========================================================================

def build_coupling_matrix(
    populations: np.ndarray,
    distances: np.ndarray,
    b: float = 0.5,
    c: float = 1.5,
    max_export: float = 0.15,
) -> np.ndarray:
    """
    Row-normalised gravity coupling matrix M.

    Convention: M[i, j] = fraction of patch i's infectious pressure
    that reaches patch j (rows=source, cols=destination).

    Construction
    ------------
    1. Raw off-diagonal weights:
           W[i, j] = N_j**b / d_ij**c     (diagonal = 0 because d_ii=inf)

    2. Scale each row so Σ_{j≠i} M[i,j] = min(raw_fraction, max_export).

    3. Set diagonal M[i,i] = 1 − Σ_{j≠i} M[i,j]   →  row sums = 1.

    FOI at patch j:
        λ_j = β(t) · (Mᵀ @ I/N)[j]   =  β(t) · Σ_i M[i,j]·(I_i/N_i)
    """
    n = len(populations)
    N_dest = populations.astype(float) ** b

    # Raw gravity weights (diagonal naturally 0: N/inf = 0)
    W = N_dest[None, :] / (distances ** c)          # (n, n)

    row_sums  = W.sum(axis=1)
    safe_sums = np.where(row_sums > 0, row_sums, 1.0)
    scale     = np.where(row_sums > 0, max_export / safe_sums, 0.0)
    M         = W * scale[:, None]                  # off-diagonal scaled

    off_diag_sums = M.sum(axis=1)                   # = max_export (or 0)
    np.fill_diagonal(M, 1.0 - off_diag_sums)        # local retention

    if not np.allclose(M.sum(axis=1), 1.0, atol=1e-12):
        raise RuntimeError("Coupling matrix rows do not sum to 1.")
    if not np.all(off_diag_sums <= max_export + 1e-9):
        raise RuntimeError("Export cap violated.")

    return M


class GravityCouplingComponent:
    """
    LASER component: builds model.M once at initialisation.
    Static for this model; recompute inside __call__ for full demography.
    """

    def __init__(self, model: "PatchModel") -> None:
        p = model.params
        model.M = build_coupling_matrix(
            populations=model.N,
            distances=model.distances,
            b=p.gravity_b,
            c=p.gravity_c,
            max_export=p.max_export_fraction,
        )

    def __call__(self, model: "PatchModel", step: int) -> None:
        pass   # static matrix; re-enable recomputation for dynamic demography


# ===========================================================================
# Extension 2 — SeasonalForcingComponent
# ===========================================================================

def seasonal_beta(
    beta0: float,
    doy: int | np.ndarray,
    peak_doy: int = 227,
    amplitude: float = 0.3,
) -> float | np.ndarray:
    """
    β(t) = β₀ · [1 + A · cos(2π(DOY − peak_doy) / 365)]

    With A=0.3:
      DOY 227 (Aug 15) → 1.30× β₀   (peak, centre of Jul–Oct window)
      DOY  44 (Feb 13) → 0.70× β₀   (trough, inside Dec–Mar window)
    """
    phase = 2.0 * np.pi * (np.asarray(doy, dtype=float) - peak_doy) / 365.0
    return beta0 * (1.0 + amplitude * np.cos(phase))


class SeasonalForcingComponent:
    """LASER component: updates model.beta each timestep."""

    def __init__(self, model: "PatchModel") -> None:
        self(model, 0)

    def __call__(self, model: "PatchModel", step: int) -> None:
        doy = _day_of_year(model.params.start_date, step)
        p   = model.params
        model.beta = float(
            seasonal_beta(p.beta0, doy, p.season_peak_doy, p.season_amplitude)
        )


# ===========================================================================
# SIR Dynamics Component
# ===========================================================================

class SIRDynamicsComponent:
    """
    LASER component: Euler step for the spatial SIR ODE.

    Reads model.M (gravity coupling) and model.beta (seasonal forcing)
    — both updated by preceding components in the pipeline.

    FOI:   λ_j = β(t) · (Mᵀ @ I/N)_j
    dS_i = births_i − μ·S_i − λ_i·S_i
    dI_i = λ_i·S_i − μ·I_i − γ·I_i
    dR_i = γ·I_i  − μ·R_i
    """

    def __call__(self, model: "PatchModel", step: int) -> None:
        p  = model.params
        dt = p.dt

        safe_N   = np.where(model.N > 0, model.N, 1.0)
        prev     = model.I / safe_N                  # infectious prevalence

        # Extensions 1 & 3: gravity-coupled effective prevalence
        eff_prev = model.M.T @ prev                  # (n,)

        # Extension 2: seasonally-adjusted FOI
        foi      = model.beta * eff_prev             # (n,)

        new_inf  = np.minimum(foi * model.S * dt, model.S)
        new_rec  = np.minimum(p.gamma * model.I * dt, model.I)

        births   = p.birth_rate * model.N * dt
        deaths_S = p.death_rate * model.S * dt
        deaths_I = p.death_rate * model.I * dt
        deaths_R = p.death_rate * model.R * dt

        model.S += births   - deaths_S - new_inf
        model.I += new_inf  - deaths_I - new_rec
        model.R += new_rec  - deaths_R

        model.S  = np.maximum(model.S, 0.0)
        model.I  = np.maximum(model.I, 0.0)
        model.R  = np.maximum(model.R, 0.0)
        model.N  = model.S + model.I + model.R


# ===========================================================================
# PatchModel — LASER-style model container
# ===========================================================================

class PatchModel:
    """
    Spatial SIR model container (mirrors laser_core.Model interface).

    Component pipeline per timestep:
        GravityCouplingComponent  →  model.M    (static after init)
        SeasonalForcingComponent  →  model.beta (updated each step)
        SIRDynamicsComponent      →  model.S/I/R/N
    """

    def __init__(self, params: Params) -> None:
        self.params = params
        n = params.n_patches

        self.N = np.full(n, float(params.initial_population))
        self.I = np.zeros(n)
        self.I[params.seed_patch] = (
            self.N[params.seed_patch] * params.initial_prevalence
        )
        self.R = np.zeros(n)
        self.S = self.N - self.I - self.R

        self.distances = _distance_matrix(n, params.patch_spacing_km)

        self._gravity  = GravityCouplingComponent(self)   # sets self.M
        self._seasonal = SeasonalForcingComponent(self)   # sets self.beta
        self._sir      = SIRDynamicsComponent()

        self.results: dict[str, np.ndarray] = {}

    def step(self, t: int) -> None:
        self._gravity(self, t)
        self._seasonal(self, t)
        self._sir(self, t)

    def run(self) -> dict[str, np.ndarray]:
        """
        Run for params.n_days steps.

        Returns dict with keys "t", "S", "I", "R", "beta".
        "S"/"I"/"R" are shape (n_days, n_patches).
        """
        n_days, n = self.params.n_days, self.params.n_patches
        t_hist    = np.empty(n_days, dtype=int)
        S_hist    = np.empty((n_days, n))
        I_hist    = np.empty((n_days, n))
        R_hist    = np.empty((n_days, n))
        beta_hist = np.empty(n_days)

        for t in range(n_days):
            self.step(t)
            t_hist[t] = t
            S_hist[t] = self.S
            I_hist[t] = self.I
            R_hist[t] = self.R
            beta_hist[t] = self.beta

        self.results = {
            "t": t_hist, "S": S_hist, "I": I_hist,
            "R": R_hist, "beta": beta_hist,
        }
        return self.results


# ===========================================================================
# Diagnostics
# ===========================================================================

def print_coupling_summary(model: PatchModel) -> None:
    M, p = model.M, model.params
    diag_vals     = np.diag(M)
    off_diag_sums = M.sum(axis=1) - diag_vals

    print("── Coupling matrix (Extensions 1 & 3) ─────────────────────────")
    print(f"  Gravity b={p.gravity_b}, c={p.gravity_c}, "
          f"max_export={p.max_export_fraction}")
    print(f"  Row sums        : min={M.sum(axis=1).min():.8f}, "
          f"max={M.sum(axis=1).max():.8f}")
    print(f"  Diagonal        : min={diag_vals.min():.6f}, "
          f"max={diag_vals.max():.6f}")
    print(f"  Off-diag export : min={off_diag_sums.min():.6f}, "
          f"max={off_diag_sums.max():.6f}")
    print(f"  Cap satisfied   : "
          f"{np.all(off_diag_sums <= p.max_export_fraction + 1e-9)}")
    print()


def print_seasonal_summary(params: Params) -> None:
    months = [
        ("Jan 15", 15),  ("Feb 13", 44),  ("Mar 15", 74),
        ("Apr 15", 105), ("May 15", 135), ("Jun 15", 166),
        ("Jul 15", 196), ("Aug 15", 227), ("Sep 15", 258),
        ("Oct 15", 288), ("Nov 15", 319), ("Dec 15", 349),
    ]
    print("── Seasonal forcing (Extension 2) ──────────────────────────────")
    print(f"  β₀={params.beta0}, A={params.season_amplitude}, "
          f"peak DOY={params.season_peak_doy}")
    for label, doy in months:
        b    = seasonal_beta(params.beta0, doy, params.season_peak_doy,
                             params.season_amplitude)
        mult = b / params.beta0
        bar  = "█" * int(round(mult * 15))
        print(f"  {label:<10} DOY={doy:3d}  β={b:.4f}  {mult:.3f}×  {bar}")
    print()


def print_epidemic_summary(model: PatchModel) -> None:
    I_ts, S_ts = model.results["I"], model.results["S"]
    T, n = I_ts.shape
    N0   = model.params.initial_population
    total_I = I_ts.sum(axis=1)
    peak_t  = int(total_I.argmax())

    print("── Epidemic summary ─────────────────────────────────────────────")
    print(f"  Patches {n}, Days {T}, Peak I={total_I[peak_t]:,.0f} (day {peak_t})")
    print(f"\n  {'Patch':>5}  {'Peak I':>10}  {'Peak day':>9}  {'Attack rate':>12}")
    for k in range(n):
        print(f"  {k:>5}  {I_ts[:,k].max():>10,.0f}  "
              f"{int(I_ts[:,k].argmax()):>9}  "
              f"{1 - S_ts[-1,k]/N0:>11.1%}")
    print()


# ===========================================================================
# Entry point
# ===========================================================================

def main() -> None:
    params = Params(n_patches=10, n_days=365 * 5)

    print("=" * 65)
    print("Extended LASER Polio Model")
    print(f"  laser_core available : {_LASER_CORE_AVAILABLE}")
    print(f"  Patches / spacing    : {params.n_patches} × {params.patch_spacing_km} km")
    print(f"  Population / patch   : {params.initial_population:,}")
    print(f"  Simulation           : {params.n_days} days ({params.n_days//365} yr)")
    print("=" * 65 + "\n")

    model = PatchModel(params)
    print_coupling_summary(model)
    print_seasonal_summary(params)

    print("Running simulation ...", flush=True)
    model.run()
    print("Done.\n")
    print_epidemic_summary(model)


if __name__ == "__main__":
    main()
```

---

## Design notes

### Extension 1 — Gravity coupling

The raw off-diagonal weight from patch `i` to patch `j` is:

```
W[i,j] = N_j**b / d_ij**c     b=0.5, c=1.5
```

Setting `distances[i,i] = np.inf` makes the diagonal zero automatically (`N/inf = 0`), so no special-casing is needed. Only the destination population exponent `b` appears here; origin-population scaling (`N_i`) is already captured implicitly by the `I_i/N_i` term in the FOI.

### Extension 2 — Seasonal forcing

A single cosine centred on **DOY 227 (Aug 15)** — the midpoint of the Jul–Oct window — gives:

| Period | Multiplier |
|--------|-----------|
| Aug 15 (DOY 227) | **1.30×** |
| Jul 15 / Sep 15 | 1.26× |
| Oct 15 | 1.15× |
| Feb 13 (DOY 44, trough) | **0.70×** |
| Jan 15 / Mar 15 | 0.74× |

The trough falls at DOY 44 (Feb 13), firmly inside the Dec–Mar window.

### Extension 3 — Row normalisation

Each row of `M` represents one source patch's total infectious-pressure budget. After scaling the off-diagonal weights to sum to `max_export = 0.15`, the diagonal is set to `1 − 0.15 = 0.85`. The matrix is row-stochastic, so:

```
FOI at patch j:  λ_j = β(t) · (Mᵀ @ I/N)[j]
```

The `Mᵀ` transpose collects contributions *into* each destination across all sources.

### LASER compatibility

Each extension is a callable class with `__call__(model, step)`. Dropping them into a `laser_core.Model` pipeline requires only adding them to the component list — the interface is identical to laser_core's Component protocol.
