The working directory is mounted read-only (`ro,fakeowner`), so I'll present the full framework directly. Here is the complete calibration framework, designed around LASER's `PropertySet`/`LaserFrame` conventions, structured as a `calibration/` package.

---

## `calibration/__init__.py`

```python
"""LASER polio calibration framework.

Calibrates beta, gravity_k, and seasonal_amplitude against AFP
surveillance data using zero-week fraction as the goodness-of-fit metric.
"""

from .loop import CalibResult, CalibrationStudy
from .params import CalibParams, ParameterSampler
from .runner import SimOutput, run_single_trial
from .scoring import (
    PARALYSIS_TO_INFECTION_RATIO,
    compute_afp_from_infections,
    daily_to_weekly,
    district_weights_by_burden,
    score_zero_weeks,
    zero_week_fraction,
)

__all__ = [
    "CalibParams",
    "CalibResult",
    "CalibrationStudy",
    "ParameterSampler",
    "SimOutput",
    "PARALYSIS_TO_INFECTION_RATIO",
    "compute_afp_from_infections",
    "daily_to_weekly",
    "district_weights_by_burden",
    "run_single_trial",
    "score_zero_weeks",
    "zero_week_fraction",
]
```

---

## `calibration/params.py`

```python
"""Parameter definitions and Latin Hypercube Sampling for LASER calibration.

Prior distributions
-------------------
beta              : U(0.15, 0.30)   – per-contact transmission probability
gravity_k         : 10^U(-3, -1)    – gravity model scaling constant
seasonal_amplitude: U(1.0, 1.5)     – multiplicative peak-season factor
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CalibParams:
    """One calibration parameter set."""

    beta: float               # Transmission probability per infectious contact
    gravity_k: float          # Gravity model constant (log-uniform prior)
    seasonal_amplitude: float # Peak-to-trough seasonal multiplier
    trial_id: int = 0

    def to_dict(self) -> dict:
        return {
            "beta": self.beta,
            "gravity_k": self.gravity_k,
            "seasonal_amplitude": self.seasonal_amplitude,
        }

    def __repr__(self) -> str:
        return (
            f"CalibParams(trial={self.trial_id}, "
            f"beta={self.beta:.4f}, "
            f"gravity_k={self.gravity_k:.3e}, "
            f"seas_amp={self.seasonal_amplitude:.3f})"
        )


class ParameterSampler:
    """Draw parameter sets via Latin Hypercube Sampling (LHS).

    LHS ensures uniform coverage of each marginal range with no clustering,
    which is more efficient than pure Monte Carlo for calibration sweeps.

    Parameters
    ----------
    n_trials : int
        Number of parameter sets to generate.
    seed : int, optional
        RNG seed for reproducibility.
    """

    # Prior bounds -------------------------------------------------------
    BETA_LO, BETA_HI = 0.15, 0.30
    LOG_K_LO, LOG_K_HI = -3.0, -1.0       # gravity_k = 10^u, u ~ U(-3,-1)
    AMP_LO, AMP_HI = 1.0, 1.5
    # --------------------------------------------------------------------

    def __init__(self, n_trials: int, seed: int | None = None) -> None:
        self.n_trials = n_trials
        self.rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _lhs(self, n_params: int) -> np.ndarray:
        """Return an (n_trials, n_params) LHS design in the unit hypercube."""
        n = self.n_trials
        # One uniform draw per stratum [j/n, (j+1)/n) per dimension
        strata = (np.arange(n)[:, None] + self.rng.uniform(size=(n, n_params))) / n
        # Shuffle each column independently to break correlation
        for col in range(n_params):
            self.rng.shuffle(strata[:, col])
        return strata  # shape: (n_trials, n_params)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def draw(self) -> list[CalibParams]:
        """Sample n_trials parameter sets and return them as CalibParams objects."""
        u = self._lhs(n_params=3)

        beta_vals   = self.BETA_LO + u[:, 0] * (self.BETA_HI - self.BETA_LO)
        log_k_vals  = self.LOG_K_LO + u[:, 1] * (self.LOG_K_HI - self.LOG_K_LO)
        amp_vals    = self.AMP_LO + u[:, 2] * (self.AMP_HI - self.AMP_LO)

        return [
            CalibParams(
                beta=float(beta_vals[i]),
                gravity_k=float(10.0 ** log_k_vals[i]),
                seasonal_amplitude=float(amp_vals[i]),
                trial_id=i,
            )
            for i in range(self.n_trials)
        ]
```

---

## `calibration/scoring.py`

```python
"""AFP comparison and fitness metrics for LASER polio calibration.

Fitness metric — zero-week fraction (ZWF)
-----------------------------------------
For each district d, the ZWF is the proportion of weeks with zero AFP
cases.  This mirrors critical community size (CCS) analysis: districts
whose effective population falls below the CCS exhibit frequent
transmission fade-outs (high ZWF), while large or well-connected
districts sustain near-endemic transmission (low ZWF).

Goodness of fit is the weighted mean squared error between the observed
and simulated per-district ZWF vectors.  A perfect calibration would
reproduce the observed spatial pattern of fade-out frequency.

AFP / infection relationship
----------------------------
Roughly 1 in 200 poliovirus infections produces clinically detectable
acute flaccid paralysis (AFP).  AFP surveillance therefore captures a
small, approximately constant fraction of total infections.  When the
LASER model reports infection incidence rather than paralysis directly,
this ratio converts incidence to expected AFP counts before comparison.
"""

from __future__ import annotations

import numpy as np

# ~1 paralysis per 200 poliovirus infections (range 1:100 – 1:1000)
PARALYSIS_TO_INFECTION_RATIO: float = 1.0 / 200.0


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def daily_to_weekly(daily: np.ndarray) -> np.ndarray:
    """Aggregate daily counts into complete ISO-week bins.

    Parameters
    ----------
    daily : ndarray, shape (n_days, n_nodes)
        Daily case counts, one column per district.

    Returns
    -------
    ndarray, shape (n_complete_weeks, n_nodes)
        Weekly sums.  Trailing days that do not form a full week are dropped.
    """
    n_days, n_nodes = daily.shape
    n_weeks = n_days // 7
    # Trim to an exact multiple of 7 then reshape and sum across day-of-week axis
    return daily[: n_weeks * 7].reshape(n_weeks, 7, n_nodes).sum(axis=1)


def compute_afp_from_infections(new_infections: np.ndarray) -> np.ndarray:
    """Scale infection incidence to expected AFP counts using the 1:200 ratio.

    Use this only when the LASER model exposes new infections but not new
    paralysis events directly.

    Parameters
    ----------
    new_infections : ndarray, shape (n_days, n_nodes)
        Daily new infection counts from the simulation.

    Returns
    -------
    ndarray, shape (n_days, n_nodes)
        Expected AFP counts (fractional; round for integer comparisons).
    """
    return new_infections * PARALYSIS_TO_INFECTION_RATIO


# ---------------------------------------------------------------------------
# Fitness metric: zero-week fraction
# ---------------------------------------------------------------------------

def zero_week_fraction(weekly: np.ndarray) -> np.ndarray:
    """Compute the proportion of weeks with zero AFP cases per district.

    Parameters
    ----------
    weekly : ndarray, shape (n_weeks, n_nodes)

    Returns
    -------
    ndarray, shape (n_nodes,)
        Values in [0, 1].  1.0 = all weeks silent, 0.0 = never silent.
    """
    return (weekly == 0).mean(axis=0)


def score_zero_weeks(
    obs_weekly: np.ndarray,
    sim_weekly: np.ndarray,
    weights: np.ndarray | None = None,
) -> float:
    """Scalar goodness-of-fit score based on per-district zero-week fractions.

    Computes the weighted mean squared error (wMSE) between the observed
    and simulated ZWF vectors.  **Lower is better.**

    Parameters
    ----------
    obs_weekly : ndarray, shape (n_weeks, n_nodes)
        Observed AFP weekly counts from surveillance.
    sim_weekly : ndarray, shape (n_weeks, n_nodes)
        Simulated AFP weekly counts from LASER.
    weights : ndarray, shape (n_nodes,), optional
        Per-district weights.  If None, all districts are weighted equally.
        Pass ``district_weights_by_burden(obs_weekly)`` to up-weight
        high-transmission districts.

    Returns
    -------
    float
        wMSE of ZWF vectors across districts.
    """
    obs_zwf = zero_week_fraction(obs_weekly)
    sim_zwf = zero_week_fraction(sim_weekly)
    squared_errors = (obs_zwf - sim_zwf) ** 2

    if weights is None:
        return float(np.mean(squared_errors))
    return float(np.average(squared_errors, weights=weights))


def district_weights_by_burden(obs_weekly: np.ndarray) -> np.ndarray:
    """Weight districts proportionally to their total observed AFP burden.

    Districts with more confirmed AFP cases receive higher weight, reflecting
    the epidemiological priority of high-transmission areas in calibration.

    Parameters
    ----------
    obs_weekly : ndarray, shape (n_weeks, n_nodes)

    Returns
    -------
    ndarray, shape (n_nodes,)  – sums to 1.0
    """
    burden = obs_weekly.sum(axis=0).astype(float)
    total = burden.sum()
    if total == 0:
        # No observed cases anywhere: fall back to uniform weights
        n_nodes = obs_weekly.shape[1]
        return np.full(n_nodes, 1.0 / n_nodes)
    return burden / total
```

---

## `calibration/runner.py`

```python
"""Execute a single LASER simulation trial and extract AFP-comparable output.

LASER conventions used here
---------------------------
* ``lp.run_sim(**kwargs)`` is the top-level entry point (matching the
  laser-polio demo pattern).
* Parameters flow through a ``PropertySet`` (``sim.pars``); the kwargs
  accepted by ``run_sim`` map directly to ``pars`` fields.
* Results live on ``sim.results`` as 2-D arrays indexed ``[timestep, node]``
  following the ``(nt, n_nodes)`` shape convention used by ``results.S``,
  ``results.new_paralyzed``, etc.
* ``new_paralyzed`` (daily new AFP-equivalent events) is preferred over
  ``new_infections`` because the LASER model explicitly tracks paralysis.
  If only infections are available the 1:200 ratio is applied instead.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .params import CalibParams
from .scoring import compute_afp_from_infections, daily_to_weekly

try:
    import laser_polio as lp
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "laser_polio must be installed.  See https://laser.idmod.org/ "
        "for installation instructions."
    ) from exc


@dataclass
class SimOutput:
    """Calibration-relevant outputs from one LASER simulation."""

    trial_id: int
    afp_weekly: np.ndarray   # shape: (n_complete_weeks, n_nodes)
    n_days: int
    n_nodes: int


def run_single_trial(
    pars: CalibParams,
    base_config: dict,
) -> SimOutput:
    """Run one LASER simulation with the given calibration parameters.

    Parameters
    ----------
    pars : CalibParams
        The sampled calibration parameters for this trial.
    base_config : dict
        Fixed simulation configuration (regions, n_days, pop_scale, seed
        schedule, vaccination coverage, etc.) passed to ``lp.run_sim()``.
        Keys in ``base_config`` that conflict with calibration parameters
        (beta, gravity_k, seasonal_amplitude) are silently overridden.

    Returns
    -------
    SimOutput
        Weekly AFP counts ready for scoring.
    """
    sim_kwargs = {
        **base_config,
        # Calibration parameters override any matching base_config keys
        "beta": pars.beta,
        "gravity_k": pars.gravity_k,
        "seasonal_amplitude": pars.seasonal_amplitude,
    }

    sim = lp.run_sim(**sim_kwargs)
    results = sim.results

    # ------------------------------------------------------------------
    # Extract AFP-equivalent daily counts
    # Prefer explicit paralysis tracking (new_paralyzed) which the LASER
    # DiseaseState_ABM component logs each timestep.  Fall back to applying
    # the 1:200 paralysis ratio to raw infection counts.
    # ------------------------------------------------------------------
    if hasattr(results, "new_paralyzed"):
        daily_afp = np.asarray(results.new_paralyzed, dtype=float)
        # Shape: (n_days, n_nodes) — LASER (nt, n_nodes) convention
    elif hasattr(results, "new_infections"):
        daily_afp = compute_afp_from_infections(
            np.asarray(results.new_infections, dtype=float)
        )
    else:
        raise AttributeError(
            "sim.results must expose 'new_paralyzed' or 'new_infections'. "
            "Verify that DiseaseState_ABM is included in the component list "
            "and that logging is enabled."
        )

    n_days, n_nodes = daily_afp.shape
    afp_weekly = daily_to_weekly(daily_afp)

    return SimOutput(
        trial_id=pars.trial_id,
        afp_weekly=afp_weekly,
        n_days=n_days,
        n_nodes=n_nodes,
    )
```

---

## `calibration/loop.py`

```python
"""Calibration loop: sample → simulate → score → rank.

Design notes
------------
* ParameterSampler uses LHS, giving better prior coverage than random
  Monte Carlo in the same number of trials.
* Scoring is based on zero-week fraction (ZWF) wMSE across districts,
  weighted by observed AFP burden so that high-transmission districts
  (where calibration matters most) drive the fit.
* Results are returned sorted ascending by score; rank 0 is the best fit.
* The loop is intentionally serial so that each LASER trial can use all
  available cores internally (Numba parallel=True).  Pass n_trials ≥ 100
  for reliable parameter space coverage.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from .params import CalibParams, ParameterSampler
from .runner import SimOutput, run_single_trial
from .scoring import (
    daily_to_weekly,
    district_weights_by_burden,
    score_zero_weeks,
)


@dataclass
class CalibResult:
    """Outcome of a single calibration trial."""

    params: CalibParams
    score: float          # wMSE of zero-week fractions — lower is better
    sim_output: SimOutput
    elapsed_s: float      # Wall-clock seconds for this trial

    # Support sorting by score
    def __lt__(self, other: "CalibResult") -> bool:
        return self.score < other.score

    def __repr__(self) -> str:
        return (
            f"CalibResult(rank=?, score={self.score:.5f}, "
            f"{self.params!r})"
        )


class CalibrationStudy:
    """Drive a full calibration sweep for the LASER polio model.

    Parameters
    ----------
    obs_afp_daily : ndarray, shape (n_days, n_nodes)
        Observed AFP surveillance counts in daily bins aligned to the
        simulation calendar.  Districts with no reported cases should
        be zero (not NaN).
    base_config : dict
        Fixed LASER model configuration forwarded unchanged to every
        trial (e.g. regions, n_days, pop_scale, seed schedule, RI
        coverage).  Calibration parameters beta / gravity_k /
        seasonal_amplitude are injected on top.
    n_trials : int
        Number of parameter sets to evaluate.  ≥100 recommended.
    seed : int, optional
        Global RNG seed.  Set for reproducibility; leave None for a
        different draw each run.
    weighted : bool
        If True, weight district scores by total observed AFP burden
        (recommended).  If False, all districts are weighted equally.

    Examples
    --------
    >>> study = CalibrationStudy(
    ...     obs_afp_daily=afp_array,          # shape (n_days, n_districts)
    ...     base_config={
    ...         "regions": ["PAKISTAN"],
    ...         "n_days": 365 * 3,
    ...         "pop_scale": 1.0,
    ...         "vx_prob_ri": 0.55,
    ...     },
    ...     n_trials=200,
    ...     seed=42,
    ... )
    >>> results = study.run()
    >>> print(study.summary_table(k=5))
    """

    def __init__(
        self,
        obs_afp_daily: np.ndarray,
        base_config: dict,
        n_trials: int = 100,
        seed: int | None = None,
        weighted: bool = True,
    ) -> None:
        # Pre-aggregate observed AFP to weekly resolution once
        self.obs_afp_weekly = daily_to_weekly(obs_afp_daily)
        self.base_config = base_config
        self.n_trials = n_trials
        self.weighted = weighted

        self._sampler = ParameterSampler(n_trials=n_trials, seed=seed)
        self._district_weights: np.ndarray | None = (
            district_weights_by_burden(self.obs_afp_weekly) if weighted else None
        )
        self._results: list[CalibResult] = []

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def results(self) -> list[CalibResult]:
        """Trials sorted best-to-worst (ascending score).  Empty before run()."""
        return self._results

    # ------------------------------------------------------------------
    # Core loop
    # ------------------------------------------------------------------

    def run(self, verbose: bool = True) -> list[CalibResult]:
        """Execute all trials and return them ranked by goodness of fit.

        Each trial:
          1. Injects sampled beta / gravity_k / seasonal_amplitude into LASER.
          2. Runs the full simulation.
          3. Computes the zero-week-fraction wMSE against observed AFP.

        Returns
        -------
        list[CalibResult]
            Sorted ascending by score (index 0 = best fit).
        """
        param_sets: list[CalibParams] = self._sampler.draw()
        self._results = []

        width = len(str(self.n_trials))

        for pars in param_sets:
            if verbose:
                print(
                    f"  trial {pars.trial_id + 1:{width}d}/{self.n_trials}"
                    f"  beta={pars.beta:.4f}"
                    f"  gravity_k={pars.gravity_k:.3e}"
                    f"  seas_amp={pars.seasonal_amplitude:.3f}",
                    end="  ",
                    flush=True,
                )

            t0 = time.perf_counter()
            sim_out = run_single_trial(pars=pars, base_config=self.base_config)
            elapsed = time.perf_counter() - t0

            score = score_zero_weeks(
                obs_weekly=self.obs_afp_weekly,
                sim_weekly=sim_out.afp_weekly,
                weights=self._district_weights,
            )

            result = CalibResult(
                params=pars,
                score=score,
                sim_output=sim_out,
                elapsed_s=elapsed,
            )
            self._results.append(result)

            if verbose:
                print(f"score={score:.5f}  ({elapsed:.1f}s)")

        # Rank: lower score = better fit
        self._results.sort()
        return self._results

    # ------------------------------------------------------------------
    # Post-hoc analysis helpers
    # ------------------------------------------------------------------

    def top_k(self, k: int = 10) -> list[CalibResult]:
        """Return the k best-scoring trials (call after run())."""
        if not self._results:
            raise RuntimeError("Call run() before accessing top_k().")
        return self._results[:k]

    def summary_table(self, k: int = 20) -> np.ndarray:
        """Structured array of the top-k results for quick inspection.

        Fields: rank, trial_id, score, beta, gravity_k,
                seasonal_amplitude, elapsed_s

        Returns
        -------
        numpy structured array, length min(k, n_trials)
        """
        top = self.top_k(k)
        dtype = np.dtype([
            ("rank",               np.int32),
            ("trial_id",           np.int32),
            ("score",              np.float64),
            ("beta",               np.float64),
            ("gravity_k",          np.float64),
            ("seasonal_amplitude", np.float64),
            ("elapsed_s",          np.float32),
        ])
        arr = np.empty(len(top), dtype=dtype)
        for rank, r in enumerate(top):
            arr[rank] = (
                rank,
                r.params.trial_id,
                r.score,
                r.params.beta,
                r.params.gravity_k,
                r.params.seasonal_amplitude,
                r.elapsed_s,
            )
        return arr

    def parameter_array(self) -> np.ndarray:
        """Return all sampled parameters as a plain (n_trials, 3) float array.

        Columns: beta, gravity_k (linear), seasonal_amplitude.
        Rows are sorted by score (best first), matching self.results order.
        """
        if not self._results:
            raise RuntimeError("Call run() before accessing parameter_array().")
        rows = [
            [r.params.beta, r.params.gravity_k, r.params.seasonal_amplitude]
            for r in self._results
        ]
        return np.array(rows, dtype=float)
```

---

## `run_calibration.py`

```python
"""Entry point: load AFP data, configure base model, and run calibration sweep.

Usage
-----
    python run_calibration.py

Outputs
-------
* Console table of top-20 parameter sets ranked by zero-week-fraction wMSE.
* ``calibration_results.npz`` — full results array for post-processing.
"""

from __future__ import annotations

import numpy as np

from calibration import CalibrationStudy


# ---------------------------------------------------------------------------
# 1. Load observed AFP surveillance data
#    Expected shape: (n_days, n_nodes) — daily counts, one column per district.
#    Replace this block with your actual data loading logic.
# ---------------------------------------------------------------------------

def load_afp_data(path: str) -> np.ndarray:
    """Load AFP surveillance data from a .npz or .npy file.

    The returned array must have shape (n_days, n_nodes) where:
      * n_days  aligns with the simulation calendar (same start date, same
                timestep resolution — daily)
      * n_nodes matches the number of districts in base_config["regions"]
      * Values are non-negative integer counts (0 for zero-incidence weeks)
    """
    data = np.load(path)
    if isinstance(data, np.ndarray):
        return data.astype(float)
    # Structured / dict-style .npz
    return data["afp_daily"].astype(float)


# ---------------------------------------------------------------------------
# 2. Fixed (non-calibrated) simulation configuration
#    These kwargs are forwarded unchanged to lp.run_sim() for every trial.
# ---------------------------------------------------------------------------

BASE_CONFIG = {
    "regions":        ["PAKISTAN"],   # ADM1 region filter
    "start_year":     2019,
    "n_days":         365 * 3,        # 3-year burn-in + calibration window
    "pop_scale":      1.0,
    "migration_method": "gravity",    # Use gravity model (gravity_k is calibrated)
    "vx_prob_ri":     0.55,           # Routine immunization coverage
    "missed_frac":    0.10,           # Fraction missed by SIA campaigns
    "save_plots":     False,
    "save_data":      False,
}

# ---------------------------------------------------------------------------
# 3. Calibration settings
# ---------------------------------------------------------------------------

N_TRIALS = 200   # LHS sample size — increase for finer coverage
SEED     = 42    # Fix for reproducibility

# ---------------------------------------------------------------------------
# 4. Run
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # -- Load data ----------------------------------------------------------
    afp_daily = load_afp_data("data/pakistan_afp_daily.npy")
    n_days, n_nodes = afp_daily.shape
    print(f"AFP data: {n_days} days × {n_nodes} districts")

    # -- Sanity check -------------------------------------------------------
    assert n_days == BASE_CONFIG["n_days"], (
        f"AFP data length ({n_days} days) must match n_days "
        f"({BASE_CONFIG['n_days']}) in BASE_CONFIG."
    )

    # -- Calibration sweep --------------------------------------------------
    study = CalibrationStudy(
        obs_afp_daily=afp_daily,
        base_config=BASE_CONFIG,
        n_trials=N_TRIALS,
        seed=SEED,
        weighted=True,   # Up-weight high-AFP-burden districts
    )

    print(f"\nRunning {N_TRIALS} calibration trials (LHS)...\n")
    results = study.run(verbose=True)

    # -- Report top results -------------------------------------------------
    print("\n--- Top 20 parameter sets (ranked by zero-week-fraction wMSE) ---")
    table = study.summary_table(k=20)
    header = (
        f"{'Rank':>4}  {'Trial':>5}  {'Score':>8}  "
        f"{'beta':>6}  {'gravity_k':>10}  {'seas_amp':>8}  {'sec':>5}"
    )
    print(header)
    print("-" * len(header))
    for row in table:
        print(
            f"{row['rank']:>4d}  {row['trial_id']:>5d}  {row['score']:>8.5f}  "
            f"{row['beta']:>6.4f}  {row['gravity_k']:>10.3e}  "
            f"{row['seasonal_amplitude']:>8.4f}  {row['elapsed_s']:>5.1f}"
        )

    # -- Save full results --------------------------------------------------
    param_arr = study.parameter_array()
    score_arr = np.array([r.score for r in results])

    np.savez(
        "calibration_results.npz",
        parameters=param_arr,         # (n_trials, 3): beta, gravity_k, seas_amp
        scores=score_arr,             # (n_trials,): wMSE, sorted ascending
        obs_afp_weekly=study.obs_afp_weekly,
    )
    print(f"\nSaved calibration_results.npz")
    print(f"Best fit: score={results[0].score:.5f}  {results[0].params!r}")
```

---

## Design rationale

**Parameter sampling — LHS over pure MC**
`ParameterSampler._lhs()` stratifies each dimension into `n_trials` equally-probable intervals and draws exactly one point per interval, then shuffles columns. This eliminates clustering and gives deterministic marginal uniformity, which is standard practice for expensive simulation calibration.

**`gravity_k` on log scale**
Because `gravity_k` spans two decades (`10⁻³` – `10⁻¹`), the LHS samples on `log₁₀(k)` and exponentiates, ensuring equal coverage of each order of magnitude rather than over-sampling values near `0.1`.

**Zero-week fraction as fitness metric**
`score_zero_weeks` computes the weighted MSE between the observed and simulated ZWF vectors. A district that truly experiences frequent fade-outs (small population, low connectivity) should have a high ZWF in both data and model. This metric is strictly more diagnostic than comparing raw incidence totals, because it captures the *temporal structure* of silence rather than just total burden.

**`new_paralyzed` vs the 1:200 ratio**
The LASER `DiseaseState_ABM` component logs `new_paralyzed` directly. `runner.py` uses that array first. The `compute_afp_from_infections` fallback (×1/200) is available when using a model configuration that only outputs infection incidence, or when running against a stripped-down component list.

**Serial outer loop**
The calibration loop is intentionally serial so each LASER trial can exploit all cores via Numba's `parallel=True` JIT compilation. Parallelising the outer loop would over-subscribe the CPU and degrade per-trial throughput.
