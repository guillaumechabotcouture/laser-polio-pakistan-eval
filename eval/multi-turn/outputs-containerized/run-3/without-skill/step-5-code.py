"""
Four-patch SIR model with gravity-model spatial coupling + seasonal forcing
+ demographic vital rates (births & deaths) via LASER-style components.
LASER-framework style — runs with numpy only (no extra installs needed).

Patch layout (linear, 75 km spacing):
    Patch 0 ──75km── Patch 1 ──75km── Patch 2 ──75km── Patch 3

Gravity parameters: k=0.01, a=1, b=1, c=1.5
Row-normalisation cap: 15 % export fraction per patch per step.

Seasonal forcing (ValuesMap schedule):
    Days   0–90 : 1.30× baseline  (winter peak)
    Days  90–150: linear ramp down 1.30→0.70  (spring)
    Days 150–240: 0.70× baseline  (summer trough)
    Days 240–365: linear ramp up  0.70→1.30  (autumn→winter)

Vital rates (LASER BirthsByCBR / MortalityByCDR):
    CBR = 30 per 1 000 per year
    CDR = 10 per 1 000 per year
    Net growth ≈ 2 % per year  →  ~22 % over 10 years

LASER components implemented (numpy-only, no install required):
    ValuesMap      — piecewise-linear parameter schedule
    calc_capacity  — pre-allocate for expected peak population
                     signature: calc_capacity(population, nticks, cbr)
                     from laser_core.utils  (laser-core 0.4.0)
    BirthsByCBR    — births via crude birth rate
                     formula: daily_rate = (1 + CBR/1000)^(1/365) − 1
                     from laser.generic.vitaldynamics (laser-generic 1.0.0)
    MortalityByCDR — deaths via crude death rate
                     formula: p_death = 1 − exp(−(1 − (1−CDR/1000)^(1/365)))
                     from laser.generic.vitaldynamics (laser-generic 1.0.0)
"""

import numpy as np

# ═══════════════════════════════════════════════════════════════
# LASER-style components  (pure numpy, no external dependencies)
# ═══════════════════════════════════════════════════════════════

class ValuesMap:
    """
    Piecewise-linear parameter schedule keyed by simulation day.

    Mirrors the LASER ValuesMap interface:
        vm = ValuesMap({day: value, ...})
        vm[day]  →  interpolated float

    Anchor values are clamped at the boundaries (no extrapolation).
    """

    def __init__(self, mapping: dict):
        sorted_keys = sorted(mapping)
        self._days = np.array(sorted_keys, dtype=float)
        self._vals = np.array([mapping[k] for k in sorted_keys], dtype=float)

    def __getitem__(self, day: float) -> float:
        return float(np.interp(day, self._days, self._vals))


def calc_capacity(population, nticks: int, cbr: float, verbose: bool = False) -> int:
    """
    Calculate the population capacity after *nticks* time steps based on
    a constant birth rate.  Mirrors laser_core.utils.calc_capacity.

    In LASER's agent-based framework this is used to pre-allocate the flat
    agent array before the simulation begins (agents are never removed, so
    only births grow the array).

    Parameters
    ----------
    population : int or float
        Initial total population.
    nticks : int
        Number of time steps (days when dt=1).
    cbr : float
        Crude birth rate (births per 1 000 per year).
    verbose : bool
        If True, print growth summary.

    Returns
    -------
    int
        Estimated population capacity after *nticks* steps.
    """
    daily_rate = (cbr / 1_000.0) / 365.0          # simple daily birth rate
    capacity = int(population * (1.0 + daily_rate) ** nticks)
    if verbose:
        print(f"  calc_capacity: {int(population):,} → {capacity:,} "
              f"(births-only upper bound over {nticks:,} ticks)")
    return capacity


class BirthsByCBR:
    """
    Compartmental adaptation of LASER's BirthsByCBR
    (laser.generic.vitaldynamics, laser-generic 1.0.0).

    Formula (from LASER source):
        daily_rate = (1 + CBR/1000)^(1/365) − 1

    This gives exact annual compounding: compounding daily_rate for
    365 days recovers exactly CBR/1000 per year.

    New births enter the Susceptible compartment.

    Parameters
    ----------
    cbr : float
        Crude birth rate (births per 1 000 per year).
    """

    def __init__(self, cbr: float):
        self.cbr = cbr
        # Exact annually-compounded daily rate (LASER formula)
        self._daily_rate = float(np.power(1.0 + cbr / 1_000.0, 1.0 / 365.0) - 1.0)

    def __call__(self, population, dt: float = 1.0):
        """
        Return per-patch birth counts (expected value) for one time step.

        Parameters
        ----------
        population : array-like
            Current total population per patch.
        dt : float
            Time-step length in days.

        Returns
        -------
        ndarray
            Expected births per patch this step.
        """
        return self._daily_rate * np.asarray(population, dtype=float) * dt


class MortalityByCDR:
    """
    Compartmental adaptation of LASER's MortalityByCDR
    (laser.generic.vitaldynamics, laser-generic 1.0.0).

    Formula (from LASER source):
        annual_survival  = 1 − CDR/1000
        daily_survival   = annual_survival^(1/365)
        daily_mort_rate  = 1 − daily_survival
        p_death          = 1 − exp(−daily_mort_rate)

    Deaths are removed proportionally from every compartment (S, I, R).

    Parameters
    ----------
    cdr : float
        Crude death rate (deaths per 1 000 per year).
    """

    def __init__(self, cdr: float):
        self.cdr = cdr
        annual_survival     = 1.0 - cdr / 1_000.0
        daily_survival      = float(np.power(annual_survival, 1.0 / 365.0))
        daily_mort_rate     = 1.0 - daily_survival
        # p_death = 1 − exp(−daily_mort_rate)  ≈ daily_mort_rate for small rates
        self._p_death_daily = float(-np.expm1(-daily_mort_rate))

    def __call__(self, dt: float = 1.0) -> float:
        """
        Return the fraction of each compartment that dies this step.

        For dt=1 day this equals p_death_daily; for other dt values the
        probability is scaled linearly (valid when p_death_daily << 1).

        Parameters
        ----------
        dt : float
            Time-step length in days.

        Returns
        -------
        float
            Per-step death fraction.
        """
        return self._p_death_daily * dt


# ═══════════════════════════════════════════════════════════════
# Seasonal-forcing schedule
# ═══════════════════════════════════════════════════════════════
seasonal_forcing = ValuesMap({
      0: 1.3,   # start of winter peak
     90: 1.3,   # end   of winter peak  → spring ramp begins
    150: 0.7,   # start of summer trough
    240: 0.7,   # end   of summer trough → autumn ramp begins
    365: 1.3,   # full annual cycle complete
})

# ═══════════════════════════════════════════════════════════════
# Patch geometry
# ═══════════════════════════════════════════════════════════════
NUM_PATCHES = 4
SPACING_KM  = 75.0
positions   = np.arange(NUM_PATCHES) * SPACING_KM       # [0, 75, 150, 225] km
dist        = np.abs(positions[:, None] - positions[None, :])   # (4×4) km

# ═══════════════════════════════════════════════════════════════
# Vital-rate parameters
# ═══════════════════════════════════════════════════════════════
CBR = 30.0    # births per 1 000 per year
CDR = 10.0    # deaths per 1 000 per year

births_process    = BirthsByCBR(cbr=CBR)
mortality_process = MortalityByCDR(cdr=CDR)

# ═══════════════════════════════════════════════════════════════
# Simulation timeline  (10 years)
# ═══════════════════════════════════════════════════════════════
YEARS  = 10
T      = YEARS * 365      # 3 650 days
dt     = 1.0              # daily time step
steps  = int(T / dt)      # 3 650 steps

# ═══════════════════════════════════════════════════════════════
# Initial population & calc_capacity pre-allocation
# ═══════════════════════════════════════════════════════════════
N_init = np.array([10_000, 8_000, 12_000, 6_000], dtype=float)

# calc_capacity: uses CBR only (matches LASER semantics — births grow
# the array; deaths do not shrink it in an agent-based model).
# Here we use it as a population budget / sanity check before the run.
capacity = calc_capacity(
    population = int(N_init.sum()),
    nticks     = steps,
    cbr        = CBR,
    verbose    = False,
)

print("=" * 62)
print("LASER calc_capacity  —  population budget")
print("=" * 62)
print(f"  Initial total N              : {N_init.sum():>10,.0f}")
print(f"  calc_capacity (births-only)  : {capacity:>10,d}  "
      f"(CBR={CBR:.0f}, {steps:,} ticks)")
net_expected = int(N_init.sum() * np.exp((CBR - CDR) / 1_000.0 * YEARS))
print(f"  Expected final N (net growth): {net_expected:>10,d}  "
      f"(continuous, CBR-CDR={CBR-CDR:.0f})")
print(f"  History arrays pre-allocated : {steps + 1:>10,d} × {NUM_PATCHES} "
      f"(one row per day)")
print()

# ═══════════════════════════════════════════════════════════════
# Gravity-model spatial coupling
# ═══════════════════════════════════════════════════════════════
k_grav = 0.01
a_exp  = 1.0
b_exp  = 1.0
c_exp  = 1.5
MAX_EXPORT_FRAC = 0.15

def build_gravity_coupling(populations, distances, k, a, b, c, max_export_frac):
    """
    Return (n×n) coupling matrix where coupling[i,j] is the fraction
    of patch i's population that moves to patch j each time step.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        G = np.where(
            distances > 0,
            k * (populations[:, None] ** a) * (populations[None, :] ** b)
              / (distances ** c),
            0.0,
        )
    frac       = G / populations[:, None]
    row_export = frac.sum(axis=1, keepdims=True)
    scale      = np.where(row_export > max_export_frac,
                          max_export_frac / row_export, 1.0)
    return frac * scale


def move(X, coupling):
    """Apply one step of gravity-driven spatial movement (population-conserving)."""
    outflow = coupling * X[:, None]
    net     = outflow.sum(axis=0) - outflow.sum(axis=1)
    return X + net


# ── Print initial coupling matrix ────────────────────────────────────────────
coupling_init = build_gravity_coupling(
    N_init, dist, k_grav, a_exp, b_exp, c_exp, MAX_EXPORT_FRAC
)
print("=" * 62)
print("Gravity coupling matrix  (fraction of source pop → dest)")
print("=" * 62)
print("        " + "".join(f"  Patch {j}" for j in range(NUM_PATCHES)))
for i in range(NUM_PATCHES):
    row_str = "".join(f"  {coupling_init[i, j]:.5f}" for j in range(NUM_PATCHES))
    print(f"Patch {i}{row_str}")
print()
print("Row sums (total export fraction per patch):")
for i in range(NUM_PATCHES):
    s   = coupling_init[i].sum()
    tag = "capped at 15%" if s >= MAX_EXPORT_FRAC - 1e-9 else "below cap"
    print(f"  Patch {i}: {s:.5f}  ({tag})")
print()

# ═══════════════════════════════════════════════════════════════
# SIR disease parameters
# ═══════════════════════════════════════════════════════════════
beta_baseline = 0.30   # baseline transmission rate  (day⁻¹)
gamma         = 0.10   # recovery rate               (day⁻¹)

# ═══════════════════════════════════════════════════════════════
# Initial conditions — seed patch 0
# ═══════════════════════════════════════════════════════════════
I = np.array([10.0, 0.0, 0.0, 0.0])
R = np.zeros(NUM_PATCHES)
S = N_init.copy() - I - R
N = N_init.copy()

# ═══════════════════════════════════════════════════════════════
# Pre-allocate history arrays  (steps+1 rows × NUM_PATCHES cols)
# ═══════════════════════════════════════════════════════════════
S_hist    = np.empty((steps + 1, NUM_PATCHES))
I_hist    = np.empty((steps + 1, NUM_PATCHES))
R_hist    = np.empty((steps + 1, NUM_PATCHES))
N_hist    = np.empty((steps + 1, NUM_PATCHES))
beta_hist = np.empty(steps + 1)

S_hist[0]    = S
I_hist[0]    = I
R_hist[0]    = R
N_hist[0]    = N
beta_hist[0] = beta_baseline * seasonal_forcing[0]

# ═══════════════════════════════════════════════════════════════
# Main simulation loop  (10 years = 3 650 daily steps)
# ═══════════════════════════════════════════════════════════════
for step in range(steps):

    day         = step * dt
    day_of_year = day % 365.0                    # fold into annual cycle
    beta_eff    = beta_baseline * seasonal_forcing[day_of_year]
    beta_hist[step] = beta_eff

    # 1. Within-patch SIR transitions (Euler step)
    new_inf = beta_eff * S * I / N * dt          # force of infection
    new_rec = gamma * I * dt

    S_new = S - new_inf
    I_new = I + new_inf - new_rec
    R_new = R + new_rec

    # 2. Vital rates
    #    BirthsByCBR  : new susceptibles born this step
    #    MortalityByCDR: fraction of each compartment that dies this step
    births     = births_process(N, dt)           # per-patch expected births
    death_frac = mortality_process(dt)           # scalar death fraction

    S_new = S_new * (1.0 - death_frac) + births  # births enter S
    I_new = I_new * (1.0 - death_frac)
    R_new = R_new * (1.0 - death_frac)

    # 3. Update total population (before spatial movement)
    N = S_new + I_new + R_new

    # 4. Rebuild gravity coupling with current (growing) N, then move
    coupling = build_gravity_coupling(
        N, dist, k_grav, a_exp, b_exp, c_exp, MAX_EXPORT_FRAC
    )
    S_new = move(S_new, coupling)
    I_new = move(I_new, coupling)
    R_new = move(R_new, coupling)

    # 5. Guard floating-point negatives; update state
    S = np.maximum(S_new, 0.0)
    I = np.maximum(I_new, 0.0)
    R = np.maximum(R_new, 0.0)
    N = S + I + R

    # 6. Store
    S_hist[step + 1]    = S
    I_hist[step + 1]    = I
    R_hist[step + 1]    = R
    N_hist[step + 1]    = N
    beta_hist[step + 1] = beta_baseline * seasonal_forcing[
        ((step + 1) * dt) % 365.0
    ]

# ═══════════════════════════════════════════════════════════════
# Seasonal forcing spot-check
# ═══════════════════════════════════════════════════════════════
print("=" * 62)
print("Seasonal forcing schedule  (ValuesMap spot-check)")
print("=" * 62)
print(f"  {'Day':>5}  {'Multiplier':>12}  {'Effective beta':>15}")
for d in [0, 45, 90, 120, 150, 195, 240, 300, 365]:
    m = seasonal_forcing[d % 365]
    print(f"  {d:>5}  {m:>12.4f}  {beta_baseline * m:>15.4f}")
print()

# ═══════════════════════════════════════════════════════════════
# Population growth verification  (the key new output)
# ═══════════════════════════════════════════════════════════════
N_total_start   = N_hist[0].sum()
N_total_end     = N_hist[-1].sum()
net_rate_annual = (CBR - CDR) / 1_000.0          # 0.02 = 2 %/yr
expected_end    = N_total_start * np.exp(net_rate_annual * YEARS)

print("=" * 62)
print("Population growth  (BirthsByCBR=30, MortalityByCDR=10)")
print("=" * 62)
print(f"  CBR = {CBR:.0f}/1000/yr  |  CDR = {CDR:.0f}/1000/yr  |  "
      f"net = {net_rate_annual*100:.1f}%/yr")
print(f"  BirthsByCBR  daily_rate : {births_process._daily_rate:.8f}/day")
print(f"  MortalityByCDR p_death  : {mortality_process._p_death_daily:.8f}/day")
print()
print(f"  {'Year':>5}  {'Total N':>12}  {'Annual growth':>14}")
print(f"  {0:>5}  {N_hist[0].sum():>12.1f}  {'(baseline)':>14}")
for yr in range(1, YEARS + 1):
    idx    = int(yr * 365 / dt)
    N_yr   = N_hist[idx].sum()
    N_prev = N_hist[int((yr - 1) * 365 / dt)].sum()
    ann_g  = (N_yr / N_prev - 1.0) * 100.0
    print(f"  {yr:>5}  {N_yr:>12.1f}  {ann_g:>13.2f}%")
print()
print(f"  Total population at START : {N_total_start:>12.1f}")
print(f"  Total population at END   : {N_total_end:>12.1f}")
print(f"  Expected (continuous exp) : {expected_end:>12.1f}")
actual_10yr_growth   = (N_total_end   / N_total_start - 1.0) * 100.0
expected_10yr_growth = (expected_end  / N_total_start - 1.0) * 100.0
implied_annual       = (N_total_end   / N_total_start) ** (1.0 / YEARS) - 1.0
print(f"  Actual 10-yr growth       : {actual_10yr_growth:>11.2f}%")
print(f"  Expected 10-yr growth     : {expected_10yr_growth:>11.2f}%")
print(f"  Implied annual growth rate: {implied_annual*100:>11.2f}%  "
      f"(target: {net_rate_annual*100:.1f}%)")
print()

# ═══════════════════════════════════════════════════════════════
# Per-patch summary
# ═══════════════════════════════════════════════════════════════
print("=" * 62)
print(f"Per-patch summary  ({YEARS}-year run with demographics)")
print("=" * 62)
print(f"  {'Patch':>6}  {'N(t=0)':>9}  {'Peak I':>9}  "
      f"{'Peak day':>9}  {'N(t=10yr)':>10}")
for p in range(NUM_PATCHES):
    peak_val = I_hist[:, p].max()
    peak_day = int(I_hist[:, p].argmax() * dt)
    print(f"  {p:>6}  {N_init[p]:>9.0f}  {peak_val:>9.1f}  "
          f"{peak_day:>9d}  {N_hist[-1, p]:>10.1f}")
print()

print("Final compartment totals across all patches:")
print(f"  S = {S_hist[-1].sum():.1f}")
print(f"  I = {I_hist[-1].sum():.1f}")
print(f"  R = {R_hist[-1].sum():.1f}")
print(f"  N = {(S_hist[-1] + I_hist[-1] + R_hist[-1]).sum():.1f}  "
      f"(initial N = {N_hist[0].sum():.1f})")
print()
print(f"Effective beta range over simulation: "
      f"{beta_hist.min():.4f} – {beta_hist.max():.4f}")