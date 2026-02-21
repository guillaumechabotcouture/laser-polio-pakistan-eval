#!/usr/bin/env python3
"""
LASER Spatial Polio Transmission Model for Pakistan
====================================================
Large-scale Array-based Simulation for Epidemiological Research (LASER)
framework implementing a stochastic SEIR compartmental model.

Features:
  - 20 Pakistani districts with heterogeneous populations (60k–2M)
  - Gravity-model spatial coupling between districts
  - Monsoon-season transmission forcing (July–October peak)
  - Per-district OPV routine immunization (pop-weighted avg ≈ 80%)
  - Supplementary Immunization Activities (SIA) every 6 months
  - Periodic importation of infections (first 5 years)
  - 20-year simulation with 10-year burn-in
  - Tau-leap stochastic transitions (Poisson / binomial draws)
"""

import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

RNG = np.random.default_rng(seed=42)

# -- Time ------------------------------------------------------------------
DAYS_PER_YEAR = 365
YEARS         = 20
BURNIN_YEARS  = 10
N_DAYS        = YEARS * DAYS_PER_YEAR          # 7 300

# -- SEIR parameters -------------------------------------------------------
R0                = 6.0
LATENT_PERIOD     = 3.0                         # days  (E → I)
INFECTIOUS_PERIOD = 28.0                        # days  (I → R)
SIGMA  = 1.0 / LATENT_PERIOD                   # 1/3   day⁻¹
GAMMA  = 1.0 / INFECTIOUS_PERIOD               # 1/28  day⁻¹
BETA0  = R0 * GAMMA                            # ≈ 0.2143 day⁻¹

# -- Demographics (Pakistan) -----------------------------------------------
BIRTH_RATE = 25.0 / 1000.0 / DAYS_PER_YEAR     # per capita per day
DEATH_RATE =  7.0 / 1000.0 / DAYS_PER_YEAR

# -- Seasonal forcing ------------------------------------------------------
SEASONAL_AMP     = 0.30                         # ±30 % around β₀
MONSOON_PEAK_DOY = 240                          # ~late August

# -- Gravity-model spatial coupling ----------------------------------------
GRAV_K       = 0.005
GRAV_ALPHA   = 0.6                              # origin-pop exponent
GRAV_BETA_P  = 0.6                              # dest-pop exponent
GRAV_GAMMA_D = 2.0                              # distance-decay exponent
MAX_EXT_FOI  = 0.10                             # cap on external FOI fraction

# -- Vaccination ------------------------------------------------------------
ROUTINE_EFF  = 0.50           # per-series OPV seroconversion (type 1, tropical)

SIA_INTERVAL = 180            # days between SIA campaigns
SIA_EFF      = 0.50           # per-dose seroconversion during SIA
SIA_COV_BASE = 0.70           # baseline SIA coverage (for 80 % routine areas)
UNDER5_FRAC  = 0.25           # fraction of susceptibles that are under-5 (SIA target)

# -- Importation ------------------------------------------------------------
IMPORT_INTERVAL = 60          # days between import events
IMPORT_N        = 3           # infections seeded per event
IMPORT_YEARS    = 5           # importation window (years from start)

# ═══════════════════════════════════════════════════════════════════════════════
# DISTRICT DATA
# ═══════════════════════════════════════════════════════════════════════════════

#           name               pop       lat    lon    rout_cov
_DISTRICTS = [
    ("Karachi",           2_000_000,  24.86, 67.01,  0.89),
    ("Lahore",            1_800_000,  31.52, 74.35,  0.91),
    ("Faisalabad",        1_200_000,  31.42, 73.08,  0.87),
    ("Rawalpindi",        1_000_000,  33.60, 73.05,  0.88),
    ("Peshawar",            900_000,  34.01, 71.58,  0.82),
    ("Multan",              800_000,  30.20, 71.47,  0.85),
    ("Hyderabad",           700_000,  25.39, 68.37,  0.87),
    ("Gujranwala",          650_000,  32.16, 74.19,  0.88),
    ("Islamabad",           600_000,  33.69, 73.04,  0.92),
    ("Quetta",              500_000,  30.18, 67.00,  0.65),
    ("Bannu",               400_000,  32.99, 70.60,  0.55),
    ("D.I. Khan",           350_000,  31.83, 70.90,  0.50),
    ("Khyber",              300_000,  34.10, 71.15,  0.45),
    ("North Waziristan",    250_000,  33.00, 69.85,  0.30),
    ("Lakki Marwat",        200_000,  32.61, 70.91,  0.50),
    ("Mohmand",             150_000,  34.50, 71.30,  0.40),
    ("South Waziristan",    120_000,  32.30, 69.65,  0.25),
    ("Tank",                100_000,  32.22, 70.38,  0.35),
    ("Chaman",               80_000,  30.92, 66.45,  0.55),
    ("Zhob",                 60_000,  31.34, 69.45,  0.45),
]

N_DIST = len(_DISTRICTS)
NAMES  = [d[0] for d in _DISTRICTS]
POP0   = np.array([d[1] for d in _DISTRICTS], dtype=np.float64)
COORDS = np.array([[d[2], d[3]] for d in _DISTRICTS])

# Per-district routine coverage (pop-weighted average ≈ 80 %)
DIST_ROUTINE_COV = np.array([d[4] for d in _DISTRICTS])
DIST_EFF_ROUTINE = DIST_ROUTINE_COV * ROUTINE_EFF   # effective fraction of births immunised

# Per-district SIA coverage — scales with routine access
DIST_SIA_COV  = np.clip(DIST_ROUTINE_COV * (SIA_COV_BASE / 0.80), 0.10, 0.90)
DIST_SIA_PROB = DIST_SIA_COV * SIA_EFF              # per-targeted-child immunisation probability


# ═══════════════════════════════════════════════════════════════════════════════
# SPATIAL STRUCTURE
# ═══════════════════════════════════════════════════════════════════════════════

def haversine_matrix(coords):
    """Pairwise great-circle distances (km) via Haversine formula."""
    lat = np.radians(coords[:, 0])
    lon = np.radians(coords[:, 1])
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = (np.sin(dlat / 2.0) ** 2
         + np.cos(lat[:, None]) * np.cos(lat[None, :])
         * np.sin(dlon / 2.0) ** 2)
    return 6_371.0 * 2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def gravity_matrix(pops, dist_km):
    """Row-normalised gravity-model coupling matrix."""
    n = len(pops)
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                d = max(dist_km[i, j], 1.0)
                C[i, j] = (GRAV_K
                           * pops[i] ** GRAV_ALPHA
                           * pops[j] ** GRAV_BETA_P
                           / d ** GRAV_GAMMA_D)
    row_sums = np.maximum(C.sum(axis=1, keepdims=True), 1e-12)
    C *= MAX_EXT_FOI / row_sums
    return C


# ═══════════════════════════════════════════════════════════════════════════════
# SIMULATION ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

def seasonal_beta(day):
    """β(t) with monsoon-season cosine forcing."""
    doy = day % DAYS_PER_YEAR
    return BETA0 * (1.0 + SEASONAL_AMP
                    * np.cos(2.0 * np.pi * (doy - MONSOON_PEAK_DOY)
                             / DAYS_PER_YEAR))


def _pois(rate):
    """Safe Poisson draw (clamp negative rates)."""
    return RNG.poisson(np.clip(rate, 0.0, None))


def _clamp(drawn, pool):
    """Ensure drawn ≤ available pool (integer)."""
    return np.minimum(drawn, np.maximum(pool, 0).astype(np.int64))


def run():
    """Execute the LASER spatial SEIR simulation; return weekly incidence."""

    # ── spatial structure ────────────────────────────────────────────────────
    dist_km = haversine_matrix(COORDS)
    C       = gravity_matrix(POP0, dist_km)

    # ── initial conditions (near endemic equilibrium S* ≈ 1/R₀) ─────────────
    S = POP0 / R0                                             # ~16.7 % susceptible
    I = np.maximum(POP0 * 5e-4, 5.0)                         # small seed
    E = np.maximum(POP0 * 2e-4, 3.0)
    R = POP0 - S - E - I
    N = POP0.copy()

    # ── output arrays ────────────────────────────────────────────────────────
    n_weeks   = YEARS * 52
    week_inc  = np.zeros((n_weeks, N_DIST))
    week_acc  = np.zeros(N_DIST)
    day_in_wk = 0
    wk        = 0

    # ── print header ─────────────────────────────────────────────────────────
    wtd_cov = np.dot(POP0, DIST_ROUTINE_COV) / POP0.sum()
    print("=" * 74)
    print("   LASER Spatial Polio Transmission Model — Pakistan")
    print("=" * 74)
    print(f"   Districts        : {N_DIST}")
    print(f"   Total population : {POP0.sum():,.0f}")
    print(f"   R₀ = {R0}   σ = 1/{LATENT_PERIOD:.0f} d   "
          f"γ = 1/{INFECTIOUS_PERIOD:.0f} d   β₀ = {BETA0:.4f}/d")
    print(f"   Routine vacc     : {wtd_cov*100:.1f}% pop-wtd cov × "
          f"{ROUTINE_EFF*100:.0f}% eff  (range "
          f"{DIST_ROUTINE_COV.min()*100:.0f}%–{DIST_ROUTINE_COV.max()*100:.0f}%)")
    print(f"   SIA campaigns    : every {SIA_INTERVAL} d, "
          f"targeting under-5 susceptibles ({UNDER5_FRAC*100:.0f}% of S)")
    print(f"   Importation      : {IMPORT_N} inf / {IMPORT_INTERVAL} d "
          f"for first {IMPORT_YEARS} y")
    print(f"   Seasonal forcing : ±{SEASONAL_AMP*100:.0f}% "
          f"(monsoon peak DOY {MONSOON_PEAK_DOY})")
    print(f"   Simulation       : {YEARS} y ({N_DAYS:,} d), "
          f"burn-in {BURNIN_YEARS} y")
    print("-" * 74)

    # ── daily loop ───────────────────────────────────────────────────────────
    for day in range(N_DAYS):

        bt = seasonal_beta(day)

        # Force of infection: local prevalence + gravity-coupled neighbours
        prev = np.where(N > 0, I / N, 0.0)
        foi  = bt * (prev + C @ prev)

        # Stochastic SEIR transitions (tau-leap)
        new_E = _clamp(_pois(foi * S), S)
        new_I = _clamp(_pois(SIGMA * E), E)
        new_R = _clamp(_pois(GAMMA * I), I)

        # Demographics
        births = _pois(BIRTH_RATE * N)
        d_S = _clamp(_pois(DEATH_RATE * S), S)
        d_E = _clamp(_pois(DEATH_RATE * E), E)
        d_I = _clamp(_pois(DEATH_RATE * I), I)
        d_R = _clamp(_pois(DEATH_RATE * R), R)

        # Routine immunization (district-specific coverage × efficacy)
        vacc_births = RNG.binomial(births, DIST_EFF_ROUTINE)
        susc_births = births - vacc_births

        # SIA pulse — targets under-5 susceptibles, once per cycle
        sia = np.zeros(N_DIST, dtype=np.int64)
        if day > 0 and day % SIA_INTERVAL == 0:
            u5_target = np.maximum(UNDER5_FRAC * S, 0.0).astype(np.int64)
            sia = RNG.binomial(u5_target, DIST_SIA_PROB)
            sia = _clamp(sia, S)

        # Importation (first IMPORT_YEARS only)
        imp = np.zeros(N_DIST, dtype=np.int64)
        if (day // DAYS_PER_YEAR) < IMPORT_YEARS and day % IMPORT_INTERVAL == 0:
            for _ in range(IMPORT_N):
                imp[RNG.integers(N_DIST)] += 1

        # ── update compartments ──────────────────────────────────────────────
        S += susc_births - new_E - d_S - sia
        E += new_E - new_I - d_E
        I += new_I - new_R - d_I + imp      # imports arrive as infectious
        R += new_R + vacc_births + sia - d_R

        S = np.maximum(S, 0.0)
        E = np.maximum(E, 0.0)
        I = np.maximum(I, 0.0)
        R = np.maximum(R, 0.0)
        N = S + E + I + R

        # ── accumulate weekly incidence ──────────────────────────────────────
        week_acc += new_I
        day_in_wk += 1
        if day_in_wk == 7:
            if wk < n_weeks:
                week_inc[wk] = week_acc
            week_acc[:] = 0.0
            day_in_wk = 0
            wk += 1

        # annual progress
        if (day + 1) % DAYS_PER_YEAR == 0:
            y   = (day + 1) // DAYS_PER_YEAR
            tag = " [burn-in]" if y <= BURNIN_YEARS else ""
            s_pct = 100.0 * S.sum() / N.sum()
            print(f"   Year {y:>2d}{tag:<11s}  "
                  f"Pop {N.sum():>12,.0f}   "
                  f"S {s_pct:>5.1f}%   "
                  f"Active inf {I.sum():>8,.0f}")

    recorded = min(wk, n_weeks)
    print("-" * 74)
    print(f"   Simulation complete — {recorded} weeks recorded.\n")
    return week_inc[:recorded]


# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS & OUTPUT
# ═══════════════════════════════════════════════════════════════════════════════

def analyse(week_inc):
    """Post-burn-in analysis of weekly incidence."""

    burnin_wk = BURNIN_YEARS * 52
    total_wk  = week_inc.shape[0]
    if burnin_wk >= total_wk:
        print("   !! Burn-in exceeds recorded data; using all weeks.")
        burnin_wk = 0
    data = week_inc[burnin_wk:]
    n_wk = data.shape[0]

    print("=" * 74)
    print(f"   POST BURN-IN RESULTS   (years {BURNIN_YEARS+1}–{YEARS},  "
          f"{n_wk} weeks)")
    print("=" * 74)

    # ── per-district table ───────────────────────────────────────────────────
    hdr = (f"   {'District':<20s} {'Cov%':>4s} {'Pop₀':>9s} "
           f"{'Total':>10s} {'Mean/wk':>8s} {'Peak/wk':>8s} "
           f"{'Zero wks':>8s} {'%Zero':>6s}")
    print(hdr)
    print("   " + "─" * 71)

    total_all    = 0.0
    sum_pct_zero = 0.0

    for j in range(N_DIST):
        col  = data[:, j]
        tot  = col.sum()
        mn   = col.mean()
        pk   = col.max()
        zw   = int(np.sum(col == 0))
        pz   = 100.0 * zw / n_wk
        total_all    += tot
        sum_pct_zero += pz
        cov  = DIST_ROUTINE_COV[j] * 100
        print(f"   {NAMES[j]:<20s} {cov:>3.0f}% {POP0[j]:>9,.0f} "
              f"{tot:>10,.0f} {mn:>8.1f} {pk:>8.0f} "
              f"{zw:>8d} {pz:>5.1f}%")

    print("   " + "─" * 71)

    agg = data.sum(axis=1)
    print(f"   {'ALL DISTRICTS':<20s}      {POP0.sum():>9,.0f} "
          f"{total_all:>10,.0f} {agg.mean():>8.1f} {agg.max():>8.0f}")

    # ── zero-incidence summary ───────────────────────────────────────────────
    agg_zero = int(np.sum(agg == 0))
    mean_pct = sum_pct_zero / N_DIST

    print()
    print(f"   National zero-incidence weeks  : "
          f"{agg_zero} / {n_wk}  ({100.0 * agg_zero / n_wk:.1f}%)")
    print(f"   Mean per-district %-zero weeks : {mean_pct:.1f}%")

    # ── yearly breakdown ─────────────────────────────────────────────────────
    print()
    print(f"   {'Year':>6s} {'Cases':>10s} "
          f"{'Districts w/ cases':>20s} {'Natl zero wks':>14s}")
    print("   " + "─" * 56)

    for y in range(BURNIN_YEARS, YEARS):
        y_off = y - BURNIN_YEARS
        w0 = y_off * 52
        w1 = min(w0 + 52, n_wk)
        if w0 >= n_wk:
            break
        yr  = data[w0:w1]
        cas = yr.sum()
        act = int(np.sum(yr.sum(axis=0) > 0))
        yz  = int(np.sum(yr.sum(axis=1) == 0))
        print(f"   {y + 1:>6d} {cas:>10,.0f} "
              f"{act:>14d} / {N_DIST:<4d} {yz:>14d} / 52")

    # ── per-district zero-incidence bar chart ────────────────────────────────
    print()
    print("   PROPORTION OF ZERO-INCIDENCE WEEKS PER DISTRICT")
    print("   " + "─" * 56)
    bar_w = 40
    for j in range(N_DIST):
        col = data[:, j]
        zw  = int(np.sum(col == 0))
        frac = zw / n_wk
        filled = int(bar_w * frac)
        bar = "█" * filled + "░" * (bar_w - filled)
        print(f"   {NAMES[j]:<20s} {bar} {100.0 * frac:>5.1f}%")

    return data


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    weekly = run()
    analyse(weekly)
