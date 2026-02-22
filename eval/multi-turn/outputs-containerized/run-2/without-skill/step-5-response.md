LASER is not installed — same situation as `ValuesMap` in the original. I'll implement all three inline, mirroring the LASER contract.

```python
#!/usr/bin/env python3
"""
SEIR 4-Patch Gravity Model with Seasonal Forcing and Vital Dynamics (10-year)
Institute for Disease Modeling – LASER framework demo

Vital dynamics:
  Births  : BirthsByCBR   – CBR = 30 per 1,000/year → enter as susceptibles
  Deaths  : MortalityByCDR – CDR = 10 per 1,000/year → proportional draw from S,E,I,R
  Net annual growth = (CBR − CDR) / 1,000 = 2% / year

Seasonal forcing (repeating annual cycle):
  Winter peak  (days 0–90 of each year)   : β × 1.3
  Summer trough(days 150–240 of each year): β × 0.7
  Linear transitions between breakpoints.

Note: ValuesMap, BirthsByCBR, MortalityByCDR, and calc_capacity are all
implemented inline below; they are not exported by the current public
laser.core API (v1.0.0).  Each interface mirrors the LASER framework contract.
"""

import numpy as np
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# ValuesMap – LASER-style piecewise-linear time-varying parameter map
# ─────────────────────────────────────────────────────────────────────────────

class ValuesMap:
    """
    Map from simulation day → parameter multiplier via piecewise-linear
    interpolation between caller-supplied breakpoints.

    Usage
    -----
    vm = ValuesMap(days=[0, 90, 150, 240], values=[1.3, 1.3, 0.7, 0.7])
    beta_today = beta_baseline * vm[t]
    """

    def __init__(self, days, values):
        self._days   = np.asarray(days,   dtype=np.float64)
        self._values = np.asarray(values, dtype=np.float64)
        if self._days.shape != self._values.shape:
            raise ValueError("days and values must have the same length")

    def __getitem__(self, day):
        return float(np.interp(day, self._days, self._values))

    def __call__(self, day):
        return self.__getitem__(day)

    def as_array(self, days):
        return np.interp(days, self._days, self._values)


# ─────────────────────────────────────────────────────────────────────────────
# BirthsByCBR – LASER-style births component
# ─────────────────────────────────────────────────────────────────────────────

class BirthsByCBR:
    """
    LASER-style births component using a crude birth rate (CBR).

    New births per tick = (CBR / 1,000 / 365) × N_patch × dt.
    All births enter the Susceptible (S) compartment.

    Parameters
    ----------
    cbr : float  Crude birth rate (births per 1,000 population per year).
    dt  : float  Simulation time step in days (default 1).
    """

    def __init__(self, cbr, dt=1.0):
        self.daily_rate = cbr / (1_000.0 * 365.0)
        self.dt = dt

    def apply(self, S, N_patch, t):
        """
        Add births to S[:, t+1].

        Parameters
        ----------
        S       : ndarray (N_patches, NTICKS+1)  susceptible state array
        N_patch : ndarray (N_patches,)           total population per patch at tick t
        t       : int                            current tick index

        Returns
        -------
        births : ndarray (N_patches,)  per-patch birth counts this tick
        """
        births = self.daily_rate * N_patch * self.dt
        S[:, t + 1] += births
        return births


# ─────────────────────────────────────────────────────────────────────────────
# MortalityByCDR – LASER-style mortality component
# ─────────────────────────────────────────────────────────────────────────────

class MortalityByCDR:
    """
    LASER-style mortality component using a crude death rate (CDR).

    Background deaths removed from each compartment X each tick:
        deaths_X = (CDR / 1,000 / 365) × X[:, t+1] × dt

    Deaths are drawn proportionally from S, E, I, R (not disease-induced).

    Parameters
    ----------
    cdr : float  Crude death rate (deaths per 1,000 population per year).
    dt  : float  Simulation time step in days (default 1).
    """

    def __init__(self, cdr, dt=1.0):
        self.daily_rate = cdr / (1_000.0 * 365.0)
        self.dt = dt

    def apply(self, S, E, I, R, t):
        """
        Reduce S, E, I, R at step t+1 by the per-tick survival fraction.

        Parameters
        ----------
        S, E, I, R : ndarray (N_patches, NTICKS+1)  compartment state arrays
        t          : int  current tick index
        """
        survival = 1.0 - self.daily_rate * self.dt
        S[:, t + 1] *= survival
        E[:, t + 1] *= survival
        I[:, t + 1] *= survival
        R[:, t + 1] *= survival


# ─────────────────────────────────────────────────────────────────────────────
# calc_capacity – LASER-style pre-allocation capacity estimator
# ─────────────────────────────────────────────────────────────────────────────

def calc_capacity(populations, years, cbr, cdr, buffer=0.05):
    """
    Estimate the maximum population per patch over the simulation horizon and
    return an integer capacity suitable for pre-allocating LASER agent arrays.

    Capacity = ceil( N_0 × (1 + net_annual_rate)^years × (1 + buffer) )

    For compartmental models this validates array headroom; in individual-based
    LASER models it determines the number of agent slots to allocate upfront.

    Parameters
    ----------
    populations : array-like  Initial population per patch.
    years       : int         Simulation duration (years).
    cbr         : float       Crude birth rate (per 1,000/year).
    cdr         : float       Crude death rate (per 1,000/year).
    buffer      : float       Fractional safety margin above projection (default 5%).

    Returns
    -------
    capacity : ndarray (int)  Pre-allocation size per patch.
    """
    net_annual_rate = (cbr - cdr) / 1_000.0
    projected_max   = np.asarray(populations, dtype=float) * (1.0 + net_annual_rate) ** years
    return np.ceil(projected_max * (1.0 + buffer)).astype(int)


# ─────────────────────────────────────────────────────────────────────────────
# Model parameters
# ─────────────────────────────────────────────────────────────────────────────

N_PATCHES       = 4
YEARS           = 10
NTICKS          = YEARS * 365            # 3,650 ticks (one per day)
DT              = 1                      # days per tick

POPULATIONS     = np.array([500_000, 200_000, 150_000, 100_000], dtype=float)
BETA_BASELINE   = 0.30                   # baseline transmission rate (/day)
SIGMA           = 1 / 5.0               # 1 / incubation period (days)
GAMMA           = 1 / 10.0              # 1 / infectious period (days)

CBR             = 30.0                   # crude birth rate  (per 1,000/year)
CDR             = 10.0                   # crude death rate  (per 1,000/year)

MAX_EXPORT_FRAC = 0.01                   # max daily export fraction per patch
GRAVITY_K       = 1.0                    # gravity model constant

I0              = 100                    # seed infectious individuals in patch 0

COLORS = {"S": "#2196F3", "E": "#FF9800", "I": "#F44336", "R": "#4CAF50"}

# ─────────────────────────────────────────────────────────────────────────────
# Pre-allocation via calc_capacity
# ─────────────────────────────────────────────────────────────────────────────

capacity = calc_capacity(POPULATIONS, YEARS, CBR, CDR)
print(f"calc_capacity ({YEARS}-yr horizon, CBR={CBR:.0f}, CDR={CDR:.0f} per 1,000/yr):")
for p in range(N_PATCHES):
    print(f"  Patch {p}: initial {POPULATIONS[p]:>9,.0f}  →  capacity {capacity[p]:>10,}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# Seasonal forcing profile using ValuesMap (annual cycle)
#
#   days   0– 90  → 1.3× (winter peak, flat)
#   days  90–150  → 1.3 → 0.7 (spring decline, linear)
#   days 150–240  → 0.7× (summer trough, flat)
#   days 240–300  → 0.7 → 1.0 (autumn recovery, linear)
#   days 300–365  → 1.0 → 1.3 (pre-winter rise, linear)
# ─────────────────────────────────────────────────────────────────────────────

seasonal_map = ValuesMap(
    days   = [  0,  90, 150, 240, 300, 365],
    values = [1.3, 1.3, 0.7, 0.7, 1.0, 1.3],
)

# ─────────────────────────────────────────────────────────────────────────────
# Gravity migration rate matrix (computed from initial populations; held fixed)
# ─────────────────────────────────────────────────────────────────────────────

lats = np.array([0.0, 0.0, 1.0, 1.0])
lons = np.array([0.0, 1.0, 0.0, 1.0])


def _haversine_km(i, j):
    R_earth = 6371.0
    phi1, phi2 = np.radians(lats[i]), np.radians(lats[j])
    dphi = np.radians(lats[j] - lats[i])
    dlam = np.radians(lons[j] - lons[i])
    a = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2) ** 2
    return 2 * R_earth * np.arcsin(np.sqrt(max(a, 1e-12)))


dist = np.array([[_haversine_km(i, j) for j in range(N_PATCHES)]
                 for i in range(N_PATCHES)])

raw_flow = np.zeros((N_PATCHES, N_PATCHES))
for i in range(N_PATCHES):
    for j in range(N_PATCHES):
        if i != j:
            raw_flow[i, j] = GRAVITY_K * POPULATIONS[i] * POPULATIONS[j] / dist[i, j] ** 2

mig_rates = np.zeros((N_PATCHES, N_PATCHES))
for i in range(N_PATCHES):
    total_out = raw_flow[i].sum()
    if total_out > 0:
        scale = min(1.0, MAX_EXPORT_FRAC * POPULATIONS[i] / total_out)
        for j in range(N_PATCHES):
            if i != j:
                mig_rates[i, j] = raw_flow[i, j] / POPULATIONS[i] * scale

mig_out_rates = mig_rates.sum(axis=1)   # shape (N_PATCHES,)

# ─────────────────────────────────────────────────────────────────────────────
# SEIR state arrays  –  shape: (N_PATCHES, NTICKS + 1)
# ─────────────────────────────────────────────────────────────────────────────

S = np.zeros((N_PATCHES, NTICKS + 1))
E = np.zeros((N_PATCHES, NTICKS + 1))
I = np.zeros((N_PATCHES, NTICKS + 1))
R = np.zeros((N_PATCHES, NTICKS + 1))

S[0, 0] = POPULATIONS[0] - I0
I[0, 0] = I0
for p in range(1, N_PATCHES):
    S[p, 0] = POPULATIONS[p]

# ─────────────────────────────────────────────────────────────────────────────
# Vital dynamics components
# ─────────────────────────────────────────────────────────────────────────────

births_component = BirthsByCBR(CBR,  dt=DT)
deaths_component = MortalityByCDR(CDR, dt=DT)

# ─────────────────────────────────────────────────────────────────────────────
# Simulation – deterministic Euler
#   Step order each tick: SEIR transitions → migration → births → deaths
#
# Births base population on N at tick t (before transitions) to keep the
# vital-dynamics step independent of the disease dynamics step.
# Deaths multiply each compartment at t+1 by the daily survival fraction.
# ─────────────────────────────────────────────────────────────────────────────

for t in range(NTICKS):
    # Seasonal β cycles annually
    day_of_year = t % 365
    beta_t = BETA_BASELINE * seasonal_map[day_of_year]

    # Total population per patch at t (FOI denominator + births denominator)
    Np = S[:, t] + E[:, t] + I[:, t] + R[:, t]
    Np = np.where(Np > 0, Np, 1.0)

    # SEIR transitions (vectorised over all patches simultaneously)
    foi   = beta_t * I[:, t] / Np
    new_E = foi   * S[:, t] * DT
    new_I = SIGMA * E[:, t] * DT
    new_R = GAMMA * I[:, t] * DT

    S[:, t + 1] = S[:, t] - new_E
    E[:, t + 1] = E[:, t] + new_E - new_I
    I[:, t + 1] = I[:, t] + new_I - new_R
    R[:, t + 1] = R[:, t] + new_R

    # Gravity migration (operator-split)
    # inflow[i]  = (mig_rates.T @ c)[i]  = Σ_j mig_rates[j,i] * c[j]
    # outflow[i] = mig_out_rates[i] * c[i]
    for comp in (S, E, I, R):
        c = comp[:, t + 1].copy()           # snapshot before in-place update
        comp[:, t + 1] += mig_rates.T @ c - mig_out_rates * c

    # Births: new susceptibles enter at CBR, based on population at tick t
    births_component.apply(S, Np, t)

    # Deaths: background mortality drawn proportionally from all compartments
    deaths_component.apply(S, E, I, R, t)

# ─────────────────────────────────────────────────────────────────────────────
# Derived totals and summary statistics
# ─────────────────────────────────────────────────────────────────────────────

years_axis      = np.arange(NTICKS + 1) / 365.0
S_total         = S.sum(axis=0)
E_total         = E.sum(axis=0)
I_total         = I.sum(axis=0)
R_total         = R.sum(axis=0)
N_total         = S_total + E_total + I_total + R_total

peak_tick       = int(np.argmax(I_total))
seasonal_factor = seasonal_map.as_array(np.arange(NTICKS + 1) % 365)

pop_start         = N_total[0]
pop_end           = N_total[-1]
actual_growth_ann = (pop_end / pop_start) ** (1.0 / YEARS) - 1.0
expected_growth   = (CBR - CDR) / 1_000.0

print(f"Total population at start (day 0)       : {pop_start:>12,.0f}")
print(f"Total population at end   (day {NTICKS:,}) : {pop_end:>12,.0f}")
print(f"Expected annual growth rate              :     {expected_growth:.2%}")
print(f"Actual   annual growth rate              :     {actual_growth_ann:.2%}")
print()
print(f"Peak infectious : {I_total[peak_tick]:>12,.0f}  "
      f"at year {peak_tick / 365:.2f}  (tick {peak_tick})")

# ─────────────────────────────────────────────────────────────────────────────
# Visualisation – 2 × 3 panel figure (x-axis in years)
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle(
    f"SEIR 4-Patch Gravity Model – {YEARS}-year run  |  "
    f"CBR={CBR:.0f}  CDR={CDR:.0f} per 1,000/yr  |  "
    "Seasonal β (winter ×1.3  days 0–90 · summer ×0.7  days 150–240)",
    fontsize=12, fontweight="bold",
)

# ── Per-patch SEIR panels ─────────────────────────────────────────────────────
patch_axes = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]
for p, ax in enumerate(patch_axes):
    ax.plot(years_axis, S[p] / 1e3, label="S", color=COLORS["S"], lw=1.5)
    ax.plot(years_axis, E[p] / 1e3, label="E", color=COLORS["E"], lw=1.5)
    ax.plot(years_axis, I[p] / 1e3, label="I", color=COLORS["I"], lw=1.5)
    ax.plot(years_axis, R[p] / 1e3, label="R", color=COLORS["R"], lw=1.5)
    # Secondary axis: seasonal β multiplier (annual cycle)
    ax2 = ax.twinx()
    ax2.plot(years_axis, seasonal_factor, color="purple", lw=0.8, ls=":", alpha=0.5)
    ax2.set_ylim(0.4, 1.8)
    ax2.set_ylabel("β scale", fontsize=7, color="purple")
    ax2.tick_params(axis="y", colors="purple", labelsize=6)
    ax.set_title(f"Patch {p}  (initial pop {POPULATIONS[p] / 1e3:.0f} k)", fontsize=10)
    ax.set_xlabel("Year")
    ax.set_ylabel("Population (thousands)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, YEARS)

# ── Aggregate panel with total-population growth curve ───────────────────────
ax_agg = axes[1, 2]
ax_agg.plot(years_axis, S_total / 1e3, label="S",        color=COLORS["S"], lw=2)
ax_agg.plot(years_axis, E_total / 1e3, label="E",        color=COLORS["E"], lw=2)
ax_agg.plot(years_axis, I_total / 1e3, label="I",        color=COLORS["I"], lw=2)
ax_agg.plot(years_axis, R_total / 1e3, label="R",        color=COLORS["R"], lw=2)
ax_agg.plot(years_axis, N_total / 1e3, label="N (total)",
            color="black", lw=1.5, ls="--", alpha=0.7)
ax_agg.axvline(peak_tick / 365.0, color="gray", ls="--", lw=1,
               label=f"I peak yr {peak_tick / 365:.1f}")
ax_agg.set_title("All Patches Combined", fontsize=11)
ax_agg.set_xlabel("Year")
ax_agg.set_ylabel("Population (thousands)")
ax_agg.legend(loc="upper left", fontsize=8)
ax_agg.grid(True, alpha=0.3)
ax_agg.set_xlim(0, YEARS)

# ── Migration rate heatmap ────────────────────────────────────────────────────
ax_mig = axes[0, 2]
im = ax_mig.imshow(mig_rates * 100, cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=MAX_EXPORT_FRAC * 100)
ax_mig.set_title("Daily Migration Rates (%)", fontsize=11)
ax_mig.set_xticks(range(N_PATCHES))
ax_mig.set_yticks(range(N_PATCHES))
ax_mig.set_xticklabels([f"P{i}" for i in range(N_PATCHES)], fontsize=9)
ax_mig.set_yticklabels([f"P{i}" for i in range(N_PATCHES)], fontsize=9)
ax_mig.set_xlabel("Destination patch", fontsize=9)
ax_mig.set_ylabel("Origin patch", fontsize=9)
plt.colorbar(im, ax=ax_mig, label="%/day")
for i in range(N_PATCHES):
    for j in range(N_PATCHES):
        val = mig_rates[i, j] * 100
        txt_color = "white" if val > MAX_EXPORT_FRAC * 50 else "black"
        ax_mig.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=8, color=txt_color)

plt.tight_layout()
plt.savefig("seir_4patch_gravity_results.png", dpi=150, bbox_inches="tight")
print("\nFigure saved: seir_4patch_gravity_results.png")
plt.show()
```

---

**Key additions and design notes:**

**Three new inline classes (mirroring LASER contracts):**

- `BirthsByCBR(cbr, dt)` — converts CBR to a per-day rate (`cbr / 1000 / 365`), then each tick adds `rate × N_patch × dt` new individuals into `S[:, t+1]`. Uses population *at tick t* (before transitions) as the denominator, keeping the vital-dynamics step independent of the disease step.

- `MortalityByCDR(cdr, dt)` — multiplies each compartment at `t+1` by `1 − rate × dt` after SEIR transitions and migration, so deaths are drawn proportionally from S, E, I, R (background mortality, not disease-induced).

- `calc_capacity(populations, years, cbr, cdr, buffer=0.05)` — projects maximum population per patch via compound growth at the net rate `(CBR−CDR)/1000`, adds a 5 % safety margin, and returns integer capacities. In LASER agent-based models this determines slot pre-allocation; here it validates headroom and is printed at startup.

**Simulation order each tick:** SEIR transitions → gravity migration → births → deaths. This is operator-split Euler; ordering effects are O(dt²) and negligible at dt = 1 day.

**Growth rate check:** Discrete Euler gives `N(t+1) = N(t) × (1 + CBR_daily) × (1 − CDR_daily)`, which compounds to ~2.02 % per year — you should see the actual printed rate land within a few basis points of the expected 2.00 %.

**Visualization:** x-axis is now in years (0–10); the dotted purple β-scale curve repeats the annual seasonal cycle on each per-patch panel; the aggregate panel adds `N (total)` as a dashed black line to make the ~22 % cumulative population growth visible directly on the plot.
