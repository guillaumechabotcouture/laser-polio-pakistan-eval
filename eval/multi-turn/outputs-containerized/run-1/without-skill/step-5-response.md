Continuing from the visualization section:

```python
one_year_ticks = np.arange(366)
beta_1yr       = beta_seasonal.as_array(one_year_ticks)

fig = plt.figure(figsize=(13, 14))
gs  = fig.add_gridspec(
    4, 2,
    height_ratios=[1, 1, 0.55, 0.45],
    hspace=0.52,
    wspace=0.30,
)

# — per-patch SEIR panels (fraction of current N) —
patch_axes = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (0, 1)]

for p, ax in enumerate(patch_axes):
    N_t        = N_ts[:, p]
    travel_pct = (1.0 - phi[p, p]) * 100.0
    for key, arr in [("S", nodes.S[:, p]), ("E", nodes.E[:, p]),
                     ("I", nodes.I[:, p]), ("R", nodes.R[:, p])]:
        ax.plot(years_axis, arr / N_t * 100, lw=1.5, label=key, color=colors[key])
    ax.set_title(f"Patch {p}  (N₀={patch_populations[p]:,})  "
                 f"[{travel_pct:.1f}% travel]", fontsize=10)
    ax.set_ylabel("Population (%)")
    ax.set_xlabel("Year")
    ax.set_xlim(0, N_YEARS)
    ax.set_ylim(0, 100)
    ax.legend(loc="center right", fontsize=9)
    ax.grid(True, alpha=0.3)

# — total population panel (full width, row 2) —
ax_pop = fig.add_subplot(gs[2, :])
for p in range(n_patches):
    ax_pop.plot(years_axis, N_ts[:, p] / 1_000, lw=1.5, color=pcolors[p],
                label=f"Patch {p}  (N₀={patch_populations[p]//1000}k)")
ax_pop.plot(years_axis, N_ts.sum(axis=1) / 1_000, lw=2.0, color="black", ls="--",
            label=f"Total  (N₀={patch_populations.sum()//1000}k)")
ax_pop.set_xlim(0, N_YEARS)
ax_pop.set_xlabel("Year")
ax_pop.set_ylabel("Population (thousands)")
ax_pop.set_title(
    f"Population growth  (CBR={CBR:.0f}, CDR={CDR:.0f} /1000/yr  →  ≈2%/yr net)",
    fontsize=10,
)
ax_pop.legend(loc="upper left", fontsize=9, ncol=2)
ax_pop.grid(True, alpha=0.3)

# — seasonal beta panel (one annual cycle, full width, row 3) —
ax_beta = fig.add_subplot(gs[3, :])
ax_beta.plot(one_year_ticks, beta_1yr, color="purple", lw=2,
             label="β(t) seasonal  (annual pattern — repeats via tick % 365)")
ax_beta.axhline(beta_baseline, color="grey", lw=1, ls="--",
                label=f"β baseline ({beta_baseline:.2f} /day)")
ax_beta.axvspan(  0,  90, alpha=0.10, color="steelblue",  label="winter peak (1.3×)")
ax_beta.axvspan(150, 240, alpha=0.10, color="tomato",     label="summer trough (0.7×)")
ax_beta.set_xlim(0, 365)
ax_beta.set_ylim(0, beta_baseline * 1.6)
ax_beta.set_xlabel("Day of year")
ax_beta.set_ylabel("β(t)  (/day)")
ax_beta.set_title("Seasonal forcing profile  (ValuesMap, one-year pattern)", fontsize=10)
ax_beta.legend(loc="upper right", fontsize=9, ncol=2)
ax_beta.grid(True, alpha=0.3)

fig.suptitle(
    f"4-Patch SEIR  |  R₀={R0}  latent={latent_period}d  "
    f"infectious={infectious_period}d  |  {N_YEARS}-year run\n"
    f"Gravity k={k_grav} a={int(a_grav)} b={int(b_grav)} c={c_grav} "
    f"cap={int(max_export*100)}%  |  "
    f"CBR={CBR:.0f}  CDR={CDR:.0f} per 1000/yr  (≈2%/yr net growth)",
    fontsize=10,
)

plt.savefig("seir_4patch_gravity_seasonal_vital.png", dpi=150, bbox_inches="tight")
print("Plot saved → seir_4patch_gravity_seasonal_vital.png")
plt.show()
```

---

Here is the **complete self-contained script** with every section assembled:

```python
"""
4-Patch SEIR Model with Gravity-Model Spatial Coupling + Seasonal Forcing
+ Vital Dynamics (Births by CBR, Deaths by CDR)
using the LASER framework (laser-core).

Disease parameters
------------------
  R0                = 5.0
  Latent period     = 4 days   → sigma = 1/4  = 0.25 /day  (E → I rate)
  Infectious period = 10 days  → gamma = 1/10 = 0.10 /day  (I → R rate)
  beta_baseline     = R0 * gamma = 0.50 /day

Vital dynamics
--------------
  CBR = 30 per 1000 per year  →  daily birth rate ≈ 8.22e-5
  CDR = 10 per 1000 per year  →  daily death rate ≈ 2.74e-5
  Net growth ≈ 2% per year  →  ~21.9% total growth over 10 yr  (e^{0.02×10})

  BirthsByCBR   : new_births_i = cbr_daily × N_i[t]     → S_i[t+1] += new_births_i
  MortalityByCDR: dead_X_i     = cdr_daily × X_i[t]     → X_i[t+1] -= dead_X_i
                                                            (X ∈ {S, E, I, R})

  calc_capacity: pre-computes expected maximum population per patch:
      cap_i = ceil( N0_i * exp(r * years) * buffer ),  r = (CBR-CDR)/1000

Seasonal forcing (ValuesMap, repeated annually via tick % 365)
--------------------------------------------------------------
  tick   0 -> 1.3x   tick  90 -> 1.3x   tick 150 -> 0.7x
  tick 240 -> 0.7x   tick 365 -> 1.3x

Patches (line, 75 km apart)
---------------------------
  Patch 0: 100 000   Patch 1: 200 000   Patch 2: 150 000   Patch 3: 80 000

Gravity coupling: G[i,j] = k*N_i^a*N_j^b/d^c, k=0.01, a=b=1, c=1.5, cap=15%

Initial conditions (all patches): S0=90%, E0=0%, I0=1%, R0=9%

Duration: 10 years (3 650 days)
"""

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point


# =========================================================================
# 0a.  LASER FRAMEWORK STUBS
#      Pure-Python replacements; drop when laser_core is available.
# =========================================================================

class PropertySet:
    """Parameter container with attribute-style access."""

    def __init__(self, d: dict) -> None:
        for k, v in d.items():
            setattr(self, k, v)

    def __repr__(self) -> str:
        pairs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"PropertySet({pairs})"


class LaserFrame:
    """Columnar state store for patch-level time series."""

    def __init__(self, count: int) -> None:
        self._count = count

    @property
    def count(self) -> int:
        return self._count

    def add_vector_property(self, name: str, size: int, dtype=np.float64) -> None:
        """Attach a (size, count) zero array as an attribute named *name*."""
        setattr(self, name, np.zeros((size, self._count), dtype=dtype))


# =========================================================================
# 0b.  VALUESMAP  --  piecewise-linear parameter time series
# =========================================================================

class ValuesMap:
    """
    Piecewise-linear parameter time series following the LASER ValuesMap API.

    Parameters
    ----------
    ticks  : array-like -- tick values at which the parameter is specified
    values : array-like -- corresponding parameter values

    Access
    ------
    vm[tick]  ->  float, linearly interpolated between the nearest knot ticks
    """

    def __init__(self, ticks, values) -> None:
        self._ticks  = np.asarray(ticks,  dtype=np.float64)
        self._values = np.asarray(values, dtype=np.float64)

    def __getitem__(self, tick) -> float:
        return float(np.interp(float(tick), self._ticks, self._values))

    def as_array(self, ticks) -> np.ndarray:
        return np.interp(np.asarray(ticks, dtype=np.float64),
                         self._ticks, self._values)


# =========================================================================
# 0c.  calc_capacity  --  pre-allocation helper (laser_core stub)
#
#     In laser_core, calc_capacity sizes the flat agent array for
#     individual-based models, accounting for population growth over
#     the full simulation horizon.  For this compartmental model it is
#     used to report expected peak populations and confirm the ~2%/yr
#     growth rate before the run begins.
#
#       r   = (CBR - CDR) / 1000  per year
#       cap = ceil( N0 * exp(r * years) * buffer )
# =========================================================================

def calc_capacity(
    populations: np.ndarray,
    cbr:         float,
    cdr:         float,
    years:       int,
    buffer:      float = 1.10,
) -> np.ndarray:
    """
    Return per-patch capacity (max expected population) for pre-allocation.

    Parameters
    ----------
    populations : (n_patches,) initial population counts
    cbr         : crude birth rate (births per 1000 per year)
    cdr         : crude death rate (deaths per 1000 per year)
    years       : simulation horizon in years
    buffer      : safety multiplier (default 1.10 = 10% headroom)
    """
    r = (cbr - cdr) / 1000.0   # net growth rate per year
    return np.ceil(
        populations.astype(np.float64) * np.exp(r * years) * buffer
    ).astype(np.int64)


# =========================================================================
# 1.  PARAMETERS
# =========================================================================

R0                = 5.0
latent_period     = 4.0    # days
infectious_period = 10.0   # days
CBR               = 30.0   # births per 1000 per year
CDR               = 10.0   # deaths per 1000 per year
N_YEARS           = 10

gamma         = 1.0 / infectious_period   # 0.10 /day
sigma         = 1.0 / latent_period       # 0.25 /day
beta_baseline = R0 * gamma               # 0.50 /day

params = PropertySet({
    "nticks":             N_YEARS * 365,
    "beta":               beta_baseline,
    "sigma":              sigma,
    "gamma":              gamma,
    "R0":                 R0,
    "latent_period":      latent_period,
    "infectious_period":  infectious_period,
    "cbr":                CBR,
    "cdr":                CDR,
    "n_years":            N_YEARS,
    "prng_seed":          20260101,
})


# =========================================================================
# 1b.  SEASONAL FORCING  --  annual ValuesMap, accessed via tick % 365
#
#     The knot ticks span one calendar year [0, 365].  For multi-year
#     runs the current day-of-year is obtained as  tick % 365  so the
#     seasonal pattern repeats identically each year.
# =========================================================================

_seasonal_ticks       = np.array([  0,   90,  150,  240,  365], dtype=np.float64)
_seasonal_multipliers = np.array([1.3,  1.3,  0.7,  0.7,  1.3], dtype=np.float64)

beta_seasonal = ValuesMap(
    ticks  = _seasonal_ticks,
    values = _seasonal_multipliers * beta_baseline,
)

print("=" * 62)
print("4-Patch SEIR  |  Gravity + Seasonal + Vital dynamics")
print("=" * 62)
print(f"\nDisease  : R0={R0}  beta_base={beta_baseline:.4f}"
      f"  sigma={sigma:.4f}  gamma={gamma:.4f}")
print(f"Vital    : CBR={CBR}/1000/yr  CDR={CDR}/1000/yr"
      f"  net={(CBR - CDR)/10:.1f}%/yr")
print(f"Run      : {N_YEARS} years  ({params.nticks} ticks)\n")


# =========================================================================
# 2.  SCENARIO  --  4-patch GeoDataFrame  (laser model convention)
# =========================================================================

patch_populations = np.array([100_000, 200_000, 150_000, 80_000], dtype=np.int64)
n_patches         = len(patch_populations)

S0    = np.floor(patch_populations * 0.90).astype(np.int64)
E0    = np.zeros(n_patches, dtype=np.int64)
I0    = np.floor(patch_populations * 0.01).astype(np.int64)
R0arr = patch_populations - S0 - E0 - I0

coords = [(-87.6, 41.8), (-88.0, 40.7), (-86.2, 39.8), (-89.4, 43.1)]

scenario = gpd.GeoDataFrame(
    {
        "nodeid":     np.arange(n_patches),
        "population": patch_populations,
        "S": S0, "E": E0, "I": I0, "R": R0arr,
    },
    geometry=[Point(lon, lat) for lon, lat in coords],
    crs="EPSG:4326",
)

# calc_capacity: confirm 10-yr growth budget before allocating arrays
capacities = calc_capacity(patch_populations, CBR, CDR, N_YEARS)
print("calc_capacity  (10-yr pre-allocation, 10% safety buffer)")
print(f"  {'Patch':>5}  {'N0':>9}  {'Capacity':>11}  {'Factor':>8}")
for p in range(n_patches):
    print(f"  {p:>5}  {patch_populations[p]:>9,}  "
          f"{capacities[p]:>11,}  {capacities[p] / patch_populations[p]:>7.3f}x")
print()

print("Initial scenario:")
print(scenario[["nodeid", "population", "S", "E", "I", "R"]].to_string(index=False))
print()


# =========================================================================
# 3.  GRAVITY-MODEL COUPLING MATRIX
#
#     Patches in a line 75 km apart -> positions [0, 75, 150, 225] km.
#     The coupling matrix is fixed at t=0 populations; vital dynamics do
#     not update phi each tick.
# =========================================================================

patch_spacing_km = 75.0
positions_km     = np.arange(n_patches, dtype=np.float64) * patch_spacing_km

k_grav     = 0.01
a_grav     = 1.0
b_grav     = 1.0
c_grav     = 1.5
max_export = 0.15

N_float   = patch_populations.astype(np.float64)
dist_km   = np.abs(positions_km[:, None] - positions_km[None, :])
dist_safe = np.where(dist_km == 0.0, np.inf, dist_km)

G_raw = (
    k_grav
    * (N_float[:, None] ** a_grav)
    * (N_float[None, :] ** b_grav)
    / (dist_safe ** c_grav)
)
np.fill_diagonal(G_raw, 0.0)

frac          = G_raw / N_float[:, None]
off_diag_sums = frac.sum(axis=1)
scale_vec     = np.where(off_diag_sums > max_export, max_export / off_diag_sums, 1.0)
frac          = frac * scale_vec[:, None]

phi = frac.copy()
np.fill_diagonal(phi, 1.0 - frac.sum(axis=1))

print("=== Gravity coupling matrix  phi[i,j] ===")
print("    (row = source patch, col = destination patch)")
print()
col_header = "         " + "".join(f"   P{j}" for j in range(n_patches))
print(col_header)
for i in range(n_patches):
    row_vals   = "".join(f"  {phi[i, j]:6.4f}" for j in range(n_patches))
    travel_pct = (1.0 - phi[i, i]) * 100.0
    print(f"  P{i}  {row_vals}   [{travel_pct:.1f}% travel]")
print()


# =========================================================================
# 4.  NODE STATE STORE  --  LaserFrame, nticks = 3 650
# =========================================================================

nticks = params.nticks   # 3 650

nodes = LaserFrame(n_patches)
for name in ("S", "E", "I", "R",
             "newly_exposed", "newly_infectious", "newly_recovered",
             "newly_born", "newly_dead"):
    nodes.add_vector_property(name, nticks + 1, dtype=np.float64)

nodes.S[0] = S0.astype(np.float64)
nodes.E[0] = E0.astype(np.float64)
nodes.I[0] = I0.astype(np.float64)
nodes.R[0] = R0arr.astype(np.float64)


# =========================================================================
# 5.  COMPONENTS
# =========================================================================

class BirthsByCBR:
    """
    Vital dynamics -- births at a constant crude birth rate.

    new_births_i = cbr_daily * N_i[t]
    S_i[t+1]    += new_births_i    (all newborns enter the susceptible pool)

    cbr_daily = CBR / (1000 * 365)
    """

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes     = nodes
        self.cbr_daily = params.cbr / (1_000.0 * 365.0)

    def step(self, tick: int) -> None:
        n          = self.nodes
        N          = n.S[tick] + n.E[tick] + n.I[tick] + n.R[tick]
        new_births = self.cbr_daily * N

        n.newly_born[tick + 1]  = new_births
        n.S[tick + 1]          += new_births


class MortalityByCDR:
    """
    Vital dynamics -- deaths at a constant crude death rate.

    Deaths are drawn proportionally from all four compartments:
        dead_X_i = cdr_daily * X_i[t]   for X in {S, E, I, R}
        X_i[t+1] -= dead_X_i

    cdr_daily = CDR / (1000 * 365)
    """

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes     = nodes
        self.cdr_daily = params.cdr / (1_000.0 * 365.0)

    def step(self, tick: int) -> None:
        n = self.nodes
        r = self.cdr_daily

        dead_S = r * n.S[tick]
        dead_E = r * n.E[tick]
        dead_I = r * n.I[tick]
        dead_R = r * n.R[tick]

        n.newly_dead[tick + 1]  = dead_S + dead_E + dead_I + dead_R
        n.S[tick + 1]          -= dead_S
        n.E[tick + 1]          -= dead_E
        n.I[tick + 1]          -= dead_I
        n.R[tick + 1]          -= dead_R


class TransmissionSE:
    """
    S -> E  :  spatially coupled, seasonally forced force of infection.

    Seasonal beta is looked up as beta_seasonal[tick % 365] so the
    annual knot-point pattern repeats correctly across all simulation years.

    FOI_i = beta(tick % 365) * (phi @ I_frac)[i]
    """

    def __init__(
        self,
        nodes:         LaserFrame,
        params:        PropertySet,
        phi:           np.ndarray,
        beta_seasonal: ValuesMap,
    ) -> None:
        self.nodes         = nodes
        self.params        = params
        self.phi           = phi
        self.beta_seasonal = beta_seasonal

    def step(self, tick: int) -> None:
        n      = self.nodes
        N      = n.S[tick] + n.E[tick] + n.I[tick] + n.R[tick]
        I_frac = n.I[tick] / N

        # Annual seasonal repeat: map tick to day-of-year for ValuesMap lookup
        beta_t = self.beta_seasonal[tick % 365]

        foi   = beta_t * (self.phi @ I_frac)
        new_E = foi * n.S[tick]

        n.newly_exposed[tick + 1]  = new_E
        n.S[tick + 1]             -= new_E
        n.E[tick + 1]             += new_E


class ExposedProgression:
    """E -> I  :  end of latent period; agents become infectious."""

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes  = nodes
        self.params = params

    def step(self, tick: int) -> None:
        n     = self.nodes
        new_I = self.params.sigma * n.E[tick]

        n.newly_infectious[tick + 1]  = new_I
        n.E[tick + 1]                -= new_I
        n.I[tick + 1]                += new_I


class InfectiousRecovery:
    """I -> R  :  recovery after the infectious period."""

    def __init__(self, nodes: LaserFrame, params: PropertySet) -> None:
        self.nodes  = nodes
        self.params = params

    def step(self, tick: int) -> None:
        n     = self.nodes
        new_R = self.params.gamma * n.I[tick]

        n.newly_recovered[tick + 1]  = new_R
        n.I[tick + 1]               -= new_R
        n.R[tick + 1]               += new_R


# =========================================================================
# 6.  MODEL ORCHESTRATOR
# =========================================================================

class SEIRModel:
    """
    Patch-level SEIR orchestrator following the laser.generic.Model design.

    _initialize_flows copies S, E, I, R from tick -> tick+1 at the start of
    every tick so each component only needs to add or subtract its own flow.
    Vital-dynamics flow arrays (newly_born, newly_dead) are written fresh
    each tick by their respective components and are not forward-filled.
    """

    def __init__(
        self,
        scenario:   gpd.GeoDataFrame,
        nodes:      LaserFrame,
        params:     PropertySet,
        components: list,
    ) -> None:
        self.scenario   = scenario
        self.nodes      = nodes
        self.params     = params
        self.components = components

    def _initialize_flows(self, tick: int) -> None:
        """Forward-fill compartment state tick -> tick+1 (laser.generic convention)."""
        n = self.nodes
        n.S[tick + 1] = n.S[tick]
        n.E[tick + 1] = n.E[tick]
        n.I[tick + 1] = n.I[tick]
        n.R[tick + 1] = n.R[tick]

    def run(self) -> None:
        for tick in range(self.params.nticks):
            self._initialize_flows(tick)
            for component in self.components:
                component.step(tick)
        print(f"Simulation complete: {self.params.nticks} days "
              f"({self.params.n_years} yr), {self.nodes.count} patches.\n")


# =========================================================================
# 7.  BUILD AND RUN
#
#     Component order (all read from [tick], accumulate into [tick+1]):
#       1. BirthsByCBR        S[t+1]         += births
#       2. MortalityByCDR     S,E,I,R[t+1]   -= deaths (proportional)
#       3. TransmissionSE     S[t+1] -= new_E;  E[t+1] += new_E
#       4. ExposedProgression E[t+1] -= new_I;  I[t+1] += new_I
#       5. InfectiousRecovery I[t+1] -= new_R;  R[t+1] += new_R
# =========================================================================

model = SEIRModel(
    scenario=scenario,
    nodes=nodes,
    params=params,
    components=[
        BirthsByCBR(nodes, params),
        MortalityByCDR(nodes, params),
        TransmissionSE(nodes, params, phi, beta_seasonal),
        ExposedProgression(nodes, params),
        InfectiousRecovery(nodes, params),
    ],
)

model.run()


# =========================================================================
# 8.  RESULTS SUMMARY  --  population growth check + final compartment state
# =========================================================================

# Total population time series:  shape (nticks+1, n_patches)
N_ts = nodes.S + nodes.E + nodes.I + nodes.R

N_start = N_ts[0].sum()
N_end   = N_ts[nticks].sum()

r_net         = (CBR - CDR) / 1000.0          # 0.020 /yr
N_theoretical = N_start * np.exp(r_net * N_YEARS)
annual_rate   = (N_end / N_start) ** (1.0 / N_YEARS) - 1.0

print("=== Population growth check (vital dynamics) ===")
print(f"  Total N  day     0 : {N_start:>14,.1f}")
print(f"  Total N  day {nticks:>4}  : {N_end:>14,.1f}")
print(f"  Simulated total growth  : {(N_end / N_start - 1) * 100:>7.3f}%")
print(f"  Theoretical  e^(rT)     : {(N_theoretical / N_start - 1) * 100:>7.3f}%"
      f"  (N_theor = {N_theoretical:,.0f})")
print(f"  Implied annual rate     : {annual_rate * 100:>7.4f}%  "
      f"(target ~2.0000%  since net = {CBR - CDR:.0f}/1000/yr)")
print()

print(f"=== Compartment fractions at day {nticks} ===")
header = (f"{'Patch':>5}  {'N0':>8}  {'N_end':>10}  "
          f"{'S %':>7}  {'E %':>6}  {'I %':>6}  {'R %':>7}")
print(header)
print("-" * len(header))
for p in range(n_patches):
    Ne = N_ts[nticks, p]
    Sf = nodes.S[nticks, p] / Ne * 100
    Ef = nodes.E[nticks, p] / Ne * 100
    If = nodes.I[nticks, p] / Ne * 100
    Rf = nodes.R[nticks, p] / Ne * 100
    print(f"{p:>5}  {patch_populations[p]:>8,}  {Ne:>10,.0f}  "
          f"{Sf:>6.2f}%  {Ef:>5.2f}%  {If:>5.2f}%  {Rf:>6.2f}%")

print()
print("=== Peak infectious (across all 10 years) ===")
for p in range(n_patches):
    peak_tick = int(nodes.I[:, p].argmax())
    peak_I    = nodes.I[peak_tick, p]
    Ne_peak   = N_ts[peak_tick, p]
    yr        = peak_tick // 365 + 1
    print(f"  Patch {p} (N0={patch_populations[p]:,}): "
          f"peak I = {peak_I:,.0f}"
          f"  ({peak_I / Ne_peak * 100:.1f}% of N)"
          f"  day {peak_tick}  (year {yr})")
print()


# =========================================================================
# 9.  VISUALISATION
#
#   Layout  (4 rows x 2 cols, gridspec):
#     Rows 0-1 : 2x2 per-patch SEIR fraction time series  (x-axis = years)
#     Row  2   : total population per patch + all-patch total
#     Row  3   : seasonal beta profile  (one annual cycle; repeats via %365)
# =========================================================================

years_axis = np.arange(nticks + 1) / 365.0
colors     = {"S": "steelblue", "E": "darkorange", "I": "firebrick", "R": "forestgreen"}
pcolors    = ["steelblue", "darkorange", "firebrick", "forestgreen"]

one_year_ticks = np.arange(366)
beta_1yr       = beta_seasonal.as_array(one_year_ticks)

fig = plt.figure(figsize=(13, 14))
gs  = fig.add_gridspec(
    4, 2,
    height_ratios=[1, 1, 0.55, 0.45],
    hspace=0.52,
    wspace=0.30,
)

# -- per-patch SEIR panels (fraction of current N) --
patch_axes = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (0, 1)]

for p, ax in enumerate(patch_axes):
    N_t        = N_ts[:, p]
    travel_pct = (1.0 - phi[p, p]) * 100.0
    for key, arr in [("S", nodes.S[:, p]), ("E", nodes.E[:, p]),
                     ("I", nodes.I[:, p]), ("R", nodes.R[:, p])]:
        ax.plot(years_axis, arr / N_t * 100, lw=1.5, label=key, color=colors[key])
    ax.set_title(f"Patch {p}  (N0={patch_populations[p]:,})  "
                 f"[{travel_pct:.1f}% travel]", fontsize=10)
    ax.set_ylabel("Population (%)")
    ax.set_xlabel("Year")
    ax.set_xlim(0, N_YEARS)
    ax.set_ylim(0, 100)
    ax.legend(loc="center right", fontsize=9)
    ax.grid(True, alpha=0.3)

# -- total population panel (full width, row 2) --
ax_pop = fig.add_subplot(gs[2, :])
for p in range(n_patches):
    ax_pop.plot(years_axis, N_ts[:, p] / 1_000, lw=1.5, color=pcolors[p],
                label=f"Patch {p}  (N0={patch_populations[p] // 1000}k)")
ax_pop.plot(years_axis, N_ts.sum(axis=1) / 1_000, lw=2.0, color="black", ls="--",
            label=f"Total  (N0={patch_populations.sum() // 1000}k)")
ax_pop.set_xlim(0, N_YEARS)
ax_pop.set_xlabel("Year")
ax_pop.set_ylabel("Population (thousands)")
ax_pop.set_title(
    f"Population growth  (CBR={CBR:.0f}, CDR={CDR:.0f} /1000/yr  ->  ~2%/yr net)",
    fontsize=10,
)
ax_pop.legend(loc="upper left", fontsize=9, ncol=2)
ax_pop.grid(True, alpha=0.3)

# -- seasonal beta panel (one annual cycle, full width, row 3) --
ax_beta = fig.add_subplot(gs[3, :])
ax_beta.plot(one_year_ticks, beta_1yr, color="purple", lw=2,
             label="beta(t) seasonal  (annual pattern -- repeats via tick % 365)")
ax_beta.axhline(beta_baseline, color="grey", lw=1, ls="--",
                label=f"beta baseline ({beta_baseline:.2f} /day)")
ax_beta.axvspan(  0,  90, alpha=0.10, color="steelblue",  label="winter peak (1.3x)")
ax_beta.axvspan(150, 240, alpha=0.10, color="tomato",     label="summer trough (0.7x)")
ax_beta.set_xlim(0, 365)
ax_beta.set_ylim(0, beta_baseline * 1.6)
ax_beta.set_xlabel("Day of year")
ax_beta.set_ylabel("beta(t)  (/day)")
ax_beta.set_title("Seasonal forcing profile  (ValuesMap, one-year pattern)", fontsize=10)
ax_beta.legend(loc="upper right", fontsize=9, ncol=2)
ax_beta.grid(True, alpha=0.3)

fig.suptitle(
    f"4-Patch SEIR  |  R0={R0}  latent={latent_period}d  "
    f"infectious={infectious_period}d  |  {N_YEARS}-year run\n"
    f"Gravity k={k_grav} a={int(a_grav)} b={int(b_grav)} c={c_grav} "
    f"cap={int(max_export * 100)}%  |  "
    f"CBR={CBR:.0f}  CDR={CDR:.0f} per 1000/yr  (~2%/yr net growth)",
    fontsize=10,
)

plt.savefig("seir_4patch_gravity_seasonal_vital.png", dpi=150, bbox_inches="tight")
print("Plot saved -> seir_4patch_gravity_seasonal_vital.png")
plt.show()
```

**Key additions explained:**

| Addition | Location | What it does |
|---|---|---|
| `calc_capacity` | Section 0c | Stub mirroring `laser_core.calc_capacity`; computes `ceil(N0 * exp(r*T) * 1.1)` per patch to show pre-allocation headroom |
| `BirthsByCBR` | Section 5 | Adds `cbr_daily * N[t]` new susceptibles each tick; `cbr_daily = 30/(1000*365)` |
| `MortalityByCDR` | Section 5 | Removes `cdr_daily * X[t]` from each of S, E, I, R each tick; `cdr_daily = 10/(1000*365)` |
| `nticks = 3650` | Sections 1 & 4 | 10 years × 365 days |
| `tick % 365` | `TransmissionSE.step` | Wraps the annual seasonal forcing correctly for years 2–10 |
| Population check | Section 8 | Prints start/end N and compares to `N0 * e^(0.02*10)` — should show ~21.9% total / ~2.00%/yr |

The growth check output will look like:
```
Total N  day     0 :        530,000.0
Total N  day 3650  :        647,xxx.x
Simulated total growth  :  22.xxx%
Theoretical  e^(rT)     :  22.140%  (N_theor = 647,341)
Implied annual rate     :   2.0xxx%  (target ~2.0000%)
```
The small deviation from theory arises because the discrete daily Euler steps approximate the continuous exponential — the agreement is tight at ~2%/yr.
