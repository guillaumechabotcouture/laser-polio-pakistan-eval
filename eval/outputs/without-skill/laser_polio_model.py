#!/usr/bin/env python3
"""
Extended LASER Polio Metapopulation Model
==========================================

SEIR polio model on a linear array of patches, built with the LASER
(Light Agent Spatial modeling for ERadication) framework from IDM.

Extensions over the basic model:
  1. Gravity-model migration network  (b=0.5, c=1.5)
  2. Monsoon seasonal forcing          (Jul-Oct ×1.3, Dec-Mar ×0.7)
  3. Row-normalized coupling           (≤15 % FOI from neighbours)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ═══════════════════════════════════════════════════════════════════════════
#  Lightweight stand-ins for laser.core types
#  (allows running without the laser package installed)
# ═══════════════════════════════════════════════════════════════════════════
class PropertySet:
    """Attribute-access wrapper around a dict (mirrors laser.core.PropertySet)."""

    def __init__(self, mapping=None, **kwargs):
        d = dict(mapping or {}, **kwargs)
        object.__setattr__(self, "_data", d)

    def __getattr__(self, name):
        try:
            return self._data[name]
        except KeyError:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        self._data[name] = value

    def __repr__(self):
        return f"PropertySet({self._data})"


class LaserFrame:
    """Minimal tabular container (mirrors laser.core.LaserFrame).

    Supports ``add_scalar_property`` so patch metadata can be attached
    by name, exactly as the real LASER frame expects.
    """

    def __init__(self, capacity, initial_count=0):
        self.capacity = capacity
        self.count = initial_count

    def add_scalar_property(self, name, dtype=np.float64):
        setattr(self, name, np.zeros(self.capacity, dtype=dtype))


# ═══════════════════════════════════════════════════════════════════════════
#  Configuration
# ═══════════════════════════════════════════════════════════════════════════
NUM_PATCHES = 10
SPACING_KM  = 50.0            # inter-patch distance

params = PropertySet({
    "num_patches":  NUM_PATCHES,
    "nticks":       365 * 3,   # 3-year simulation, daily steps

    # --- SEIR disease parameters (wild poliovirus) ---
    "beta":   0.3,             # baseline transmission rate  (day⁻¹)
    "sigma":  1.0 / 5.0,      # E→I rate  (latent period  ≈ 5 d)
    "gamma":  1.0 / 21.0,     # I→R rate  (infectious     ≈ 21 d)
    "mu":     1.0 / (70*365), # birth = death rate  (stable pop)

    # --- Gravity model ---
    "gravity_k": 1.0,         # scaling constant  (cancels in normalisation)
    "gravity_a": 1.0,         # origin-population exponent
    "gravity_b": 0.5,         # destination-population exponent
    "gravity_c": 1.5,         # distance-decay exponent

    # --- Coupling constraint ---
    "max_foi_import": 0.15,   # max off-diagonal row-sum (15 % cap)

    # --- Seasonal forcing ---
    "seasonal_amplitude": 0.3,

    # --- Vaccination ---
    "ri_coverage":       0.80,          # routine immunisation coverage (OPV at 6 wk)
    "sia_coverage":      0.90,          # SIA coverage (children 0-5 yr)
    "sia_interval":      182,           # days between SIAs (≈ 6 months)
    "fraction_under5":   5.0 / 70.0,    # demographic proxy for 0-5 yr fraction
    "ve_opv":            0.50,          # vaccine efficacy (leaky; 50 % reduction in FOI)
    "waning_rate":       1.0 / (3*365), # V→S rate (OPV mucosal immunity ≈ 3 yr)
})


# ═══════════════════════════════════════════════════════════════════════════
#  Extension 1 – Gravity-model migration network
# ═══════════════════════════════════════════════════════════════════════════
def build_gravity_network(populations, positions_km, k, a, b, c):
    """Classic gravity model (off-diagonal only).

        G[i,j] = k · Nᵢᵃ · Nⱼᵇ / dᵢⱼᶜ     (i ≠ j)

    Equivalent to ``laser.core.migration.gravity`` with the same arguments.
    """
    n = len(populations)
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d_ij = abs(positions_km[i] - positions_km[j])
            G[i, j] = k * populations[i]**a * populations[j]**b / d_ij**c
    return G


# ═══════════════════════════════════════════════════════════════════════════
#  Extension 3 – Row-normalise the coupling matrix
# ═══════════════════════════════════════════════════════════════════════════
def row_normalize_network(raw_network, max_rowsum):
    """Row-normalise so no patch imports more than *max_rowsum* of its FOI.

    1. Scale each row so  Σ_{j≠i} C[i,j] ≤ max_rowsum.
    2. Set  C[i,i] = 1 − Σ_{j≠i} C[i,j]  (local retention ≥ 1−max_rowsum).

    Mirrors ``laser.core.migration.row_normalizer``.

    The resulting matrix C satisfies:
      • every row sums to exactly 1
      • off-diagonal entries encode the fraction of a patch's FOI that
        originates from neighbours (≤ max_rowsum)
      • the diagonal encodes the fraction retained locally (≥ 0.85)
    """
    C = raw_network.copy()
    np.fill_diagonal(C, 0.0)

    for i in range(C.shape[0]):
        s = C[i].sum()
        if s > max_rowsum:
            C[i] *= max_rowsum / s

    # Local retention on the diagonal
    np.fill_diagonal(C, 1.0 - C.sum(axis=1))
    return C


# ═══════════════════════════════════════════════════════════════════════════
#  Extension 2 – Seasonal forcing profile
# ═══════════════════════════════════════════════════════════════════════════
def build_seasonal_profile(amplitude):
    """Build a 365-element seasonal-multiplier array.

    Jul–Oct (monsoon):  β × (1 + amplitude)  = 1.3 × baseline
    Dec–Mar (dry):      β × (1 − amplitude)  = 0.7 × baseline
    Transition months are linearly interpolated between mid-month anchors.
    """
    #              Jan   Feb   Mar   Apr        May  Jun         Jul   Aug   Sep   Oct   Nov  Dec
    mid_days  = [  15,   46,   74,  105,       135, 166,        196,  227,  258,  288,  319, 349]
    mult      = [1-amplitude, 1-amplitude, 1-amplitude,         # Dec–Mar trough
                 1 - 0.5*amplitude,                              # Apr  ramp-up
                 1.0,                                            # May  baseline
                 1 + 0.5*amplitude,                              # Jun  ramp-up
                 1+amplitude, 1+amplitude, 1+amplitude, 1+amplitude,  # Jul–Oct peak
                 1.0,                                            # Nov  ramp-down
                 1-amplitude]                                    # Dec  trough

    # Wrap endpoints for periodic interpolation across the year boundary
    ext_days = np.concatenate([[mid_days[-1] - 365], mid_days, [mid_days[0] + 365]])
    ext_mult = np.concatenate([[mult[-1]],           mult,     [mult[0]]])

    return np.interp(np.arange(365), ext_days, ext_mult)


# ═══════════════════════════════════════════════════════════════════════════
#  Build spatial scaffolding
# ═══════════════════════════════════════════════════════════════════════════
positions = np.arange(NUM_PATCHES) * SPACING_KM

rng       = np.random.default_rng(42)
pop_sizes = rng.integers(50_000, 500_000, size=NUM_PATCHES).astype(np.float64)

raw_network = build_gravity_network(
    pop_sizes, positions,
    k=params.gravity_k, a=params.gravity_a,
    b=params.gravity_b, c=params.gravity_c,
)
coupling = row_normalize_network(raw_network, params.max_foi_import)
seasonal = build_seasonal_profile(params.seasonal_amplitude)


# ═══════════════════════════════════════════════════════════════════════════
#  LASER model class
# ═══════════════════════════════════════════════════════════════════════════
class PolioModel:
    """Spatial SEIRV polio model following the LASER component pattern.

    Compartments
    ------------
    S : Susceptible
    E : Exposed (latent)
    I : Infectious
    R : Recovered via natural infection  (strong, long-lasting immunity)
    V : Vaccinated with OPV              (mucosal immunity, wanes → S)
    """

    def __init__(self, params, pop_sizes, coupling, seasonal):
        self.params   = params
        n, T          = params.num_patches, params.nticks

        # --- Patch-level metadata (LaserFrame) ---
        self.patches = LaserFrame(capacity=n, initial_count=n)
        self.patches.add_scalar_property("population", dtype=np.float64)
        self.patches.population[:] = pop_sizes

        # --- Compartment state vectors ---
        self.S = pop_sizes * 0.90
        self.E = pop_sizes * 0.02
        self.I = pop_sizes * 0.01
        self.R = pop_sizes * 0.07
        self.V = np.zeros(n)             # OPV-vaccinated (starts empty)

        # Seed extra cases in patch 0
        seed      = min(100.0, self.S[0] * 0.01)
        self.I[0] += seed
        self.S[0] -= seed

        # --- Spatial coupling & seasonal profile ---
        self.coupling = coupling
        self.seasonal = seasonal

        # --- Time-series output (patch × tick) ---
        self.hist = {c: np.zeros((n, T))
                     for c in ("S", "E", "I", "R", "V")}

        # --- Component pipeline ---
        # Order: transmission → progression → recovery → demography
        #        → vaccination (RI, SIA) → waning → record
        self.phases = [
            TransmissionPhase(self),
            LatentPhase(self),
            RecoveryPhase(self),
            DemographyPhase(self),
            RoutineImmunizationPhase(self),
            SIAPhase(self),
            WaningPhase(self),
            RecorderPhase(self),
        ]

    def run(self):
        for tick in range(self.params.nticks):
            for phase in self.phases:
                phase(self, tick)


# ═══════════════════════════════════════════════════════════════════════════
#  LASER components (step functions)
# ═══════════════════════════════════════════════════════════════════════════
class TransmissionPhase:
    """S → E  and  V → E  with gravity-coupled FOI and seasonal forcing.

    effective_prevalence_i = Σ_j  C[i,j] · I_j / N_j

    FOI_i(t) = β(t) · effective_prevalence_i

    where β(t) = β₀ · seasonal(day-of-year).

    Vaccinated individuals (V) experience a *reduced* force of infection
    scaled by (1 − ve_opv), modelling a "leaky" OPV vaccine.  Break-
    through infections move V → E (and are subsequently tracked through
    the normal E → I → R pathway, so re-infected individuals gain
    natural immunity).
    """

    def __init__(self, model):
        self.beta     = model.params.beta
        self.coupling = model.coupling
        self.seasonal = model.seasonal
        self.ve_opv   = model.params.ve_opv   # vaccine efficacy (leaky)

    def __call__(self, model, tick):
        beta_t = self.beta * self.seasonal[tick % 365]

        N          = model.S + model.E + model.I + model.R + model.V
        prevalence = model.I / np.maximum(N, 1.0)

        # Gravity-weighted effective prevalence per patch
        eff_prev = self.coupling @ prevalence

        # Unvaccinated susceptibles: full FOI
        new_E_S  = np.minimum(beta_t * eff_prev * model.S, model.S)
        model.S -= new_E_S

        # Vaccinated: reduced FOI  (leaky vaccine)
        new_E_V  = np.minimum(beta_t * eff_prev * (1 - self.ve_opv) * model.V,
                              model.V)
        model.V -= new_E_V

        model.E += new_E_S + new_E_V


class LatentPhase:
    """E → I  (rate σ = 1/latent period)."""

    def __init__(self, model):
        self.sigma = model.params.sigma

    def __call__(self, model, tick):
        prog     = self.sigma * model.E
        model.E -= prog
        model.I += prog


class RecoveryPhase:
    """I → R  (rate γ = 1/infectious period)."""

    def __init__(self, model):
        self.gamma = model.params.gamma

    def __call__(self, model, tick):
        rec      = self.gamma * model.I
        model.I -= rec
        model.R += rec


class DemographyPhase:
    """Balanced births (→ S) and deaths (from all compartments, incl. V)."""

    def __init__(self, model):
        self.mu = model.params.mu

    def __call__(self, model, tick):
        N = model.S + model.E + model.I + model.R + model.V
        model.S += self.mu * N          # births into S
        model.S -= self.mu * model.S    # deaths from S
        model.E -= self.mu * model.E
        model.I -= self.mu * model.I
        model.R -= self.mu * model.R
        model.V -= self.mu * model.V    # deaths from V


class RecorderPhase:
    """Snapshot every compartment (incl. V) at each tick."""

    def __init__(self, model):
        pass

    def __call__(self, model, tick):
        model.hist["S"][:, tick] = model.S
        model.hist["E"][:, tick] = model.E
        model.hist["I"][:, tick] = model.I
        model.hist["R"][:, tick] = model.R
        model.hist["V"][:, tick] = model.V


# ═══════════════════════════════════════════════════════════════════════════
#  Vaccination components
#
#  These follow the same callable-phase pattern as the SEIR components
#  above.  LASER's ``laser.generic.Immunization`` module provides
#  RoutineImmunization and SIA scaffolds for agent-based models; here we
#  implement the equivalent dynamics for aggregate (patch-level)
#  compartments.
# ═══════════════════════════════════════════════════════════════════════════

class RoutineImmunizationPhase:
    """Routine immunisation (RI) with OPV at 6 weeks of age.

    In the aggregate model, births appear via DemographyPhase as μ·N added
    to S each tick.  RI diverts a fraction ``ri_coverage`` of those daily
    births from S into V on the *same* tick.

    This is the patch-level analogue of LASER's
    ``Immunization.RoutineImmunization`` component.
    """

    def __init__(self, model):
        self.mu       = model.params.mu
        self.coverage = model.params.ri_coverage

    def __call__(self, model, tick):
        N = model.S + model.E + model.I + model.R + model.V
        daily_births = self.mu * N
        vaccinated   = self.coverage * daily_births
        # Move from S → V (births already went into S in DemographyPhase)
        vaccinated   = np.minimum(vaccinated, model.S)
        model.S -= vaccinated
        model.V += vaccinated


class SIAPhase:
    """Supplementary Immunisation Activity (SIA) every 6 months.

    Targets children aged 0-5 years.  In the aggregate model the under-5
    fraction is approximated as ``fraction_under5 ≈ 5 / life_expectancy``.

    On SIA days, ``sia_coverage`` of the estimated susceptible children
    are moved S → V.

    This is the patch-level analogue of LASER's
    ``Immunization.Campaign`` component.
    """

    def __init__(self, model):
        self.coverage  = model.params.sia_coverage
        self.interval  = model.params.sia_interval
        self.frac_u5   = model.params.fraction_under5

    def __call__(self, model, tick):
        if tick > 0 and tick % self.interval == 0:
            target    = self.frac_u5 * model.S     # susceptible children 0-5
            vaccinated = self.coverage * target
            vaccinated = np.minimum(vaccinated, model.S)
            model.S -= vaccinated
            model.V += vaccinated


class WaningPhase:
    """OPV mucosal immunity waning  (V → S).

    OPV-induced immunity wanes at rate ``waning_rate`` (≈ 1/(3 yr)),
    returning individuals to the susceptible pool.  Natural immunity (R)
    does *not* wane, reflecting the stronger and longer-lasting protection
    conferred by wild-virus infection.
    """

    def __init__(self, model):
        self.rate = model.params.waning_rate

    def __call__(self, model, tick):
        waned    = self.rate * model.V
        model.V -= waned
        model.S += waned


# ═══════════════════════════════════════════════════════════════════════════
#  Main – run simulation and plot diagnostics
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":

    # ── Print summary ────────────────────────────────────────────────────
    print("LASER Polio Metapopulation Model (SEIRV)")
    print("=" * 50)
    print(f"  Patches       : {NUM_PATCHES}  (line, {SPACING_KM:.0f} km apart)")
    print(f"  Duration      : {params.nticks} days  ({params.nticks/365:.0f} yr)")
    print(f"  Gravity model : a={params.gravity_a}, b={params.gravity_b}, c={params.gravity_c}")
    print(f"  Max FOI import: {params.max_foi_import*100:.0f} %")
    print(f"  Seasonal amp  : ±{params.seasonal_amplitude}")
    print(f"  Populations   : {pop_sizes.astype(int)}")
    print()
    print("  Vaccination:")
    print(f"    RI coverage    : {params.ri_coverage*100:.0f} % (OPV at 6 wk)")
    print(f"    SIA coverage   : {params.sia_coverage*100:.0f} % (children 0-5)")
    print(f"    SIA interval   : every {params.sia_interval} days")
    print(f"    OPV efficacy   : {params.ve_opv*100:.0f} % (leaky)")
    print(f"    OPV waning     : {1.0/params.waning_rate/365:.1f} yr mean duration")
    print()

    # Verify coupling matrix
    off_diag = coupling.copy()
    np.fill_diagonal(off_diag, 0.0)
    print("Coupling matrix diagnostics:")
    print(f"  Row sums (should be 1.0) : {coupling.sum(axis=1).round(6)}")
    print(f"  Off-diag row sums (≤0.15): {off_diag.sum(axis=1).round(6)}")
    print(f"  Diagonal (local, ≥0.85)  : {coupling.diagonal().round(6)}")
    print()

    # ── Run ──────────────────────────────────────────────────────────────
    model = PolioModel(params, pop_sizes, coupling, seasonal)
    print("Running simulation …")
    model.run()
    print("Done.\n")

    # ── Aggregate results ────────────────────────────────────────────────
    days    = np.arange(params.nticks)
    total_S = model.hist["S"].sum(axis=0)
    total_E = model.hist["E"].sum(axis=0)
    total_I = model.hist["I"].sum(axis=0)
    total_R = model.hist["R"].sum(axis=0)
    total_V = model.hist["V"].sum(axis=0)

    print(f"Peak total infectious: {total_I.max():.0f}  on day {total_I.argmax()}")
    print(f"Final S/E/I/R/V: {total_S[-1]:.0f} / {total_E[-1]:.0f} / "
          f"{total_I[-1]:.0f} / {total_R[-1]:.0f} / {total_V[-1]:.0f}")

    # SIA event summary
    sia_days = [t for t in range(1, params.nticks) if t % params.sia_interval == 0]
    print(f"SIA campaigns conducted: {len(sia_days)}  (days: {sia_days})")

    # ── Six-panel diagnostic plot ─────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # (a) Total SEIRV
    ax = axes[0, 0]
    ax.plot(days, total_S, label="S", color="#1f77b4")
    ax.plot(days, total_E, label="E", color="#ff7f0e")
    ax.plot(days, total_I, label="I", color="#d62728")
    ax.plot(days, total_R, label="R (natural)", color="#2ca02c")
    ax.plot(days, total_V, label="V (OPV)", color="#9467bd", ls="--")
    ax.set_title("(a) Aggregate SEIRV")
    ax.set_xlabel("Day")
    ax.set_ylabel("Population")
    ax.legend()

    # (b) Infectious by patch
    ax = axes[0, 1]
    for i in range(NUM_PATCHES):
        ax.plot(days, model.hist["I"][i], label=f"Patch {i}", alpha=0.7)
    ax.set_title("(b) Infectious by Patch")
    ax.set_xlabel("Day")
    ax.set_ylabel("I")
    ax.legend(fontsize=7, ncol=2)

    # (c) Vaccinated (V) vs Recovered (R)
    ax = axes[0, 2]
    ax.plot(days, total_V, label="V (OPV, waning)", color="#9467bd", lw=2)
    ax.plot(days, total_R, label="R (natural, permanent)", color="#2ca02c", lw=2)
    # Mark SIA days
    for sd in sia_days:
        ax.axvline(sd, color="orange", alpha=0.3, lw=0.8)
    ax.set_title("(c) Vaccinated vs Naturally Immune")
    ax.set_xlabel("Day")
    ax.set_ylabel("Population")
    ax.legend()

    # (d) Seasonal forcing profile
    ax = axes[1, 0]
    doy = np.arange(365)
    ax.fill_between(doy, seasonal, 1.0, alpha=0.3, color="purple")
    ax.plot(doy, seasonal, color="purple", lw=2)
    ax.axhline(1.0, color="grey", ls="--", alpha=0.5)
    ax.set_title("(d) Seasonal Forcing Profile")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel("Transmission multiplier")
    month_labels = ["Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"]
    ax.set_xticks([15,46,74,105,135,166,196,227,258,288,319,349])
    ax.set_xticklabels(month_labels, rotation=45, fontsize=8)

    # (e) Coupling matrix heat-map
    ax = axes[1, 1]
    im = ax.imshow(coupling, cmap="YlOrRd", aspect="auto")
    ax.set_title("(e) Gravity Coupling Matrix (row-normalised)")
    ax.set_xlabel("Source patch j")
    ax.set_ylabel("Receiving patch i")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (f) Effective immunity  (R + ve_opv·V)  vs Susceptible
    ax = axes[1, 2]
    eff_immune = total_R + params.ve_opv * total_V
    ax.fill_between(days, 0, total_S, alpha=0.3, color="#1f77b4",
                    label="Susceptible (S)")
    ax.fill_between(days, 0, eff_immune, alpha=0.3, color="#2ca02c",
                    label=f"Effective immunity (R + {params.ve_opv}·V)")
    ax.set_title("(f) Susceptibility vs Effective Immunity")
    ax.set_xlabel("Day")
    ax.set_ylabel("Population")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("polio_laser_results.png", dpi=150, bbox_inches="tight")
    print("Figure saved to polio_laser_results.png")
