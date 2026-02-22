The workspace filesystem is read-only, so the file is at `/tmp/pakistan_polio_model.py`. Here is the complete, runnable code:

---

```python
#!/usr/bin/env python3
"""
Pakistan Spatial Polio Transmission Model — LASER Framework (v1.0.0)

Specifications:
  - 20 districts, heterogeneous populations (50k-2M)
  - SEIR: R0~6, 3-day latent, 28-day infectious period
  - Gravity-model spatial coupling
  - Monsoon seasonal forcing (peak Jul-Oct, cosine profile)
  - OPV routine immunization at 80% coverage (weekly RI, age 270d)
  - SIA campaigns every 6 months at 80% coverage (ages 0-5yr)
  - Importation: 3 infections every 60 days, years 0-5
  - 20-year simulation, 10-year burn-in
  - Output: weekly incidence per district, proportion of zero-incidence weeks
"""

import sys
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point

from laser.core.propertyset import PropertySet
from laser.core.demographics import AliasedDistribution
import laser.core.distributions as dists
from laser.generic import SEIR, Model
from laser.generic.utils import ValuesMap
from laser.generic.vitaldynamics import BirthsByCBR, MortalityByCDR
from laser.core.migration import gravity, row_normalizer, distance as haversine

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# SECTION 1: Custom Components (fully inline)
# =============================================================================

class PatchImportation:
    """Seeds infections in specific patches every `period` days.

    Targets susceptible-only agents. Models cross-border importation
    from endemic reservoirs (Afghanistan-Pakistan corridor).
    """
    def __init__(self, model, infdurdist, patches,
                 period=60, count=1, end_tick=None):
        self.model      = model
        self.infdurdist = infdurdist
        self.patches    = np.asarray(patches, dtype=np.int32)
        self.period     = period
        self.count      = count
        self.end_tick   = (end_tick if end_tick is not None
                           else model.params.nticks)

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0 or tick >= self.end_tick:
            return
        people   = self.model.people
        nodes    = self.model.nodes
        n_active = people.count
        for pid in self.patches:
            susc = np.nonzero(
                (people.state[:n_active] == SEIR.State.SUSCEPTIBLE.value) &
                (people.nodeid[:n_active] == pid)
            )[0]
            if len(susc) == 0:
                continue
            n      = min(self.count, len(susc))
            chosen = np.random.choice(susc, size=n, replace=False)
            people.state[chosen] = SEIR.State.INFECTIOUS.value
            durations = dists.sample_floats(
                self.infdurdist, np.zeros(n, np.float32))
            durations = np.maximum(np.round(durations), 1).astype(
                people.itimer.dtype)
            people.itimer[chosen] = durations
            nodes.S[tick + 1, pid] = max(0, nodes.S[tick + 1, pid] - n)
            nodes.I[tick + 1, pid] += n


class PerPatchVaccination:
    """OPV routine immunization + biannual SIA with correlated missedness.

    CRITICAL: Sets state = RECOVERED (not susceptibility = 0).
    LASER's TransmissionSE checks state == SUSCEPTIBLE (int8 == 0).
    Setting susceptibility = 0 without changing state has zero effect
    on transmission — the most common silent vaccination bug.

    Correlated Missedness:
        Per-agent `reachable` flag (1=reachable, 0=unreachable) set at
        birth from the district's unreachable_frac. These agents are never
        vaccinated by RI or SIA, modelling persistent zero-dose pockets.
        Prevents independent Bernoulli draws from overestimating cumulative
        coverage across rounds.
    """
    def __init__(self, model, ri_coverage, sia_coverage, unreachable_frac,
                 ri_period=7, ri_age=270, sia_period=180, sia_max_age=5*365):
        self.model            = model
        self.ri_coverage      = np.asarray(ri_coverage,      dtype=np.float32)
        self.sia_coverage     = np.asarray(sia_coverage,     dtype=np.float32)
        self.unreachable_frac = np.asarray(unreachable_frac, dtype=np.float32)
        self.ri_period        = ri_period
        self.ri_age           = ri_age
        self.sia_period       = sia_period
        self.sia_max_age      = sia_max_age
        # add_scalar_property initialises all agents to default=1 (reachable)
        model.people.add_scalar_property("reachable", dtype=np.int8, default=1)
        self._assign(0, model.people.count)

    def _assign(self, istart, iend):
        people = self.model.people
        n = iend - istart
        if n <= 0:
            return
        thresh = self.unreachable_frac[people.nodeid[istart:iend]]
        draws  = np.random.random(n).astype(np.float32)
        people.reachable[istart:iend] = (draws >= thresh).astype(np.int8)

    def on_birth(self, istart, iend, tick):
        self._assign(istart, iend)

    def _vaccinate(self, indices, tick):
        people  = self.model.people
        nodes   = self.model.nodes
        people.state[indices] = SEIR.State.RECOVERED.value
        by_node = np.bincount(
            people.nodeid[indices], minlength=nodes.count
        ).astype(nodes.S.dtype)
        nodes.S[tick + 1] -= by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.R[tick + 1] += by_node

    def step(self, tick):
        people   = self.model.people
        n_active = people.count

        # Routine Immunization: weekly check, 9-month target window
        if tick > 0 and tick % self.ri_period == 0:
            ages = tick - people.dob[:n_active]
            elig = np.nonzero(
                (people.state[:n_active]     == SEIR.State.SUSCEPTIBLE.value) &
                (people.reachable[:n_active]  == 1) &
                (ages >= self.ri_age) &
                (ages <  self.ri_age + self.ri_period)
            )[0]
            if len(elig) > 0:
                cov  = self.ri_coverage[people.nodeid[elig]]
                hits = elig[np.random.random(len(elig)).astype(np.float32) < cov]
                if len(hits) > 0:
                    self._vaccinate(hits, tick)

        # SIA: every 6 months, children 0-5 years
        if tick > 0 and tick % self.sia_period == 0:
            ages = tick - people.dob[:n_active]
            elig = np.nonzero(
                (people.state[:n_active]     == SEIR.State.SUSCEPTIBLE.value) &
                (people.reachable[:n_active]  == 1) &
                (ages >= 0) &
                (ages <  self.sia_max_age)
            )[0]
            if len(elig) > 0:
                cov  = self.sia_coverage[people.nodeid[elig]]
                hits = elig[np.random.random(len(elig)).astype(np.float32) < cov]
                if len(hits) > 0:
                    self._vaccinate(hits, tick)


# =============================================================================
# SECTION 2: District Data
# =============================================================================
# unreachable_frac > 1/R0 (~17%) sustains local transmission chains even
# at high nominal coverage — the structural driver of persistent polio.

DISTRICTS = pd.DataFrame([
    # name,             province,         pop,      lat,     lon,  unreachable
    ("Lahore",         "Punjab",    1_900_000, 31.5204, 74.3587, 0.03),
    ("Faisalabad",     "Punjab",    1_800_000, 31.4504, 73.1350, 0.03),
    ("Rawalpindi",     "Punjab",    1_400_000, 33.5651, 73.0169, 0.03),
    ("Multan",         "Punjab",    1_500_000, 30.1575, 71.5249, 0.04),
    ("Bahawalpur",     "Punjab",      900_000, 29.3544, 71.6911, 0.05),
    ("Islamabad",      "ICT",         450_000, 33.6844, 73.0479, 0.02),
    ("Karachi",        "Sindh",     2_000_000, 24.8607, 67.0011, 0.15),
    ("Hyderabad",      "Sindh",     1_200_000, 25.3960, 68.3578, 0.12),
    ("Sukkur",         "Sindh",       500_000, 27.7244, 68.8228, 0.18),
    ("Jacobabad",      "Sindh",       100_000, 28.2769, 68.4514, 0.25),
    ("Peshawar",       "KP",        1_700_000, 34.0150, 71.5249, 0.18),
    ("Mardan",         "KP",          700_000, 34.1988, 72.0404, 0.15),
    ("Swat",           "KP",          400_000, 35.2220, 72.3528, 0.20),
    ("DI_Khan",        "KP",        1_000_000, 31.8320, 70.9015, 0.40),
    ("Bannu",          "KP",          600_000, 32.9860, 70.6027, 0.40),
    ("N_Waziristan",   "KP",          350_000, 32.3045, 69.8597, 0.55),
    ("S_Waziristan",   "KP",          300_000, 32.0711, 69.5310, 0.55),
    ("Quetta",         "Balochistan",  800_000, 30.1798, 66.9750, 0.40),
    ("Zhob",           "Balochistan",  200_000, 31.3413, 69.4494, 0.45),
    ("Chaman",         "Balochistan",   50_000, 30.9210, 66.4536, 0.55),
], columns=["name", "province", "population", "lat", "lon", "unreachable_frac"])


# =============================================================================
# SECTION 3: Model Parameters
# =============================================================================

PARAMS = PropertySet({
    "prng_seed":              42,
    "nticks":                 20 * 365,
    # Transmission: R0 = beta * D_inf = (6/28) * 28 = 6
    "beta":                   6.0 / 28.0,   # ~0.2143 per day
    # Exposed: Gamma(3, 1) -> mean 3 days latent
    "exp_shape":              3,
    "exp_scale":              1.0,
    # Infectious: Normal(28, 3) -> mean 28 days
    "inf_mean":               28,
    "inf_sigma":              3,
    # Vital dynamics (per-1000/year; LASER divides by 1000 internally)
    "cbr":                    29.0,
    # Gravity (no gravity_a -> no auto-setup; manual network below)
    "gravity_k":              0.005,
    "gravity_b":              0.5,
    "gravity_c":              1.5,
    "capacity_safety_factor": 2.5,
})

CDR          = 7.0    # crude death rate per-1000/year
BURNIN_YEARS = 10
SIM_YEARS    = 20
NTICKS       = SIM_YEARS * 365


# =============================================================================
# SECTION 4: Setup Functions
# =============================================================================

def build_scenario():
    """GeoDataFrame near endemic equilibrium: S/N~17%, R/N~83%."""
    geometry = [Point(lon, lat)
                for lat, lon in zip(DISTRICTS.lat, DISTRICTS.lon)]
    scenario = gpd.GeoDataFrame(DISTRICTS.copy(), geometry=geometry,
                                crs="EPSG:4326")
    scenario["nodeid"] = range(len(scenario))
    pop = scenario["population"]
    scenario["I"] = np.maximum(np.round(0.001  * pop).astype(np.int32), 1)
    scenario["E"] = np.round(0.0003 * pop).astype(np.int32)
    scenario["R"] = np.round(0.83   * pop).astype(np.int32)
    scenario["S"] = (pop - scenario["E"] - scenario["I"]
                     - scenario["R"]).astype(np.int32)
    assert (scenario.S + scenario.E + scenario.I + scenario.R == pop).all()
    assert (scenario.S >= 0).all()
    return scenario


def build_age_pyramid():
    """Pakistan stable age distribution (~65-year life expectancy)."""
    ages    = np.arange(100)
    weights = np.array(1000 * np.exp(-ages / 65.0), dtype=np.int64)
    return AliasedDistribution(np.maximum(weights, 1))


def build_monsoon_seasonality(nticks, nnodes):
    """Cosine forcing peaked during Pakistan monsoon (Jul-Oct).

    season(t) = 1.0 + 0.30 * cos(2*pi*(t-245)/365)
      peak 1.30 at day 245 (early Sep)  |  trough 0.70 at day 62 (Mar)
      mean = 1.0 (cosine over full period integrates to zero)
    """
    days    = np.arange(365)
    profile = 1.0 + 0.30 * np.cos(2 * np.pi * (days - 245) / 365)
    assert abs(profile.mean() - 1.0) < 0.001
    tiled      = np.tile(profile, nticks // 365 + 1)[:nticks]
    seasonality = ValuesMap.from_timeseries(tiled, nnodes)
    return seasonality, profile


def initialize_ages(model, pyramid):
    """Assign realistic initial ages to prevent all-newborn cohort artifact.

    Without this, all agents have dob=0. At tick 270 every agent hits the RI
    age window simultaneously, instantly depleting most susceptibles.
    """
    count     = model.people.count
    ages_yrs  = pyramid.sample(count=count, dtype=np.int32)
    ages_yrs  = np.minimum(ages_yrs, 99)
    ages_days = ages_yrs * 365 + np.random.randint(0, 365, size=count)
    model.people.dob[:count] = -ages_days.astype(model.people.dob.dtype)


def build_distance_matrix(lats, lons):
    """Pairwise Haversine great-circle distance matrix (km)."""
    n = len(lats)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = haversine(lats[i], lons[i], lats[j], lons[j])
    return D


# =============================================================================
# SECTION 5: Build and Run the Model
# =============================================================================

def run_model():
    scenario = build_scenario()
    nnodes   = len(scenario)

    expdurdist = dists.gamma(shape=PARAMS.exp_shape, scale=PARAMS.exp_scale)
    infdurdist = dists.normal(loc=PARAMS.inf_mean,   scale=PARAMS.inf_sigma)

    seasonality, season_365 = build_monsoon_seasonality(NTICKS, nnodes)

    # Vital dynamics rates — UNIT GUARD (wrong units = silent birth failure)
    birthrates = np.full((NTICKS, nnodes), PARAMS.cbr, dtype=np.float32)
    deathrates = np.full((NTICKS, nnodes), CDR,         dtype=np.float32)
    assert np.all(birthrates >= 1) and np.all(birthrates <= 60), \
        "CBR must be per-1000/yr (10-50 typical)"

    pyramid = build_age_pyramid()
    model   = Model(scenario, PARAMS, birthrates=birthrates)

    # Gravity network: M_{i,j} = k * p_j^b / d_{ij}^c  (a=0 by convention)
    lats       = np.array(scenario.lat)
    lons       = np.array(scenario.lon)
    pops       = np.array(scenario.population, dtype=np.float64)
    D          = build_distance_matrix(lats, lons)
    network    = gravity(pops, D, 1.0, 0, PARAMS.gravity_b, PARAMS.gravity_c)
    avg_exp    = network.sum(axis=1).mean()
    if avg_exp > 0:
        network = network / avg_exp * PARAMS.gravity_k
    network       = row_normalizer(network, 0.2)
    model.network = network

    assert model.network.sum() > 0, "Network is all-zero"
    assert model.network.sum(axis=1).max() < 0.3

    ri_cov      = np.full(nnodes, 0.80, dtype=np.float32)
    sia_cov     = np.full(nnodes, 0.80, dtype=np.float32)
    unreachable = np.array(scenario.unreachable_frac, dtype=np.float32)

    # 3 patches x 1 infection/event = 3 total per 60-day importation cycle
    n2i            = {n: i for i, n in enumerate(scenario.name)}
    import_patches = [n2i["Quetta"], n2i["Peshawar"], n2i["N_Waziristan"]]

    # Component ordering: Susceptible/Recovered bracket transitions;
    # Transmission runs last on current-tick counts.
    model.components = [
        SEIR.Susceptible(model),
        SEIR.Exposed(model, expdurdist, infdurdist),
        SEIR.Infectious(model, infdurdist),
        SEIR.Recovered(model),
        BirthsByCBR(model, birthrates=birthrates, pyramid=pyramid),
        MortalityByCDR(model, mortalityrates=deathrates),
        PerPatchVaccination(model, ri_cov, sia_cov, unreachable,
                            ri_period=7, ri_age=270, sia_period=180),
        PatchImportation(model, infdurdist, import_patches,
                         period=60, count=1, end_tick=5 * 365),
        SEIR.Transmission(model, expdurdist, seasonality=seasonality),
    ]

    initialize_ages(model, pyramid)   # Must run BEFORE model.run()

    R0 = PARAMS.beta * PARAMS.inf_mean
    print(f"\nPakistan Polio 20-District SEIR Model")
    print(f"  {nnodes} districts  |  pop={int(pops.sum()):,}")
    print(f"  R0={R0:.2f}  beta={PARAMS.beta:.4f}/d  "
          f"D_lat={PARAMS.exp_shape:.0f}d  D_inf={PARAMS.inf_mean}d")
    print(f"  Seasonality: monsoon peak day 245, +/-30%")
    print(f"  RI 80% (age 270d, weekly)  |  SIA 80% (every 180d, 0-5yr)")
    print(f"  Unreachable: {np.average(unreachable, weights=pops):.1%} pop-wtd avg")
    print(f"  Importation: 3/60d in Quetta, Peshawar, N_Waziristan (yr 0-5)")
    print(f"  {SIM_YEARS}yr sim ({BURNIN_YEARS}yr burn-in)  |  {NTICKS} ticks")

    model.run("Pakistan Polio 20-District SEIR")
    return model, scenario, season_365, D, network


# =============================================================================
# SECTION 6: Post-Run Verification
# =============================================================================

def verify_simulation(model):
    """Five health checks to catch silent failures after model.run()."""
    nticks = model.params.nticks
    nodes  = model.nodes
    failures = []
    print("\n=== Post-Run Verification ===")

    def pop_at(t):
        return int(nodes.S[t].sum() + nodes.E[t].sum()
                   + nodes.I[t].sum() + nodes.R[t].sum())

    # 1. Population growth (CBR=29 > CDR=7)
    ratio   = pop_at(nticks - 1) / max(pop_at(1), 1)
    growing = ratio > 1.01
    print(f"  [{'OK  ' if growing else 'WARN'}] Population ratio={ratio:.3f} "
          f"({'growing' if growing else 'STATIC — check CBR units'})")
    if not growing:
        failures.append("population")

    # 2. Non-negativity
    st = np.linspace(1, nticks - 1, 10, dtype=int)
    mins = {c: int(getattr(nodes, c)[st].min()) for c in ["S", "E", "I", "R"]}
    ok   = all(v >= 0 for v in mins.values())
    print(f"  [{'OK  ' if ok else 'FAIL'}] Non-negativity: "
          f"min S/E/I/R = {mins['S']}/{mins['E']}/{mins['I']}/{mins['R']}")
    if not ok:
        failures.append("nonnegativity")

    # 3. Epidemic dynamics
    total    = int(nodes.newly_infected[:nticks].sum())
    nw       = nticks // 7
    wkly     = (nodes.newly_infected[:nw * 7]
                .reshape(nw, 7, nodes.count).sum(axis=(1, 2)))
    act_frac = float(np.count_nonzero(wkly)) / nw
    late_act = int(np.count_nonzero(wkly[3 * nw // 4:]))
    ok       = total > 0 and act_frac > 0.05 and late_act > 0
    print(f"  [{'OK  ' if ok else 'FAIL'}] Epidemic: "
          f"total={total:,}, active={act_frac:.0%} weeks, late={late_act} wks")
    if not ok:
        failures.append("epidemic")

    # 4. Vaccination effect
    def sfrac(t):
        s = nodes.S[t].sum()
        n = s + nodes.E[t].sum() + nodes.I[t].sum() + nodes.R[t].sum()
        return s / max(n, 1)
    early = float(np.mean([sfrac(t)
                            for t in np.linspace(1, nticks//4, 5, dtype=int)]))
    late  = float(np.mean([sfrac(t)
                            for t in np.linspace(3*nticks//4, nticks-1, 5, dtype=int)]))
    ok    = abs(early - late) > 0.001 or late < 0.40
    print(f"  [{'OK  ' if ok else 'WARN'}] S/N: {early:.3f} -> {late:.3f}")

    # 5. Network
    ok = float(model.network.sum()) > 0
    print(f"  [{'OK  ' if ok else 'FAIL'}] Network: "
          f"sum={model.network.sum():.4f}, "
          f"max_row={model.network.sum(axis=1).max():.4f}")
    if not ok:
        failures.append("network")

    if failures:
        print(f"\n  FAILED: {failures}")
        if set(failures) & {"nonnegativity", "epidemic", "network"}:
            raise RuntimeError(f"Critical failures: {failures}")
    else:
        print("\n  All checks passed.")
    print("=" * 30)


# =============================================================================
# SECTION 7: Output Analysis
# =============================================================================

def compute_weekly_incidence(model, scenario):
    """Post-burn-in weekly incidence per district.

    Returns
    -------
    df  : DataFrame (week, year, district, province, population, incidence)
    arr : ndarray (n_weeks, nnodes)
    """
    nticks   = model.params.nticks
    nnodes   = len(scenario)
    post_inc = model.nodes.newly_infected[BURNIN_YEARS * 365:nticks]
    n_weeks  = post_inc.shape[0] // 7
    arr      = post_inc[:n_weeks * 7].reshape(n_weeks, 7, nnodes).sum(axis=1)
    rows = []
    for w in range(n_weeks):
        for d in range(nnodes):
            rows.append({
                "week":       w + 1,
                "year":       round(BURNIN_YEARS + w * 7 / 365.0, 3),
                "district":   scenario.name.iloc[d],
                "province":   scenario.province.iloc[d],
                "population": int(scenario.population.iloc[d]),
                "incidence":  int(arr[w, d]),
            })
    return pd.DataFrame(rows), arr


def analyze_zero_incidence(weekly_df, scenario):
    """Print per-district proportion of zero-incidence weeks."""
    n_weeks = weekly_df["week"].nunique()
    nnodes  = len(scenario)
    print(f"\n{'='*76}")
    print(f"ZERO-INCIDENCE WEEK ANALYSIS  ({n_weeks} post-burn-in weeks)")
    print(f"{'='*76}")
    print(f"\n{'District':<18} {'Province':<14} {'Pop':>10}  "
          f"{'Total':>8}  {'Active':>7}  {'Silent':>7}  {'%Silent':>8}")
    print("-" * 76)
    zf = np.zeros(nnodes)
    for d in range(nnodes):
        name   = scenario.name.iloc[d]
        data   = weekly_df[weekly_df["district"] == name]
        total  = int(data["incidence"].sum())
        n_act  = int((data["incidence"] > 0).sum())
        n_zero = n_weeks - n_act
        pct    = 100.0 * n_zero / n_weeks
        zf[d]  = pct / 100.0
        print(f"{name:<18} {scenario.province.iloc[d]:<14} "
              f"{scenario.population.iloc[d]:>10,}  "
              f"{total:>8,}  {n_act:>7}  {n_zero:>7}  {pct:>7.1f}%")
    print("-" * 76)
    print(f"{'OVERALL (unweighted avg)':<52}  "
          f"{'':{7}}  {'':{7}}  {zf.mean():>7.1%}")
    print(f"{'OVERALL (pop-weighted avg)':<52}  "
          f"{'':{7}}  {'':{7}}  "
          f"{np.average(zf, weights=scenario.population.values):>7.1%}")
    return zf


def print_summary(model, scenario):
    """Per-district annual incidence and compartment integrity."""
    nticks = model.params.nticks
    nyears = nticks // 365
    inc    = model.nodes.newly_infected
    pops   = np.array(scenario.population, dtype=np.float64)
    R0     = PARAMS.beta * PARAMS.inf_mean
    nnodes = len(scenario)
    print(f"\n{'='*76}")
    print(f"SIMULATION SUMMARY  (post burn-in, years {BURNIN_YEARS}-{SIM_YEARS})")
    print(f"{'='*76}")
    print(f"\n{'District':<18} {'Pop':>10}  {'Med/yr':>9}  "
          f"{'Avg/yr':>9}  {'S/N(final)':>11}")
    print("-" * 65)
    total_post = 0
    for i in range(nnodes):
        annual = np.array([inc[y*365:(y+1)*365, i].sum()
                           for y in range(BURNIN_YEARS, nyears)])
        total_post += int(annual.sum())
        t = nticks - 1
        N = max(int(model.nodes.S[t, i] + model.nodes.E[t, i]
                    + model.nodes.I[t, i] + model.nodes.R[t, i]), 1)
        print(f"{scenario.name.iloc[i]:<18} {pops[i]:>10,.0f}  "
              f"{np.median(annual):>9,.0f}  {annual.mean():>9,.0f}  "
              f"{model.nodes.S[t, i] / N:>11.3f}")
    ayrs  = nyears - BURNIN_YEARS
    total = int(pops.sum())
    print(f"\n  Total pop: {total:,}  |  Post-burn-in infections: {total_post:,} "
          f"over {ayrs}yr")
    print(f"  Avg annual: {total_post/ayrs:,.0f}  |  "
          f"Per million/yr: {total_post/ayrs/total*1e6:.0f}  |  R0={R0:.1f}")
    print(f"\nCompartment integrity (S+E+I+R=N):")
    for label, t in [("Start", 1), ("Mid", nticks//2), ("End", nticks-1)]:
        S, E, I, R = (int(model.nodes.S[t].sum()), int(model.nodes.E[t].sum()),
                      int(model.nodes.I[t].sum()), int(model.nodes.R[t].sum()))
        print(f"  {label:<5} t={t:>5}: S={S:>10,} E={E:>6,} "
              f"I={I:>5,} R={R:>10,} N={S+E+I+R:>10,}")
    print(f"\nPopulation trajectory:")
    for y in range(0, nyears+1, 5):
        t = min(y * 365, nticks - 1)
        N = int(model.nodes.S[t].sum() + model.nodes.E[t].sum()
                + model.nodes.I[t].sum() + model.nodes.R[t].sum())
        print(f"  Year {y:>2}: {N:>12,}")


# =============================================================================
# SECTION 8: Diagnostic Plots
# =============================================================================

def plot_diagnostics(model, scenario, season_365, weekly_arr, zf, outdir):
    """Six-panel diagnostic figure.

    (A) Monsoon seasonal forcing    (B) Annual incidence per million
    (C) R_eff = R0*S/N over time   (D) Annual incidence heatmap
    (E) Weekly incidence time series (F) Zero-incidence bar chart
    """
    nticks = model.params.nticks
    nnodes = len(scenario)
    nyears = nticks // 365
    names  = list(scenario.name)
    pops   = np.array(scenario.population, dtype=np.float64)
    R0     = PARAMS.beta * PARAMS.inf_mean
    inc    = model.nodes.newly_infected

    prov_colors = {
        "Punjab": "#2ecc71", "ICT": "#9b59b6",
        "Sindh":  "#3498db", "KP": "#e67e22",
        "Balochistan": "#e74c3c",
    }
    dc = [prov_colors.get(p, "gray") for p in scenario.province]

    annual_inc = np.array([inc[y*365:(y+1)*365].sum(axis=0)
                           for y in range(nyears)])
    s30  = np.arange(0, nticks, 30)
    t30  = s30 / 365.0
    N30  = np.maximum(
        model.nodes.S[s30].astype(np.float64) +
        model.nodes.E[s30] + model.nodes.I[s30] + model.nodes.R[s30], 1.0)
    Reff = R0 * model.nodes.S[s30].astype(np.float64) / N30

    fig, axes = plt.subplots(2, 3, figsize=(21, 12))

    # (A) Seasonal forcing
    ax = axes[0, 0]
    ax.plot(np.arange(365), season_365, "b-", lw=2)
    ax.axhline(1.0, color="gray", ls="--", alpha=0.6)
    ax.axvspan(182, 304, alpha=0.15, color="steelblue", label="Monsoon Jul-Oct")
    ax.axvspan(335, 365, alpha=0.10, color="orange")
    ax.axvspan(0,   89, alpha=0.10, color="orange", label="Dry Dec-Mar")
    ax.set(xlabel="Day of Year", ylabel="Seasonal Multiplier",
           title="(A) Monsoon Seasonal Forcing", ylim=(0.60, 1.44))
    ax.legend(fontsize=9)

    # (B) Annual incidence rate
    ax = axes[0, 1]
    for i in range(nnodes):
        rate = annual_inc[:, i] / pops[i] * 1e6
        if rate.max() >= 1:
            ax.plot(np.arange(nyears) + 0.5, rate, color=dc[i],
                    alpha=0.75, lw=1, label=names[i])
    ax.axvline(BURNIN_YEARS, color="black", ls=":", lw=1.5, alpha=0.5)
    ax.set_yscale("symlog", linthresh=10)
    ax.set(xlabel="Year", ylabel="Annual Incidence per Million",
           title="(B) Incidence Rate by District")
    ax.legend(fontsize=4.5, ncol=3, loc="upper right")

    # (C) R_eff
    ax = axes[0, 2]
    for i in range(nnodes):
        ax.plot(t30, Reff[:, i], color=dc[i], alpha=0.6, lw=0.8)
    ax.axhline(1.0, color="black", ls="-", lw=2, alpha=0.35, label="R_eff=1")
    ax.axvline(BURNIN_YEARS, color="gray", ls=":", alpha=0.5)
    ax.set(xlabel="Year", ylabel="R_eff = R0 x S/N",
           title="(C) Effective Reproduction Number",
           ylim=(0, min(5, float(Reff.max()) * 1.1 + 0.5)))
    ax.legend(fontsize=9)

    # (D) Heatmap
    ax = axes[1, 0]
    log_rate = np.log10(annual_inc / pops[None, :] * 1e6 + 1)
    im = ax.imshow(log_rate.T, aspect="auto", cmap="YlOrRd", origin="lower",
                   extent=[0, nyears, -0.5, nnodes - 0.5])
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=6.5)
    ax.axvline(BURNIN_YEARS, color="white", ls="--", lw=1.5, alpha=0.85)
    ax.set(xlabel="Year", title="(D) Annual Incidence per Million (log10+1)")
    fig.colorbar(im, ax=ax, shrink=0.85, label="log10(cases/M + 1)")

    # (E) Weekly post-burn-in
    ax      = axes[1, 1]
    n_weeks = weekly_arr.shape[0]
    wx      = np.arange(n_weeks) / 52.177 + BURNIN_YEARS
    for i in range(nnodes):
        if weekly_arr[:, i].max() > 0:
            ax.plot(wx, weekly_arr[:, i], color=dc[i],
                    alpha=0.7, lw=0.65, label=names[i])
    for yr in range(BURNIN_YEARS, SIM_YEARS):
        ax.axvspan(yr + 182/365, yr + 304/365, alpha=0.04, color="steelblue")
    ax.set(xlabel="Year", ylabel="Weekly Cases",
           title="(E) Post-Burn-in Weekly Incidence  (blue=monsoon)")
    ax.legend(fontsize=4.5, ncol=3, loc="upper right")

    # (F) Zero-incidence proportions
    ax = axes[1, 2]
    ax.barh(range(nnodes), zf * 100, color=[dc[i] for i in range(nnodes)],
            edgecolor="white", linewidth=0.3)
    ax.set_yticks(range(nnodes))
    ax.set_yticklabels(names, fontsize=6.5)
    ax.axvline((1 - 1/R0) * 100, color="gray", ls="--", lw=1, alpha=0.6,
               label=f"HIT = {1 - 1/R0:.0%}")
    ax.set(xlabel="% Weeks with Zero Incidence",
           title="(F) Zero-Incidence Weeks  (post burn-in)", xlim=(0, 108))
    for i, v in enumerate(zf * 100):
        ax.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=6)
    ax.legend(fontsize=8)

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, label=p) for p, c in prov_colors.items()]
    fig.legend(handles=handles, title="Province", loc="lower center",
               ncol=5, fontsize=8, bbox_to_anchor=(0.5, -0.01))

    plt.suptitle(
        f"Pakistan Polio 20-District SEIR  |  "
        f"R0={R0:.0f}, {SIM_YEARS}yr sim ({BURNIN_YEARS}yr burn-in)\n"
        f"Gravity k={PARAMS.gravity_k}, b={PARAMS.gravity_b}, "
        f"c={PARAMS.gravity_c}  |  RI+SIA 80%  |  "
        f"3 importations/60d  (years 0-5)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    outpath = outdir / "pakistan_polio_diagnostics.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nDiagnostic plot saved to {outpath}")
    plt.close(fig)


# =============================================================================
# SECTION 9: Main
# =============================================================================

if __name__ == "__main__":

    model, scenario, season_365, D, network = run_model()
    verify_simulation(model)
    print_summary(model, scenario)

    weekly_df, weekly_arr = compute_weekly_incidence(model, scenario)

    outdir = Path("eval/outputs")
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "pakistan_polio_weekly_incidence.csv"
    weekly_df.to_csv(csv_path, index=False)
    n_weeks = weekly_df["week"].nunique()
    print(f"\nWeekly incidence CSV: {csv_path}")
    print(f"  {len(weekly_df):,} rows  ({n_weeks} weeks x {len(scenario)} districts)")
    print(f"  Columns: {list(weekly_df.columns)}")

    zf = analyze_zero_incidence(weekly_df, scenario)
    plot_diagnostics(model, scenario, season_365, weekly_arr, zf, outdir)
    print("\nDone.")
```

---

**Key design decisions and LASER-specific details:**

**Beta calculation:**  `R0 / D_inf = 6/28 ≈ 0.2143/day`. Each tick is one day; `TransmissionSE` computes `FOI = beta * I/N` then converts rate→probability via `-expm1(-FOI)`.

**Exposed/Infectious distributions:**  `Gamma(3, 1)` → mean 3-day latent; `Normal(28, 3)` → mean 28-day infectious. These parameterize LASER's stochastic E→I and I→R timers.

**Vaccination correctness:**  `PerPatchVaccination` sets `state = RECOVERED` (not `susceptibility = 0`). LASER's `TransmissionSE` only checks `state == SUSCEPTIBLE`; the susceptibility field is irrelevant for transmission.

**Importation:**  3 patches × 1 case/event every 60 days = 3 total. Stops at tick 5×365. After year 5, endemic persistence depends on endogenous transmission in high-unreachable districts (N_Waziristan, S_Waziristan, Chaman with 55% unreachable > 1/R0 ≈ 17%).

**Age initialization:**  `initialize_ages()` must be called before `model.run()`. Without it, all 20M agents have `dob=0` and trigger a single massive RI event at tick 270.

**Birthrate units:**  Must be `per-1000/year` (here: 29.0). Wrong units (daily per-capita ≈ 0.00008) cause `calc_capacity` to see near-zero growth → no free slots → zero births silently.
