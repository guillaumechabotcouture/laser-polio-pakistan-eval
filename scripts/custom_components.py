#!/usr/bin/env python3
"""
Custom LASER components for Pakistan polio SEIRV model.

Components:
    VaccinatedCompartment - Tracks vaccinated (V) agents separately from naturally
                            recovered (R), with waning OPV immunity (V→S).

    PerPatchVaccinationSEIRV - Routine immunization at 6 weeks + SIA every 6 months.
                                Sets state=VACCINATED (value 4) to distinguish from
                                natural recovery (state=RECOVERED, value 3).
                                Includes correlated missedness via per-agent reachable flag.

    PerPatchVaccination - (Legacy) Original vaccination setting state=RECOVERED.

    PulseImportation    - Introduces infections into patch 0 at regular intervals
                          (default: 5 infections every 90 ticks).

    PatchImportation    - Seeds infections in endemic corridor patches only
                          to represent cross-border and persistent transmission foci.

Usage:
    from custom_components import (VaccinatedCompartment, PerPatchVaccinationSEIRV,
                                   PatchImportation)
"""

import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR

# New state value for vaccinated agents (distinct from RECOVERED=3).
# TransmissionSE only infects agents with state==SUSCEPTIBLE(0), so any
# non-zero state value prevents infection.
VACCINATED = np.int8(4)


class VaccinatedCompartment:
    """Tracks vaccinated agents with waning OPV immunity.

    OPV provides mucosal immunity that wanes over time, distinct from the
    stronger, longer-lasting immunity conferred by natural poliovirus infection.
    When vaccine immunity wanes (vtimer reaches 0), agents return to SUSCEPTIBLE.

    Requires Model to be constructed with additional_states=["V"] so that
    V[t+1] is automatically initialized from V[t] each tick.

    Parameters:
        model: LASER Model instance (must have additional_states=["V"])
        wandurdist: Distribution for waning duration sampling (e.g., gamma with
                    mean ~3 years for OPV type 1 immunity)
        wandurmin: Minimum waning duration in days (default: 365)
    """

    def __init__(self, model, wandurdist, wandurmin=365):
        self.model = model
        self.wandurdist = wandurdist
        self.wandurmin = wandurmin

        # Per-agent vaccine immunity countdown timer
        model.people.add_scalar_property("vtimer", dtype=np.int16, default=0)

        # Node-level V count — created by Model(additional_states=["V"])
        # which also registers V in model.states so _initialize_flows copies
        # V[t+1]=V[t] each tick. Only create V here if Model didn't.
        try:
            _ = model.nodes.V
        except AttributeError:
            model.nodes.add_vector_property("V", model.params.nticks + 1, dtype=np.int32)

        model.nodes.add_vector_property("newly_vax_waned", model.params.nticks + 1, dtype=np.int32)

    def sample_waning_duration(self, n):
        """Sample n waning durations from the distribution."""
        samples = dists.sample_floats(self.wandurdist, np.zeros(n, dtype=np.float32))
        samples = np.maximum(np.round(samples), self.wandurmin).astype(np.int16)
        return samples

    def step(self, tick):
        """Decrement vtimer for vaccinated agents; waned agents return to S."""
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # Find vaccinated agents with active timers
        vaccinated_mask = (people.state[:count] == VACCINATED)
        active_timer = vaccinated_mask & (people.vtimer[:count] > 0)
        active_idx = np.nonzero(active_timer)[0]

        if len(active_idx) == 0:
            return

        # Decrement timers
        people.vtimer[active_idx] -= 1

        # Agents whose timer just hit 0: wane back to susceptible
        waned_idx = active_idx[people.vtimer[active_idx] == 0]

        if len(waned_idx) == 0:
            return

        people.state[waned_idx] = SEIR.State.SUSCEPTIBLE.value

        waned_by_node = np.bincount(
            people.nodeid[waned_idx], minlength=nodes.count
        ).astype(nodes.V.dtype)

        nodes.V[tick + 1] -= waned_by_node
        np.maximum(nodes.V[tick + 1], 0, out=nodes.V[tick + 1])
        nodes.S[tick + 1] += waned_by_node
        nodes.newly_vax_waned[tick] = waned_by_node


class PerPatchVaccinationSEIRV:
    """Per-district RI + SIA vaccination for the SEIRV model.

    Sets state = VACCINATED (value 4) to track vaccine-derived immunity
    separately from natural infection immunity (state = RECOVERED, value 3).

    Routine Immunization (RI):
        - Administered at age 6 weeks (42 days) with configurable coverage
        - Only susceptible, reachable agents are eligible
        - Checked weekly (every ri_period days)

    Supplementary Immunization Activities (SIA):
        - Campaigns every sia_period days (default: 180 = biannual)
        - Targets all children aged 0-5 years (0-1825 days)
        - Susceptible agents → VACCINATED (new vaccination)
        - Already-vaccinated agents get vtimer reset (booster effect)
        - Only reachable agents are eligible

    Correlated missedness:
        Each agent has a permanent `reachable` flag (0=unreachable, 1=reachable)
        set at birth based on the district's unreachable fraction. Unreachable
        agents are never vaccinated by either RI or SIA, modeling the real-world
        dynamic where the same zero-dose children are missed by every campaign.

    Parameters:
        model: LASER Model instance
        vaccinated_compartment: VaccinatedCompartment instance (for waning timer sampling)
        ri_coverage: float or 1D array — RI coverage among reachable (default: 0.80)
        sia_coverage: float or 1D array — SIA coverage among reachable (default: 0.90)
        unreachable_frac: 1D array of per-patch unreachable fraction (0-1)
        ri_period: Days between RI checks (default: 7)
        ri_age: Age in days at which RI is given (default: 42 = 6 weeks)
        sia_period: Days between SIA campaigns (default: 180 = biannual)
        sia_max_age: Max age targeted by SIA in days (default: 1825 = 5 years)
    """

    def __init__(self, model, vaccinated_compartment, unreachable_frac,
                 ri_coverage=0.80, sia_coverage=0.90,
                 ri_period=7, ri_age=42, sia_period=180, sia_max_age=1825):
        self.model = model
        self.vax_compartment = vaccinated_compartment
        self.ri_period = ri_period
        self.ri_age = ri_age
        self.sia_period = sia_period
        self.sia_max_age = sia_max_age

        nnodes = model.nodes.count
        # Broadcast scalar coverage to per-patch arrays
        if np.isscalar(ri_coverage):
            self.ri_coverage = np.full(nnodes, ri_coverage, dtype=np.float32)
        else:
            self.ri_coverage = np.asarray(ri_coverage, dtype=np.float32)
        if np.isscalar(sia_coverage):
            self.sia_coverage = np.full(nnodes, sia_coverage, dtype=np.float32)
        else:
            self.sia_coverage = np.asarray(sia_coverage, dtype=np.float32)

        self.unreachable_frac = np.asarray(unreachable_frac, dtype=np.float32)

        # Add reachable property (0=unreachable, 1=reachable) if not already present
        if not hasattr(model.people, "reachable"):
            model.people.add_property("reachable", dtype=np.int8)

        # Initialize reachability for all existing agents
        self._set_reachability(0, model.people.count)

        # Track doses administered for monitoring
        model.nodes.add_vector_property("ri_doses_given", model.params.nticks + 1, dtype=np.int32)
        model.nodes.add_vector_property("sia_doses_given", model.params.nticks + 1, dtype=np.int32)

    def _set_reachability(self, istart, iend):
        """Assign reachable flag based on district-level unreachable fraction."""
        people = self.model.people
        n = iend - istart
        if n <= 0:
            return
        nodeids = people.nodeid[istart:iend]
        thresholds = self.unreachable_frac[nodeids]
        draws = np.random.random(n).astype(np.float32)
        people.reachable[istart:iend] = (draws >= thresholds).astype(np.int8)

    def on_birth(self, istart, iend, tick):
        """Set reachability for newborn agents (called by BirthsByCBR)."""
        self._set_reachability(istart, iend)

    def _vaccinate_susceptible(self, indices, tick):
        """Move susceptible agents to VACCINATED state and assign waning timer."""
        people = self.model.people
        nodes = self.model.nodes

        people.state[indices] = VACCINATED
        people.vtimer[indices] = self.vax_compartment.sample_waning_duration(len(indices))

        vax_by_node = np.bincount(
            people.nodeid[indices], minlength=nodes.count
        ).astype(nodes.S.dtype)
        nodes.S[tick + 1] -= vax_by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.V[tick + 1] += vax_by_node

    def _boost_vaccinated(self, indices):
        """Reset waning timer for already-vaccinated agents (booster dose)."""
        self.model.people.vtimer[indices] = self.vax_compartment.sample_waning_duration(
            len(indices)
        )

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # --- Routine Immunization (weekly check) ---
        if tick > 0 and tick % self.ri_period == 0:
            ages = tick - people.dob[:count]
            ri_eligible = np.nonzero(
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
                (people.reachable[:count] == 1) &
                (ages >= self.ri_age) &
                (ages < self.ri_age + self.ri_period)
            )[0]

            if len(ri_eligible) > 0:
                agent_coverage = self.ri_coverage[people.nodeid[ri_eligible]]
                draws = np.random.random(len(ri_eligible)).astype(np.float32)
                vaccinated = ri_eligible[draws < agent_coverage]

                if len(vaccinated) > 0:
                    self._vaccinate_susceptible(vaccinated, tick)
                    ri_by_node = np.bincount(
                        people.nodeid[vaccinated], minlength=nodes.count
                    ).astype(nodes.ri_doses_given.dtype)
                    nodes.ri_doses_given[tick] += ri_by_node

        # --- SIA Campaigns (every sia_period days) ---
        if tick > 0 and tick % self.sia_period == 0:
            ages = tick - people.dob[:count]
            age_eligible = (
                (people.reachable[:count] == 1) &
                (ages >= 0) &
                (ages < self.sia_max_age)
            )

            # Susceptible children → new vaccination
            susc_eligible = np.nonzero(
                age_eligible &
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value)
            )[0]

            if len(susc_eligible) > 0:
                agent_coverage = self.sia_coverage[people.nodeid[susc_eligible]]
                draws = np.random.random(len(susc_eligible)).astype(np.float32)
                newly_vaxed = susc_eligible[draws < agent_coverage]
                if len(newly_vaxed) > 0:
                    self._vaccinate_susceptible(newly_vaxed, tick)
                    sia_by_node = np.bincount(
                        people.nodeid[newly_vaxed], minlength=nodes.count
                    ).astype(nodes.sia_doses_given.dtype)
                    nodes.sia_doses_given[tick] += sia_by_node

            # Already-vaccinated children → booster (reset waning timer)
            vax_eligible = np.nonzero(
                age_eligible &
                (people.state[:count] == VACCINATED)
            )[0]

            if len(vax_eligible) > 0:
                agent_coverage = self.sia_coverage[people.nodeid[vax_eligible]]
                draws = np.random.random(len(vax_eligible)).astype(np.float32)
                boosted = vax_eligible[draws < agent_coverage]
                if len(boosted) > 0:
                    self._boost_vaccinated(boosted)
                    boost_by_node = np.bincount(
                        people.nodeid[boosted], minlength=nodes.count
                    ).astype(nodes.sia_doses_given.dtype)
                    nodes.sia_doses_given[tick] += boost_by_node


class PerPatchVaccination:
    """Per-district routine immunization and SIA with correlated missedness.

    Models two populations per district:
        - REACHABLE agents: can be vaccinated by both RI and SIA
        - UNREACHABLE agents: never receive RI or SIA (zero-dose children
          in hard-to-reach/refusal communities)

    This models the key real-world dynamic where the same children are
    consistently missed by every campaign, rather than independent Bernoulli
    draws which overestimate cumulative SIA effectiveness.

    The `reachable` property (int8: 0=unreachable, 1=reachable) is set at
    agent creation and persists for life. Newborns inherit reachability via
    on_birth() callback from BirthsByCBR.

    CRITICAL: Sets state = RECOVERED (value 3) rather than susceptibility = 0,
    because the built-in TransmissionSE only checks state == SUSCEPTIBLE (0)
    when selecting agents for infection.

    Parameters:
        model: LASER Model instance
        ri_coverage: 1D array of per-patch RI coverage rates AMONG REACHABLE
        sia_coverage: 1D array of per-patch SIA coverage rates AMONG REACHABLE
        unreachable_frac: 1D array of per-patch unreachable fraction (0-1)
        ri_period: Days between RI checks (default: 7)
        ri_age: Age in days at which RI is administered (default: 270 = 9 months)
        sia_period: Days between SIA campaigns (default: 180 = ~2/year)
        sia_max_age: Maximum age targeted by SIA in days (default: 1825 = 5 years)
    """

    def __init__(self, model, ri_coverage, sia_coverage, unreachable_frac,
                 ri_period=7, ri_age=270, sia_period=180, sia_max_age=1825):
        self.model = model
        self.ri_coverage = np.asarray(ri_coverage, dtype=np.float32)
        self.sia_coverage = np.asarray(sia_coverage, dtype=np.float32)
        self.unreachable_frac = np.asarray(unreachable_frac, dtype=np.float32)
        self.ri_period = ri_period
        self.ri_age = ri_age
        self.sia_period = sia_period
        self.sia_max_age = sia_max_age

        # Add reachable property to people (0=unreachable, 1=reachable)
        model.people.add_property("reachable", dtype=np.int8)

        # Initialize reachability for all existing agents
        self._set_reachability(0, model.people.count)

    def _set_reachability(self, istart, iend):
        """Assign reachable flag based on district-level unreachable fraction."""
        people = self.model.people
        n = iend - istart
        if n <= 0:
            return
        nodeids = people.nodeid[istart:iend]
        thresholds = self.unreachable_frac[nodeids]
        draws = np.random.random(n).astype(np.float32)
        # reachable=1 if draw >= unreachable_frac, else reachable=0
        people.reachable[istart:iend] = (draws >= thresholds).astype(np.int8)

    def on_birth(self, istart, iend, tick):
        """Set reachability for newborn agents (called by BirthsByCBR)."""
        self._set_reachability(istart, iend)

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # --- Routine Immunization (weekly) ---
        if tick > 0 and tick % self.ri_period == 0:
            ages = tick - people.dob[:count]
            ri_eligible = np.nonzero(
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
                (people.reachable[:count] == 1) &  # Only reachable agents
                (ages >= self.ri_age) &
                (ages < self.ri_age + self.ri_period)
            )[0]

            if len(ri_eligible) > 0:
                agent_coverage = self.ri_coverage[people.nodeid[ri_eligible]]
                draws = np.random.random(len(ri_eligible)).astype(np.float32)
                vaccinated = ri_eligible[draws < agent_coverage]

                if len(vaccinated) > 0:
                    people.state[vaccinated] = SEIR.State.RECOVERED.value
                    vax_by_node = np.bincount(
                        people.nodeid[vaccinated],
                        minlength=nodes.count
                    ).astype(nodes.S.dtype)
                    nodes.S[tick + 1] -= vax_by_node
                    np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
                    nodes.R[tick + 1] += vax_by_node

        # --- SIA Campaigns (every sia_period days) ---
        if tick > 0 and tick % self.sia_period == 0:
            ages = tick - people.dob[:count]
            sia_eligible = np.nonzero(
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
                (people.reachable[:count] == 1) &  # Only reachable agents
                (ages >= 0) &
                (ages < self.sia_max_age)
            )[0]

            if len(sia_eligible) > 0:
                agent_coverage = self.sia_coverage[people.nodeid[sia_eligible]]
                draws = np.random.random(len(sia_eligible)).astype(np.float32)
                vaccinated = sia_eligible[draws < agent_coverage]

                if len(vaccinated) > 0:
                    people.state[vaccinated] = SEIR.State.RECOVERED.value
                    vax_by_node = np.bincount(
                        people.nodeid[vaccinated],
                        minlength=nodes.count
                    ).astype(nodes.S.dtype)
                    nodes.S[tick + 1] -= vax_by_node
                    np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
                    nodes.R[tick + 1] += vax_by_node


class PulseImportation:
    """Introduces infections into patch 0 at regular intervals.

    Every `period` ticks, selects up to `count` susceptible agents in patch 0
    and sets their state to INFECTIOUS with sampled infectious durations.

    Parameters:
        model: LASER Model instance
        infdurdist: Distribution for infectious duration sampling
        period: Ticks between importation pulses (default: 90)
        count: Number of infections per pulse (default: 5)
    """

    def __init__(self, model, infdurdist, period=90, count=5):
        self.model = model
        self.infdurdist = infdurdist
        self.period = period
        self.count = count

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0:
            return

        people = self.model.people
        nodes = self.model.nodes
        active = people.count

        susceptible_in_patch = np.nonzero(
            (people.state[:active] == SEIR.State.SUSCEPTIBLE.value) &
            (people.nodeid[:active] == 0)
        )[0]

        if len(susceptible_in_patch) == 0:
            return

        n_infect = min(self.count, len(susceptible_in_patch))
        chosen = np.random.choice(susceptible_in_patch, size=n_infect, replace=False)

        people.state[chosen] = SEIR.State.INFECTIOUS.value

        samples = dists.sample_floats(
            self.infdurdist, np.zeros(n_infect, np.float32)
        )
        samples = np.maximum(np.round(samples), 1).astype(people.itimer.dtype)
        people.itimer[chosen] = samples

        nodes.S[tick + 1, 0] -= n_infect
        nodes.S[tick + 1, 0] = max(nodes.S[tick + 1, 0], 0)
        nodes.I[tick + 1, 0] += n_infect


class PatchImportation:
    """Seeds infections in endemic corridor patches to represent cross-border
    and persistent reservoir transmission.

    Only infects susceptible agents (state == SUSCEPTIBLE) in specified patches.

    Parameters:
        model: LASER Model instance
        infdurdist: Distribution for infectious duration sampling
        endemic_patches: List of patch indices to receive importations
        period: Ticks between importation events (default: 30)
        count: Number of infections per patch per event (default: 1)
        end_tick: Stop importation after this tick (default: None = run entire sim)
    """

    def __init__(self, model, infdurdist, endemic_patches,
                 period=30, count=1, end_tick=None):
        self.model = model
        self.infdurdist = infdurdist
        self.endemic_patches = np.asarray(endemic_patches, dtype=np.int32)
        self.period = period
        self.count = count
        self.end_tick = end_tick or model.params.nticks

    def step(self, tick):
        if tick <= 0 or tick % self.period != 0 or tick >= self.end_tick:
            return

        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        for patch_id in self.endemic_patches:
            susceptible_in_patch = np.nonzero(
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
                (people.nodeid[:count] == patch_id)
            )[0]

            if len(susceptible_in_patch) == 0:
                continue

            n_infect = min(self.count, len(susceptible_in_patch))
            chosen = np.random.choice(susceptible_in_patch, size=n_infect, replace=False)

            people.state[chosen] = SEIR.State.INFECTIOUS.value

            samples = dists.sample_floats(
                self.infdurdist, np.zeros(n_infect, np.float32)
            )
            samples = np.maximum(np.round(samples), 1).astype(people.itimer.dtype)
            people.itimer[chosen] = samples

            nodes.S[tick + 1, patch_id] -= n_infect
            nodes.S[tick + 1, patch_id] = max(nodes.S[tick + 1, patch_id], 0)
            nodes.I[tick + 1, patch_id] += n_infect
