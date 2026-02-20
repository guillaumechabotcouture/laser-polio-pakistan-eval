#!/usr/bin/env python3
"""
Custom LASER components for Pakistan polio SEIR model.

Components:
    PerPatchVaccination - Routine immunization + SIA campaigns with per-district
                          coverage rates and CORRELATED MISSEDNESS via per-agent
                          reachable flag. Sets state=RECOVERED (not susceptibility=0)
                          because the built-in TransmissionSE kernel only checks
                          state==SUSCEPTIBLE for infection eligibility.

    PatchImportation    - Seeds infections in endemic corridor patches only
                          to represent cross-border and persistent transmission foci.

Usage:
    from custom_components import PerPatchVaccination, PatchImportation
"""

import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR


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
