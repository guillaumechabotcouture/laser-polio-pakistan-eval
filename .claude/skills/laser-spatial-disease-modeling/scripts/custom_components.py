#!/usr/bin/env python3
"""
Custom LASER components for spatial SEIR modeling with seasonal transmission
and disease importation.

Components:
    Importation - Seeds infections in *susceptible* agents periodically to
                  sustain epidemics in sub-critical-community-size populations.
                  Unlike the built-in Infect_Random_Agents (laser.generic.importation),
                  this only infects susceptible agents, which is more epidemiologically
                  precise. Use the built-in class for simpler setups where state
                  filtering is not required.

    SeasonalTransmission - Extends SEIR.Transmission with time-varying beta
                           using a seasonal forcing profile and spatial coupling.
                           NOTE (v1.0.0): The built-in TransmissionSE (aliased as
                           SEIR.Transmission) now accepts a `seasonality` parameter
                           (ValuesMap or ndarray) and handles spatial network coupling
                           internally. For most use cases, pass seasonality directly:
                             seasonality = ValuesMap.from_timeseries(beta_season_tiled, nnodes)
                             SEIR.Transmission(model, expdurdist, seasonality=seasonality)
                           Use this custom class only when you need non-standard behavior
                           (e.g., tick % 365 cycling, per-node profiles, modified coupling).

Usage:
    from custom_components import Importation, SeasonalTransmission
"""

import numpy as np
import numba as nb
import laser.core.distributions as dists
from laser.generic import SEIR


class Importation:
    """Seeds infections in susceptible agents periodically to sustain epidemics
    in sub-CCS populations.

    Unlike the built-in Infect_Random_Agents, this class filters agents by
    susceptible state before infecting, preventing wasted importation events
    on already-infected or recovered agents.

    Parameters:
        model: LASER Model instance
        infdurdist: Distribution for infectious duration sampling
        infdurmin: Minimum infectious duration in ticks (default: 1)
        period: Ticks between importation events (default: 30)
        count: Number of agents to infect per event (default: 3)
        end_tick: Stop importation after this tick (default: 10*365)
    """

    def __init__(self, model, infdurdist, infdurmin=1, period=30,
                 count=3, end_tick=10*365):
        self.model = model
        self.infdurdist = infdurdist
        self.infdurmin = infdurmin
        self.period = period
        self.count = count
        self.end_tick = end_tick or model.params.nticks
        self.model.nodes.add_vector_property(
            "imports", model.params.nticks + 1, dtype=np.uint32, default=0
        )

    def step(self, tick):
        if tick > 0 and tick % self.period == 0 and tick < self.end_tick:
            i_susceptible = np.nonzero(
                self.model.people.state == SEIR.State.SUSCEPTIBLE.value
            )[0]
            if len(i_susceptible) > 0:
                count = min(self.count, len(i_susceptible))
                i_infect = np.random.choice(i_susceptible, size=count, replace=False)
                self.model.people.state[i_infect] = SEIR.State.INFECTIOUS.value
                samples = dists.sample_floats(
                    self.infdurdist, np.zeros(count, np.float32)
                )
                samples = np.maximum(
                    np.round(samples), self.infdurmin
                ).astype(self.model.people.itimer.dtype)
                self.model.people.itimer[i_infect] = samples
                inf_by_node = np.bincount(
                    self.model.people.nodeid[i_infect],
                    minlength=len(self.model.nodes)
                ).astype(self.model.nodes.S.dtype)
                self.model.nodes.S[tick + 1] -= inf_by_node
                self.model.nodes.I[tick + 1] += inf_by_node
                self.model.nodes.imports[tick] = inf_by_node


class SeasonalTransmission(SEIR.Transmission):
    """Extends SEIR.Transmission with time-varying beta via seasonal forcing.

    ADVANCED CUSTOMIZATION EXAMPLE. For most use cases, prefer the built-in
    seasonality parameter on SEIR.Transmission (TransmissionSE):
        seasonality = ValuesMap.from_timeseries(beta_season_tiled, nnodes)
        SEIR.Transmission(model, expdurdist, seasonality=seasonality)

    This custom class differs from the built-in in that it:
    - Uses tick % 365 cycling (auto-repeats the 365-day profile)
    - Reads beta_season from model.params rather than a ValuesMap

    Expects model.params to contain:
        beta: Base transmission rate
        beta_season: 365-element array of seasonal modulation factors

    Expects model.network to be set (gravity coupling matrix).
    """

    def step(self, tick):
        ft = self.model.nodes.forces[tick]
        N = (self.model.nodes.S[tick] + self.model.nodes.E[tick] +
             (I := self.model.nodes.I[tick]))
        if hasattr(self.model.nodes, "R"):
            N += self.model.nodes.R[tick]

        # Seasonal beta modulation
        ft[:] = (self.model.params.beta * I / N *
                 self.model.params.beta_season[tick % 365])

        # Spatial coupling via network
        transfer = ft[:, None] * self.model.network
        ft += transfer.sum(axis=0)
        ft -= transfer.sum(axis=1)
        ft = -np.expm1(-ft)  # Convert rate to probability

        newly_infected_by_node = np.zeros(
            (nb.get_num_threads(), self.model.nodes.count), dtype=np.int32
        )
        self.nb_transmission_step(
            self.model.people.state, self.model.people.nodeid, ft,
            newly_infected_by_node, self.model.people.etimer,
            self.expdurdist, self.expdurmin, tick,
        )
        newly_infected_by_node = newly_infected_by_node.sum(axis=0).astype(
            self.model.nodes.S.dtype
        )
        self.model.nodes.S[tick + 1] -= newly_infected_by_node
        self.model.nodes.E[tick + 1] += newly_infected_by_node
        self.model.nodes.newly_infected[tick] = newly_infected_by_node
