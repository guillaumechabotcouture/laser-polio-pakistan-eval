#!/usr/bin/env python3
"""
Custom LASER components for Chad guinea worm dual-host SEIS model.

Guinea worm (Dracunculus medinensis) differs from standard SEIR diseases:
  - NO lasting immunity: recovered hosts return to susceptible (SEIS)
  - Water-mediated transmission: larvae in stagnant water via copepods
  - Dual-host: dogs (~95% of cases in Chad) and humans share water sources
  - Long pre-patent period: ~365 days (E), short emergence: ~21 days (I)
  - Dry-season peak: stagnant pools concentrate copepods (Jan-Mar)

Core components:
    SEISRecovery                  - I→S transition (no lasting immunity)
    build_chad_dry_season         - Chad dry-season transmission profile
    DualHostWaterTransmission     - Cross-species FOI via shared water sources
    PatchImportation              - Susceptible-targeted infection seeding

Prevention interventions (no vaccine exists for guinea worm):
    ABATELarvicide                - Monthly temephos application to water sources
    WaterFilterDistribution       - Cloth/pipe filters reducing copepod ingestion
    CaseContainment               - Human case detection and water-access prevention
    DogTethering                  - Infected dog tethering to prevent water access

Intervention-aware transmission:
    IntervenedWaterTransmission       - Replaces SEIR.Transmission with FOI modifiers
    IntervenedDualHostTransmission    - Replaces DualHostWaterTransmission with FOI modifiers
"""

import numpy as np
import laser.core.distributions as dists
from laser.generic import SEIR


class SEISRecovery:
    """Replaces SEIR.Recovered to implement SEIS dynamics.

    Standard SEIR routes I→R (permanent immunity). Guinea worm confers
    NO lasting immunity, so recovered hosts return to susceptible (I→S).

    This component:
    1. Detects agents whose itimer has expired (state was set to RECOVERED
       by SEIR.Infectious)
    2. Immediately resets them to SUSCEPTIBLE
    3. Updates node-level S and R counts accordingly

    Must be placed AFTER SEIR.Infectious and SEIR.Recovered in the
    component list so that I→R transitions have already occurred.
    """

    def __init__(self, model):
        self.model = model

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # Find agents that just transitioned to RECOVERED this tick.
        # SEIR.Infectious sets state=RECOVERED when itimer expires.
        # SEIR.Recovered then counts them in nodes.R[tick+1].
        # We now move them back to SUSCEPTIBLE.
        recovered_mask = (people.state[:count] == SEIR.State.RECOVERED.value)
        recovered_idx = np.nonzero(recovered_mask)[0]

        if len(recovered_idx) == 0:
            return

        # Set state back to SUSCEPTIBLE
        people.state[recovered_idx] = SEIR.State.SUSCEPTIBLE.value

        # Update node counts: move from R to S
        r_by_node = np.bincount(
            people.nodeid[recovered_idx], minlength=nodes.count
        ).astype(nodes.R.dtype)

        nodes.R[tick + 1] -= r_by_node
        np.maximum(nodes.R[tick + 1], 0, out=nodes.R[tick + 1])
        nodes.S[tick + 1] += r_by_node


def build_chad_dry_season(nticks, nnodes):
    """Build 365-day seasonal forcing for Chad's dry-season guinea worm peak.

    Transmission biology:
    - Dry season (Oct-May): Stagnant pools concentrate copepods carrying
      guinea worm larvae. Fewer water sources force sharing between hosts.
    - Peak (Jan-Mar): Maximum stagnation, 1.8x baseline transmission
    - Decline (Apr-May): Pools begin refilling, copepod density drops
    - Rainy season (Jun-Sep): Flowing water disperses larvae, dilutes
      copepods, new water sources reduce sharing. 0.4x baseline.

    Returns:
        seasonality_values: ndarray shape (nticks, nnodes) for ValuesMap
        season_365: 365-day profile (for plotting)
    """
    from laser.generic.utils import ValuesMap

    days = np.arange(365)
    season_365 = np.ones(365, dtype=np.float64)

    for d in range(365):
        month = 1 + (d * 12) // 365  # approximate month (1-12)

        if month in (1, 2, 3):       # Jan-Mar: peak dry season
            season_365[d] = 1.8
        elif month in (4, 5):         # Apr-May: declining dry
            # Linear decline from 1.8 to 1.0
            frac = (d - 90) / 61.0    # days 90-151
            frac = np.clip(frac, 0, 1)
            season_365[d] = 1.8 - 0.8 * frac
        elif month in (6, 7, 8, 9):   # Jun-Sep: rainy season
            season_365[d] = 0.4
        elif month == 10:             # Oct: transition to dry
            # Ramp from 0.4 to 1.0
            frac = (d - 273) / 30.0
            frac = np.clip(frac, 0, 1)
            season_365[d] = 0.4 + 0.6 * frac
        else:                         # Nov-Dec: early dry season
            # Ramp from 1.0 to 1.8
            frac = (d - 304) / 61.0
            frac = np.clip(frac, 0, 1)
            season_365[d] = 1.0 + 0.8 * frac

    # Normalize to mean 1.0 so beta represents the average transmission rate
    season_365 /= season_365.mean()
    assert abs(season_365.mean() - 1.0) < 0.01, \
        f"Seasonal profile mean={season_365.mean():.3f}, must be ~1.0"

    season_tiled = np.tile(season_365, nticks // 365 + 1)[:nticks]
    seasonality = ValuesMap.from_timeseries(season_tiled, nnodes)

    return seasonality, season_365


class DualHostWaterTransmission:
    """Cross-species transmission via shared water sources.

    Guinea worm transmission is water-mediated:
    1. Infected host enters water → larvae released
    2. Copepods ingest larvae in stagnant water
    3. New host drinks unfiltered water containing infected copepods
    4. ~365-day pre-patent period before worm emergence

    Both dogs and humans contaminate and drink from the SAME water bodies.
    In Chad, dogs are the dominant reservoir (~95% of detected infections),
    so dog-to-water-to-human transmission is the primary pathway.

    Cross-species coupling is asymmetric:
    - dog_to_human: fraction of dog FOI that applies to humans (high: 0.8)
    - human_to_dog: fraction of human FOI that applies to dogs (low: 0.2)
    This reflects dogs' much higher contribution to water contamination
    and their greater propensity to enter water sources.

    Each species also has within-species transmission via shared water,
    handled by the standard gravity-coupled TransmissionSE on each model.

    This component runs AFTER TransmissionSE on each model and adds the
    cross-species contribution. It applies new infections to susceptible
    agents via Bernoulli trials.

    Parameters:
        human_model: LASER Model for human population
        dog_model: LASER Model for dog population
        expdurdist: Distribution for exposed duration (shared across species)
        dog_to_human: Cross-species coupling coefficient (default: 0.8)
        human_to_dog: Cross-species coupling coefficient (default: 0.2)
        seasonality_values: ndarray shape (nticks,) of seasonal multipliers
    """

    def __init__(self, human_model, dog_model, expdurdist,
                 dog_to_human=0.8, human_to_dog=0.2,
                 seasonality_values=None):
        self.human_model = human_model
        self.dog_model = dog_model
        self.expdurdist = expdurdist
        self.dog_to_human = dog_to_human
        self.human_to_dog = human_to_dog
        self.seasonality_values = seasonality_values

        # Both models must have the same number of patches
        assert human_model.nodes.count == dog_model.nodes.count, \
            "Human and dog models must have the same number of patches"

        self.nnodes = human_model.nodes.count

    def _get_seasonal_multiplier(self, tick):
        """Get seasonal multiplier for the current tick."""
        if self.seasonality_values is None:
            return 1.0
        return self.seasonality_values[tick % len(self.seasonality_values)]

    def _compute_cross_foi(self, source_model, tick):
        """Compute per-patch force of infection from the source species.

        Returns FOI as I_source / N_source per patch (prevalence-based).
        """
        nodes = source_model.nodes
        I = nodes.I[tick].astype(np.float64)
        N = (nodes.S[tick] + nodes.E[tick] + nodes.I[tick] + nodes.R[tick]).astype(np.float64)
        N = np.maximum(N, 1.0)  # prevent division by zero
        return I / N

    def _apply_cross_infections(self, target_model, cross_foi, coupling, tick):
        """Apply cross-species infections to susceptible agents in target.

        cross_foi: per-patch prevalence from the source species
        coupling: cross-species coupling coefficient
        """
        people = target_model.people
        nodes = target_model.nodes
        count = people.count

        season = self._get_seasonal_multiplier(tick)

        # Per-patch infection probability from cross-species source
        beta = target_model.params.beta
        p_cross = 1.0 - np.exp(-beta * coupling * season * cross_foi)

        # Find susceptible agents
        susc_mask = (people.state[:count] == SEIR.State.SUSCEPTIBLE.value)
        susc_idx = np.nonzero(susc_mask)[0]

        if len(susc_idx) == 0:
            return

        # Per-agent infection probability based on their patch
        agent_probs = p_cross[people.nodeid[susc_idx]]

        # Bernoulli trials
        draws = np.random.random(len(susc_idx))
        infected_mask = draws < agent_probs
        newly_infected = susc_idx[infected_mask]

        if len(newly_infected) == 0:
            return

        # Transition S → E
        people.state[newly_infected] = SEIR.State.EXPOSED.value

        # Sample exposed durations and set etimer
        etimer_samples = dists.sample_floats(
            self.expdurdist, np.zeros(len(newly_infected), dtype=np.float32)
        )
        etimer_samples = np.maximum(np.round(etimer_samples), 1).astype(people.etimer.dtype)
        people.etimer[newly_infected] = etimer_samples

        # Update node counts
        inf_by_node = np.bincount(
            people.nodeid[newly_infected], minlength=nodes.count
        ).astype(nodes.S.dtype)

        nodes.S[tick + 1] -= inf_by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.E[tick + 1] += inf_by_node

    def step(self, tick):
        """Apply cross-species transmission for the current tick.

        Dog prevalence → infections in humans (strong coupling)
        Human prevalence → infections in dogs (weak coupling)
        """
        # Dog-to-human: dog prevalence drives human infections
        dog_foi = self._compute_cross_foi(self.dog_model, tick)
        self._apply_cross_infections(
            self.human_model, dog_foi, self.dog_to_human, tick
        )

        # Human-to-dog: human prevalence drives dog infections
        human_foi = self._compute_cross_foi(self.human_model, tick)
        self._apply_cross_infections(
            self.dog_model, human_foi, self.human_to_dog, tick
        )


class PatchImportation:
    """Seeds infections in endemic patches targeting susceptible agents only.

    Parameters:
        model: LASER Model instance
        infdurdist: Distribution for infectious duration (not used for SEIS
                    where I→S, but needed for timer initialization)
        expdurdist: Distribution for exposed duration sampling
        endemic_patches: List of patch indices to receive importations
        period: Ticks between importation events (default: 90)
        count: Number of infections per event total (distributed across patches)
        end_tick: Stop importation after this tick (default: None = entire sim)
    """

    def __init__(self, model, expdurdist, endemic_patches,
                 period=90, count=1, end_tick=None):
        self.model = model
        self.expdurdist = expdurdist
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

        # Distribute importations across endemic patches
        for patch_id in self.endemic_patches:
            susceptible_in_patch = np.nonzero(
                (people.state[:count] == SEIR.State.SUSCEPTIBLE.value) &
                (people.nodeid[:count] == patch_id)
            )[0]

            if len(susceptible_in_patch) == 0:
                continue

            n_infect = min(self.count, len(susceptible_in_patch))
            chosen = np.random.choice(susceptible_in_patch, size=n_infect, replace=False)

            # Infect as EXPOSED (will progress through E→I→S cycle)
            people.state[chosen] = SEIR.State.EXPOSED.value

            samples = dists.sample_floats(
                self.expdurdist, np.zeros(n_infect, np.float32)
            )
            samples = np.maximum(np.round(samples), 1).astype(people.etimer.dtype)
            people.etimer[chosen] = samples

            nodes.S[tick + 1, patch_id] -= n_infect
            nodes.S[tick + 1, patch_id] = max(nodes.S[tick + 1, patch_id], 0)
            nodes.E[tick + 1, patch_id] += n_infect


# ====================================================================
# Prevention-Based Interventions
# ====================================================================
#
# Guinea worm has no vaccine. All control relies on breaking the
# water transmission cycle through four interventions:
#
#   1. ABATE larvicide — kills copepods in water sources (patch-level)
#   2. Water filters — removes copepods from drinking water (agent-level)
#   3. Case containment — prevents infectious humans from contaminating water
#   4. Dog tethering — prevents infectious dogs from accessing water
#
# These modify the force of infection through IntervenedWaterTransmission
# (within-species) and IntervenedDualHostTransmission (cross-species).
#
# Full FOI formula per susceptible agent j in patch i:
#
#   I_eff_i = sum_k(contribution_k)  where contribution is:
#       1.0        if agent k is not contained/tethered
#       1-efficacy if agent k is contained (human) or tethered (dog)
#
#   prev_i = I_eff_i / N_i
#   coupled_prev_i = (1 - export_i) * prev_i + sum_k(network_ik * prev_k)
#   FOI_i = beta * season(t) * coupled_prev_i * abate_factor_i
#   P(infect j) = (1 - exp(-FOI_i)) * (1 - has_filter_j * filter_efficacy)
#
# ====================================================================


class ABATELarvicide:
    """ABATE (temephos) larvicide applied monthly to known water sources.

    During transmission season, ABATE is applied to stagnant water bodies,
    killing copepods that carry guinea worm larvae. This reduces the
    per-patch force of infection by coverage * efficacy.

    Sets nodes.abate_factor[tick+1] each tick: 1.0 = no effect,
    (1 - coverage * efficacy) during active season.

    Must be added to BOTH human and dog models (same water sources).

    Parameters:
        model: LASER Model instance
        coverage_by_patch: 1D array of per-patch ABATE coverage (0-1),
            or scalar for uniform coverage. Typical range 0.20-0.80.
        efficacy: Fraction reduction in copepod-borne transmission (default: 0.80)
        season_start_day: Day of year when application begins (default: 0 = Jan 1)
        season_end_day: Day of year when application ends (default: 150 = ~May 30)
    """

    def __init__(self, model, coverage_by_patch, efficacy=0.80,
                 season_start_day=0, season_end_day=150):
        self.model = model
        nnodes = model.nodes.count

        if np.isscalar(coverage_by_patch):
            self.coverage = np.full(nnodes, coverage_by_patch, dtype=np.float32)
        else:
            self.coverage = np.asarray(coverage_by_patch, dtype=np.float32)

        assert len(self.coverage) == nnodes, \
            f"ABATE coverage length {len(self.coverage)} != {nnodes} patches"
        assert np.all((self.coverage >= 0) & (self.coverage <= 1)), \
            f"ABATE coverage must be in [0, 1], got [{self.coverage.min():.2f}, {self.coverage.max():.2f}]"

        self.efficacy = efficacy
        self.season_start = season_start_day
        self.season_end = season_end_day

        # Per-patch transmission multiplier: 1.0 = no effect, <1.0 = reduced
        if not hasattr(model.nodes, "abate_factor"):
            model.nodes.add_vector_property(
                "abate_factor", model.params.nticks + 1, dtype=np.float32
            )
            model.nodes.abate_factor[:] = 1.0

    def step(self, tick):
        day_of_year = tick % 365

        if self.season_start <= day_of_year <= self.season_end:
            # Active season: factor = 1 - (coverage * efficacy)
            self.model.nodes.abate_factor[tick + 1] = \
                1.0 - self.coverage * self.efficacy
        else:
            # Off-season: no ABATE effect
            self.model.nodes.abate_factor[tick + 1] = 1.0


class WaterFilterDistribution:
    """Cloth/pipe water filter distribution to households.

    Filters remove copepods from drinking water, reducing individual
    ingestion risk. Each agent is assigned a filter status at creation
    based on the adoption rate, representing household-level distribution.

    Filter protection is applied during transmission (by
    IntervenedWaterTransmission): agents with filters have their
    per-tick infection probability multiplied by (1 - filter_efficacy).

    Only added to the human model (dogs don't use filters).

    Parameters:
        model: LASER Model instance (human)
        adoption_rate: Fraction of population with filters (scalar or
            per-patch array). Default: 0.60 (60% nationally).
        filter_efficacy: Fraction risk reduction per filter (default: 0.95)
    """

    def __init__(self, model, adoption_rate=0.60, filter_efficacy=0.95):
        self.model = model
        nnodes = model.nodes.count

        if np.isscalar(adoption_rate):
            self.adoption_rate = np.full(nnodes, adoption_rate, dtype=np.float32)
        else:
            self.adoption_rate = np.asarray(adoption_rate, dtype=np.float32)

        self.filter_efficacy = filter_efficacy

        if not hasattr(model.people, "has_filter"):
            model.people.add_scalar_property("has_filter", dtype=np.int8, default=0)

        # Assign filters to existing agents
        self._assign_filters(0, model.people.count)

    def _assign_filters(self, istart, iend):
        """Assign filter status based on patch-level adoption rate."""
        n = iend - istart
        if n <= 0:
            return
        nodeids = self.model.people.nodeid[istart:iend]
        rates = self.adoption_rate[nodeids]
        draws = np.random.random(n).astype(np.float32)
        self.model.people.has_filter[istart:iend] = (draws < rates).astype(np.int8)

    def on_birth(self, istart, iend, tick):
        """Newborns inherit household filter status probabilistically."""
        self._assign_filters(istart, iend)

    def step(self, tick):
        # Filters are passive — protection applied during transmission.
        pass


class CaseContainment:
    """Case containment for human guinea worm infections.

    When a worm emerges (agent becomes INFECTIOUS), surveillance detects
    the case with probability detection_rate. Detected cases are prevented
    from entering water sources for the ~30-day emergence period, reducing
    their contribution to water contamination by containment_efficacy.

    Containment status uses three states:
        0 = pending (newly infectious, no decision yet)
        1 = contained (detected, water access prevented)
        2 = uncontained (escaped surveillance)

    Containment resets when the agent leaves INFECTIOUS state
    (I→R→S via SEISRecovery), so reinfection triggers fresh detection.

    Only added to the human model.

    Parameters:
        model: LASER Model instance (human)
        detection_rate: Fraction of cases detected (default: 0.70)
        containment_efficacy: Reduction in water contamination from
            contained cases (default: 0.90)
    """

    def __init__(self, model, detection_rate=0.70, containment_efficacy=0.90):
        self.model = model
        self.detection_rate = detection_rate
        self.containment_efficacy = containment_efficacy

        if not hasattr(model.people, "contained"):
            model.people.add_scalar_property("contained", dtype=np.int8, default=0)

        # Track containment counts per patch per tick
        model.nodes.add_vector_property(
            "cases_contained", model.params.nticks + 1, dtype=np.int32
        )

    def on_birth(self, istart, iend, tick):
        """Newborns start with no containment status."""
        self.model.people.contained[istart:iend] = 0

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # Reset containment for agents who are no longer infectious
        # (handles SEIS cycle: I→R→S resets the flag for future reinfection)
        not_infectious = (people.state[:count] != SEIR.State.INFECTIOUS.value)
        people.contained[:count][not_infectious] = 0

        # Find infectious agents needing a containment decision
        pending = np.nonzero(
            (people.state[:count] == SEIR.State.INFECTIOUS.value) &
            (people.contained[:count] == 0)
        )[0]

        if len(pending) == 0:
            return

        # Attempt detection and containment
        draws = np.random.random(len(pending))
        detected = pending[draws < self.detection_rate]
        escaped = pending[draws >= self.detection_rate]

        people.contained[detected] = 1   # successfully contained
        people.contained[escaped] = 2    # escaped surveillance

        if len(detected) > 0:
            contained_by_node = np.bincount(
                people.nodeid[detected], minlength=nodes.count
            ).astype(np.int32)
            nodes.cases_contained[tick] += contained_by_node


class DogTethering:
    """Dog tethering/management for infected dogs.

    In villages with known dog infections, infected dogs are tethered
    or managed to prevent water access during worm emergence, reducing
    their contribution to water contamination.

    Mechanics mirror CaseContainment but for the dog population:
        0 = pending, 1 = tethered, 2 = not tethered

    Status resets when the dog leaves INFECTIOUS state (SEIS cycle).

    Only added to the dog model.

    Parameters:
        model: LASER Model instance (dog)
        village_coverage: Fraction of infected dogs tethered (default: 0.50)
        tether_efficacy: Reduction in water contamination from
            tethered dogs (default: 0.85)
    """

    def __init__(self, model, village_coverage=0.50, tether_efficacy=0.85):
        self.model = model
        self.village_coverage = village_coverage
        self.tether_efficacy = tether_efficacy

        if not hasattr(model.people, "tethered"):
            model.people.add_scalar_property("tethered", dtype=np.int8, default=0)

        model.nodes.add_vector_property(
            "dogs_tethered", model.params.nticks + 1, dtype=np.int32
        )

    def on_birth(self, istart, iend, tick):
        """Newborn dogs start untethered."""
        self.model.people.tethered[istart:iend] = 0

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count

        # Reset for dogs no longer infectious (SEIS cycle)
        not_infectious = (people.state[:count] != SEIR.State.INFECTIOUS.value)
        people.tethered[:count][not_infectious] = 0

        # Find infectious dogs needing a tethering decision
        pending = np.nonzero(
            (people.state[:count] == SEIR.State.INFECTIOUS.value) &
            (people.tethered[:count] == 0)
        )[0]

        if len(pending) == 0:
            return

        draws = np.random.random(len(pending))
        caught = pending[draws < self.village_coverage]
        free = pending[draws >= self.village_coverage]

        people.tethered[caught] = 1   # successfully tethered
        people.tethered[free] = 2     # not tethered

        if len(caught) > 0:
            tethered_by_node = np.bincount(
                people.nodeid[caught], minlength=nodes.count
            ).astype(np.int32)
            nodes.dogs_tethered[tick] += tethered_by_node


# ====================================================================
# Intervention-Aware Transmission
# ====================================================================


class IntervenedWaterTransmission:
    """Within-species water-borne transmission with intervention effects.

    Replaces SEIR.Transmission (TransmissionSE) for guinea worm models
    with active prevention interventions. Computes gravity-coupled FOI
    with four intervention modifiers:

    1. Containment/tethering: reduces effective infectious count
       I_eff = sum(1.0 if not contained/tethered,
                   1-efficacy if contained/tethered)

    2. ABATE larvicide: patch-level FOI reduction
       FOI *= nodes.abate_factor[tick+1]  (set by ABATELarvicide)

    3. Seasonality: standard seasonal forcing
       FOI *= season(tick)

    4. Water filters: per-agent protection
       P(infection) *= (1 - filter_efficacy) if agent has_filter

    Auto-detects which intervention properties exist on the model
    (contained, tethered, has_filter, abate_factor) via hasattr checks,
    so the same class works for both human and dog models.

    Parameters:
        model: LASER Model instance
        expdurdist: Distribution for exposed duration sampling
        seasonality_values: 1D ndarray of shape (nticks,) or (>=nticks,)
            with seasonal transmission multipliers. None = no seasonality.
        containment_efficacy: FOI reduction for contained humans (default: 0.90)
        tether_efficacy: FOI reduction for tethered dogs (default: 0.85)
        filter_efficacy: Per-agent risk reduction from filters (default: 0.95)
    """

    def __init__(self, model, expdurdist, seasonality_values=None,
                 containment_efficacy=0.90, tether_efficacy=0.85,
                 filter_efficacy=0.95):
        self.model = model
        self.expdurdist = expdurdist
        self.seasonality_values = seasonality_values
        self.containment_efficacy = containment_efficacy
        self.tether_efficacy = tether_efficacy
        self.filter_efficacy = filter_efficacy

        # Track new infections (normally added by TransmissionSE)
        if not hasattr(model.nodes, "newly_infected"):
            model.nodes.add_vector_property(
                "newly_infected", model.params.nticks + 1, dtype=np.int32
            )

    def step(self, tick):
        people = self.model.people
        nodes = self.model.nodes
        count = people.count
        nnodes = nodes.count

        # --- Effective prevalence per patch ---
        N = (nodes.S[tick].astype(np.float64) +
             nodes.E[tick].astype(np.float64) +
             nodes.I[tick].astype(np.float64) +
             nodes.R[tick].astype(np.float64))
        N = np.maximum(N, 1.0)

        infectious_idx = np.nonzero(
            people.state[:count] == SEIR.State.INFECTIOUS.value
        )[0]

        if len(infectious_idx) == 0:
            return

        # Per-agent infectiousness contribution (reduced if contained/tethered)
        contributions = np.ones(len(infectious_idx), dtype=np.float64)

        if hasattr(people, "contained"):
            contained = (people.contained[infectious_idx] == 1)
            contributions[contained] *= (1.0 - self.containment_efficacy)

        if hasattr(people, "tethered"):
            tethered = (people.tethered[infectious_idx] == 1)
            contributions[tethered] *= (1.0 - self.tether_efficacy)

        # Sum effective I by patch
        I_eff = np.zeros(nnodes, dtype=np.float64)
        np.add.at(I_eff, people.nodeid[infectious_idx], contributions)

        prev = I_eff / N

        # --- Gravity-coupled prevalence ---
        network = self.model.network
        coupled_prev = np.zeros(nnodes, dtype=np.float64)
        for i in range(nnodes):
            local_frac = 1.0 - network[i].sum()
            coupled_prev[i] = local_frac * prev[i] + np.dot(network[i], prev)

        # --- Seasonal multiplier ---
        if self.seasonality_values is not None:
            season = self.seasonality_values[tick % len(self.seasonality_values)]
        else:
            season = 1.0

        # --- ABATE patch-level reduction ---
        if hasattr(nodes, "abate_factor"):
            abate = nodes.abate_factor[tick + 1].astype(np.float64)
        else:
            abate = np.ones(nnodes, dtype=np.float64)

        # --- Per-patch FOI ---
        beta = self.model.params.beta
        foi = beta * season * coupled_prev * abate

        # --- Per-agent infection probability ---
        susc_idx = np.nonzero(
            people.state[:count] == SEIR.State.SUSCEPTIBLE.value
        )[0]

        if len(susc_idx) == 0:
            return

        agent_foi = foi[people.nodeid[susc_idx]]
        p_inf = 1.0 - np.exp(-agent_foi)

        # Filter protection (human model only)
        if hasattr(people, "has_filter"):
            filtered = people.has_filter[susc_idx].astype(bool)
            p_inf[filtered] *= (1.0 - self.filter_efficacy)

        # --- Bernoulli infection trials ---
        draws = np.random.random(len(susc_idx))
        newly_infected = susc_idx[draws < p_inf]

        if len(newly_infected) == 0:
            return

        # Transition S → E
        people.state[newly_infected] = SEIR.State.EXPOSED.value

        etimer_samples = dists.sample_floats(
            self.expdurdist, np.zeros(len(newly_infected), dtype=np.float32)
        )
        etimer_samples = np.maximum(
            np.round(etimer_samples), 1
        ).astype(people.etimer.dtype)
        people.etimer[newly_infected] = etimer_samples

        # Update node counts
        inf_by_node = np.bincount(
            people.nodeid[newly_infected], minlength=nnodes
        ).astype(nodes.S.dtype)

        nodes.S[tick + 1] -= inf_by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.E[tick + 1] += inf_by_node
        nodes.newly_infected[tick] += inf_by_node


class IntervenedDualHostTransmission:
    """Cross-species transmission via shared water with intervention effects.

    Extends DualHostWaterTransmission to account for prevention interventions
    in the cross-species force of infection:

    1. Source containment/tethering: reduces effective I in the source
       species' prevalence calculation
    2. ABATE in target patch: reduces copepod density at the water source
    3. Target agent filters: reduces individual ingestion probability

    Parameters:
        human_model: LASER Model for human population
        dog_model: LASER Model for dog population
        expdurdist: Distribution for exposed duration sampling
        dog_to_human: Cross-species coupling coefficient (default: 0.15)
        human_to_dog: Cross-species coupling coefficient (default: 0.05)
        seasonality_values: 1D ndarray of seasonal multipliers, or None
        containment_efficacy: Human case containment efficacy (default: 0.90)
        tether_efficacy: Dog tethering efficacy (default: 0.85)
        filter_efficacy: Water filter efficacy (default: 0.95)
    """

    def __init__(self, human_model, dog_model, expdurdist,
                 dog_to_human=0.15, human_to_dog=0.05,
                 seasonality_values=None,
                 containment_efficacy=0.90, tether_efficacy=0.85,
                 filter_efficacy=0.95):
        self.human_model = human_model
        self.dog_model = dog_model
        self.expdurdist = expdurdist
        self.dog_to_human = dog_to_human
        self.human_to_dog = human_to_dog
        self.seasonality_values = seasonality_values
        self.containment_efficacy = containment_efficacy
        self.tether_efficacy = tether_efficacy
        self.filter_efficacy = filter_efficacy

        assert human_model.nodes.count == dog_model.nodes.count, \
            "Human and dog models must have the same number of patches"
        self.nnodes = human_model.nodes.count

    def _get_seasonal_multiplier(self, tick):
        if self.seasonality_values is None:
            return 1.0
        return self.seasonality_values[tick % len(self.seasonality_values)]

    def _compute_effective_prevalence(self, model, tick):
        """Compute I_eff / N per patch, accounting for containment/tethering."""
        people = model.people
        nodes = model.nodes
        count = people.count
        nnodes = nodes.count

        N = (nodes.S[tick].astype(np.float64) + nodes.E[tick].astype(np.float64) +
             nodes.I[tick].astype(np.float64) + nodes.R[tick].astype(np.float64))
        N = np.maximum(N, 1.0)

        infectious_idx = np.nonzero(
            people.state[:count] == SEIR.State.INFECTIOUS.value
        )[0]

        if len(infectious_idx) == 0:
            return np.zeros(nnodes, dtype=np.float64)

        contributions = np.ones(len(infectious_idx), dtype=np.float64)

        if hasattr(people, "contained"):
            contained = (people.contained[infectious_idx] == 1)
            contributions[contained] *= (1.0 - self.containment_efficacy)

        if hasattr(people, "tethered"):
            tethered = (people.tethered[infectious_idx] == 1)
            contributions[tethered] *= (1.0 - self.tether_efficacy)

        I_eff = np.zeros(nnodes, dtype=np.float64)
        np.add.at(I_eff, people.nodeid[infectious_idx], contributions)

        return I_eff / N

    def _apply_cross_infections(self, target_model, cross_prev, coupling, tick):
        """Apply cross-species infections with intervention modifiers."""
        people = target_model.people
        nodes = target_model.nodes
        count = people.count
        nnodes = nodes.count

        season = self._get_seasonal_multiplier(tick)
        beta = target_model.params.beta

        # ABATE reduction on target patch water sources
        if hasattr(nodes, "abate_factor"):
            abate = nodes.abate_factor[tick + 1].astype(np.float64)
        else:
            abate = np.ones(nnodes, dtype=np.float64)

        # Per-patch cross-species infection probability
        p_cross = 1.0 - np.exp(-beta * coupling * season * cross_prev * abate)

        susc_idx = np.nonzero(
            people.state[:count] == SEIR.State.SUSCEPTIBLE.value
        )[0]

        if len(susc_idx) == 0:
            return

        agent_probs = p_cross[people.nodeid[susc_idx]]

        # Filter protection on target agents (human model only)
        if hasattr(people, "has_filter"):
            filtered = people.has_filter[susc_idx].astype(bool)
            agent_probs[filtered] *= (1.0 - self.filter_efficacy)

        # Bernoulli infection trials
        draws = np.random.random(len(susc_idx))
        newly_infected = susc_idx[draws < agent_probs]

        if len(newly_infected) == 0:
            return

        # Transition S → E
        people.state[newly_infected] = SEIR.State.EXPOSED.value

        etimer_samples = dists.sample_floats(
            self.expdurdist, np.zeros(len(newly_infected), dtype=np.float32)
        )
        etimer_samples = np.maximum(
            np.round(etimer_samples), 1
        ).astype(people.etimer.dtype)
        people.etimer[newly_infected] = etimer_samples

        # Update node counts
        inf_by_node = np.bincount(
            people.nodeid[newly_infected], minlength=nnodes
        ).astype(nodes.S.dtype)

        nodes.S[tick + 1] -= inf_by_node
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])
        nodes.E[tick + 1] += inf_by_node

        # Track cross-species infections in newly_infected
        if hasattr(nodes, "newly_infected"):
            nodes.newly_infected[tick] += inf_by_node

    def step(self, tick):
        """Apply cross-species transmission for the current tick.

        Dog prevalence → infections in humans (primary pathway)
        Human prevalence → infections in dogs (weak pathway)
        """
        # Dog-to-human: dog prevalence drives human infections
        dog_prev = self._compute_effective_prevalence(self.dog_model, tick)
        self._apply_cross_infections(
            self.human_model, dog_prev, self.dog_to_human, tick
        )

        # Human-to-dog: human prevalence drives dog infections
        human_prev = self._compute_effective_prevalence(self.human_model, tick)
        self._apply_cross_infections(
            self.dog_model, human_prev, self.human_to_dog, tick
        )
