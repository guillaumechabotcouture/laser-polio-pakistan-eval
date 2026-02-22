Now I have everything I need. Here is a clear, self-contained example with both the wrong approach and the correct fix:

```python
"""
Vaccination in LASER: correct state-based approach.

THE BUG: ImmunizationCampaign / RoutineImmunization set susceptibility = 0.
         TransmissionSE (and TransmissionSI, TransmissionSIx) never check
         susceptibility. They check:
             state == SUSCEPTIBLE  (i.e., state == 0)
         So setting susceptibility = 0 has zero effect on transmission.

THE FIX: Set state = RECOVERED (value 3). TransmissionSE will skip any agent
         whose state is not SUSCEPTIBLE, so vaccinated agents cannot be infected.
         You must also update the node-level S and R counts so the force-of-
         infection denominator (N = S + E + I + R) stays correct.
"""

import numpy as np
from laser.generic import SEIR


# ── WRONG ─────────────────────────────────────────────────────────────────────
# This is what ImmunizationCampaign does. It looks plausible but has no effect.

class WrongVaccination:
    """Sets susceptibility = 0.  TransmissionSE ignores susceptibility → BROKEN."""

    def __init__(self, model, coverage=0.80):
        self.model = model
        self.coverage = coverage
        self._done = False

    def step(self, tick):
        if self._done:
            return
        self._done = True

        people = self.model.people
        count = people.count

        susceptible = np.nonzero(
            people.state[:count] == SEIR.State.SUSCEPTIBLE.value
        )[0]

        n_vaccinate = int(round(self.coverage * len(susceptible)))
        chosen = np.random.choice(susceptible, size=n_vaccinate, replace=False)

        # ← WRONG: TransmissionSE does NOT check this field.
        people.susceptibility[chosen] = 0


# ── CORRECT ───────────────────────────────────────────────────────────────────
# Set state = RECOVERED.  TransmissionSE will skip these agents entirely.

class OneTimeVaccination:
    """Vaccinates `coverage` fraction of all susceptible agents at tick 0.

    TransmissionSE only infects agents where state == SUSCEPTIBLE (int8 value 0).
    Any other state value makes an agent immune.  RECOVERED (value 3) is the
    standard choice for permanent vaccine-derived immunity.

    Node-level S and R counts are updated at tick+1 so that the force-of-
    infection denominator (N = S+E+I+R summed over the patch) stays correct.

    Parameters
    ----------
    model    : LASER Model instance
    coverage : float — fraction of susceptible agents to vaccinate (default 0.80)
    """

    def __init__(self, model, coverage=0.80):
        self.model = model
        self.coverage = coverage
        self._done = False

    def step(self, tick):
        if self._done:
            return
        self._done = True

        people = self.model.people
        nodes  = self.model.nodes
        count  = people.count

        # 1. Find all living susceptible agents.
        susceptible = np.nonzero(
            people.state[:count] == SEIR.State.SUSCEPTIBLE.value
        )[0]

        if len(susceptible) == 0:
            return

        # 2. Draw 80 % without replacement.
        n_vaccinate = int(round(self.coverage * len(susceptible)))
        chosen = np.random.choice(susceptible, size=n_vaccinate, replace=False)

        # 3. THE KEY LINE: set state, not susceptibility.
        #    TransmissionSE checks  `state == SUSCEPTIBLE` (value 0).
        #    Any non-zero state is immune.  RECOVERED (value 3) is correct for
        #    permanent sterilising immunity.
        people.state[chosen] = SEIR.State.RECOVERED.value

        # 4. Update node-level counts so N = S+E+I+R stays consistent.
        #    Components write to tick+1 during their step() call.
        vax_by_node = np.bincount(
            people.nodeid[chosen], minlength=nodes.count
        ).astype(nodes.S.dtype)

        nodes.S[tick + 1] -= vax_by_node
        nodes.R[tick + 1] += vax_by_node

        # Guard against negative counts from rounding or race conditions.
        np.maximum(nodes.S[tick + 1], 0, out=nodes.S[tick + 1])

        total_vaccinated = int(vax_by_node.sum())
        print(
            f"[OneTimeVaccination] tick={tick}: vaccinated {total_vaccinated} "
            f"of {len(susceptible)} susceptible agents "
            f"({100*total_vaccinated/len(susceptible):.1f}%)"
        )


# ── USAGE ─────────────────────────────────────────────────────────────────────
#
# model.components = [
#     SEIR.Susceptible(model),
#     SEIR.Exposed(model, expdurdist, infdurdist),
#     SEIR.Infectious(model, infdurdist),
#     SEIR.Recovered(model),
#     OneTimeVaccination(model, coverage=0.80),   # ← fires once at tick 0
#     SEIR.Transmission(model, expdurdist),
#     BirthsByCBR(model, birthrates=birthrate_map.values, pyramid=pyramid),
#     MortalityByEstimator(model, estimator=survival),
# ]
#
# Component ordering: place OneTimeVaccination BEFORE Transmission so that
# the state change takes effect before the first force-of-infection calculation.


# ── VERIFICATION ──────────────────────────────────────────────────────────────
#
# After model.run(), confirm vaccination actually reduced transmission:
#
#   s_initial = scenario["S"].sum()
#   r_initial = scenario["R"].sum()
#   s_after_vax = int(model.nodes.S[1].sum())   # tick 1, after tick-0 step
#   r_after_vax = int(model.nodes.R[1].sum())
#
#   newly_immune = r_after_vax - r_initial
#   assert newly_immune > 0, "Vaccination had no effect — check state assignment"
#   assert abs(newly_immune / s_initial - 0.80) < 0.01, \
#       f"Expected ~80% coverage, got {newly_immune/s_initial:.1%}"
#
# Or use the built-in check from verification_checks.py:
#   from verification_checks import verify_model_health
#   verify_model_health(model)   # check_vaccination_effect will catch the bug
```

**Why `susceptibility = 0` silently fails:**

`TransmissionSE`'s inner Numba kernel selects infection candidates with a single condition:

```python
if state[i] == SUSCEPTIBLE:   # int8 comparison — that's ALL it checks
    if random() < force[nodeid[i]]:
        state[i] = EXPOSED
        etimer[i] = sample(expdurdist)
```

The `susceptibility` field is never read here. It exists for a different, optional pathway that most models don't use. So an agent with `state=0, susceptibility=0` is **fully susceptible** from the transmission kernel's perspective.

**The three rules for correct vaccination:**

| What | Why |
|------|-----|
| `state = SEIR.State.RECOVERED.value` | The only field `TransmissionSE` checks |
| `nodes.S[tick+1] -= vax_by_node` | Keeps N = S+E+I+R correct for FOI calculation |
| `nodes.R[tick+1] += vax_by_node` | Keeps the compartment sum balanced |
