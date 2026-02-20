# A/B Test Prompts: LASER Skill Evaluation

## Protocol

For each prompt below:
1. Open a NEW Claude Code session in the project directory
2. Run session A (with skill): `cd laser-polio-pakistan-eval && claude`
3. Run session B (without skill): `cd laser-polio-pakistan-eval && claude --deny "Skill(laser-spatial-disease-modeling)"`
4. Paste the EXACT same prompt into each session
5. Let Claude generate its response WITHOUT follow-up questions (one-shot)
6. Save the full output (code + explanation) to `eval/outputs/{with,without}-skill/prompt-N.md`

## Prompts

### Prompt 1: Basic Model Setup
**Tests:** LASER API knowledge, correct imports, component signatures

```
Using the LASER framework (laser-generic package), set up a basic spatial SEIR
model for poliovirus transmission across 10 patches representing districts in
Pakistan. Each patch has a population of 100,000. Use a transmission rate
corresponding to R0≈6, a latent period of 3 days, and an infectious period of
28 days. Initialize with 95% recovered (immune) and 5 infectious per patch.
Write the complete Python code to configure and run a 10-year simulation.
Do not install any packages.
```

### Prompt 2: Gravity Network + Seasonal Forcing
**Tests:** Migration model setup, seasonal transmission, spatial coupling

```
Extend the basic LASER polio model to include:
1. A gravity-model migration network between patches, with distance decay
   exponent c=1.5 and destination population exponent b=0.5
2. Seasonal forcing that peaks during July-October (monsoon season) with
   amplitude 1.3x baseline and troughs during December-March at 0.7x baseline
3. Row-normalize the network so no patch exports more than 15% of its
   force of infection

Write complete Python code using the LASER framework. Assume patches are
arranged in a line 50km apart. Do not install any packages.
```

### Prompt 3: Vaccination Components
**Tests:** Immunization API, domain adaptation (OPV/IPV), custom component design

```
Add vaccination to the LASER polio model:
1. Routine immunization at 6 weeks of age with 80% coverage (OPV)
2. Supplementary Immunization Activities (SIAs) every 6 months targeting
   children aged 0-5 years with 90% coverage
3. The model should track vaccinated individuals separately from naturally
   recovered, since OPV provides mucosal immunity while natural infection
   provides stronger, longer-lasting immunity

Show the complete component setup using LASER's immunization classes where
possible, and custom components where needed. Do not install any packages.
```

### Prompt 4: Calibration Framework
**Tests:** Calibration metric design, adaptation from measles CCS to polio AFP

```
Design a calibration framework for the LASER polio model that:
1. Compares simulated incidence to AFP (Acute Flaccid Paralysis) surveillance
   data, accounting for the ~1:200 paralysis-to-infection ratio
2. Uses a fitness metric based on the proportion of zero-incidence weeks per
   district (similar to critical community size analysis)
3. Samples from these parameter ranges:
   - beta: U(0.15, 0.30)
   - gravity_k: 10^U(-3, -1)
   - seasonal_amplitude: U(1.0, 1.5)
4. Ranks simulations by goodness of fit

Write the calibration loop and scoring functions. Use LASER framework
conventions. Do not install any packages.
```

### Prompt 5: Full Integration (Hardest)
**Tests:** End-to-end model assembly, all components working together

```
Build a complete spatial polio transmission model for Pakistan using the LASER
framework with these specifications:
- 20 districts with heterogeneous populations (range 50k-2M)
- SEIR dynamics with R0≈6, 3-day latent period, 28-day infectious period
- Gravity-model spatial coupling
- Monsoon-season transmission forcing (peak Jul-Oct)
- OPV routine immunization at 80% coverage
- SIA campaigns every 6 months
- Importation of 3 infections every 60 days for the first 5 years
- 20-year simulation with 10-year burn-in
- Output: weekly incidence per district, proportion of zero-incidence weeks

Write complete, runnable Python code. Do not install any packages.
```
