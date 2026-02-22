The simulation ran successfully in about 5 seconds. Here's a summary of the results:

**Model Configuration:**
- 4 patches in the Republic of Eastland (530,000 total population)
- 10-year simulation (3,650 days)
- R0=5.0, beta=0.50, latent period mean=4 days, infectious period mean=10 days
- Gravity network (k=0.01, a=1, b=1, c=1.5) with row normalization at 0.15
- Seasonal forcing: cosine cycle peaking at 1.3x in winter, troughing at 0.7x in summer
- Vital dynamics: CBR=30/1000/yr, CDR=10/1000/yr

**Results:**
- Population grew from 530,000 to 643,920 (net +113,920 from 173k births vs 60k deaths)
- Initial epidemic peaked at day 17 with ~32,600 new infections/day
- Total infections over 10 years: 475,836
- Final state shows endemic equilibrium with ~166k susceptible and ~478k recovered
- The disease burned through the population quickly (high R0=5) and subsequent seasonal waves replenished susceptibles through births

The complete code is in `eastland_seir.py`.
