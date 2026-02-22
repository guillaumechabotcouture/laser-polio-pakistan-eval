The model ran successfully. Here's a summary of what was built and the results:

**Model Configuration:**
- 4 patches (Eastport 100K, Northton 200K, Centerburg 150K, Westham 80K)
- Pairwise distances: ~73 km between adjacent patches, up to ~220 km end-to-end
- Gravity network row sums all capped at exactly 0.15 (15% max FOI export)
- R0=5, beta=0.5 (= R0 / mean infectious period = 5/10)

**Verification Results:**
- Population grew from 530K to 644K (21.5% over 10yr) — confirms CBR=30 > CDR=10 units are correct
- All compartments non-negative
- All 4 patches had infections (spatial coupling working)
- 470,612 total infections over 10 years (~89K per million per year per patch)

**Dynamics:** The initial epidemic (1% I with 90% S at R0=5) produces a large first-year outbreak that depletes susceptibles. Susceptible fraction then rebuilds slowly through births. This is the expected pattern for a high-R0 disease without ongoing importation.

The complete code is at `eastland_respiratory_model.py` and the diagnostic plot at `eval/outputs/eastland_respiratory_seir.png`.
