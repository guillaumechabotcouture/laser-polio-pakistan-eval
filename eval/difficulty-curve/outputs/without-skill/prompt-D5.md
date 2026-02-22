The model runs successfully. Here's a summary of the complete code in `eastland_seir.py`:

**Model structure:**

| Component | Details |
|---|---|
| **Patches** | 4 nodes: 100k, 200k, 150k, 80k population |
| **Coordinates** | (30-31.5N, 50-51.5E), ~73 km between neighbors |
| **Disease** | SEIR, beta=0.5 (R0~5), latent gamma(4,1), infectious gamma(5,2) |
| **Initial state** | 90% S, 0% E, 1% I, 9% R |
| **Spatial coupling** | Gravity network (k=0.01, a=1, b=1, c=1.5), row-normalized |
| **Vital dynamics** | BirthsByCBR (CBR=30), MortalityByCDR (CDR=10) |
| **Duration** | 10 years (3,650 daily ticks) |

**Key results after 10 years:**
- Population grew from 530,000 to 644,604 (net +21.6% from CBR=30 vs CDR=10)
- 173,424 births, 58,820 deaths
- 470,392 total infections across 4 patches
- Epidemic burned through, leaving all patches with 0 active E/I at year 10
