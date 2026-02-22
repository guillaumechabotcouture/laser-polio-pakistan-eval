It works. Here's the script and explanation:

## `age_pyramid_sample.py`

**Key concepts:**

1. **`AliasedDistribution(counts)`** takes an array of weights (counts or proportions) — one per bin. It builds a Vose alias table internally for O(1) sampling.

2. **`dist.sample(count=n)`** returns **bin indices** (0–4 in this case), not ages. You must map indices to actual ages yourself.

3. The script draws a uniform random age within each selected bin to produce realistic integer ages.

**How the bins map:**

| Bin | Age range | Weight |
|-----|-----------|--------|
| 0   | 0–4       | 30%    |
| 1   | 5–14      | 25%    |
| 2   | 15–29     | 20%    |
| 3   | 30–49     | 15%    |
| 4   | 50+       | 10%    |

The sampled proportions (27.3%, 28.0%, 19.6%, 14.6%, 10.5%) closely match the targets, as expected from 1000 draws with stochastic variation.

Sources:
- [LASER Population Pyramids Documentation](https://docs.idmod.org/projects/laser/en/latest/pyramids.html)
- [laser_core.demographics API Reference](https://docs.idmod.org/projects/laser/en/v0.4.1/source/laser_core.demographics.html)
- [LASER GitHub Repository](https://github.com/InstituteforDiseaseModeling/laser)
