Here's how to compute a gravity-based migration network for 4 patches and assign it to a LASER model:

```python
import numpy as np
import polars as pl
from laser.core.migration import gravity, distance

# --- 1. Define 4 patches with populations and coordinates ---
populations = np.array([100_000, 50_000, 20_000, 75_000], dtype=np.float64)
lats = np.array([0.0, 0.5, 1.0, 0.5])
lons = np.array([0.0, 0.5, 0.0, 1.0])

# --- 2. Compute pairwise distance matrix (km, via Haversine) ---
n = len(populations)
distances = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j:
            distances[i, j] = distance(lats[i], lons[i], lats[j], lons[j])

# --- 3. Compute the raw gravity network ---
#   formula: network[i,j] = k * pop_i^a * pop_j^b / dist_ij^c
network = gravity(
    pops=populations,
    distances=distances,
    k=0.01,   # scaling constant (average trip probability)
    a=1.0,    # exponent on origin population
    b=1.0,    # exponent on destination population
    c=1.5,    # distance decay exponent
)

# --- 4. Normalize into a mixing matrix (rows sum to 1) ---
#   Diagonal = probability of staying; off-diagonal = probability of traveling
row_sums = network.sum(axis=1, keepdims=True)
mixing_matrix = network / row_sums          # normalize off-diagonal
np.fill_diagonal(mixing_matrix, 0)
outflow = mixing_matrix.sum(axis=1)
np.fill_diagonal(mixing_matrix, 1.0 - outflow)  # staying = 1 - leaving

print("Gravity network (raw flows):\n", network)
print("\nMixing matrix (rows sum to 1):\n", mixing_matrix)
```

### Assigning to a model via `GravityMixing`

The idiomatic LASER way uses the `GravityMixing` wrapper, which handles the distance computation and normalization automatically:

```python
from laser.measles.compartmental.model import CompartmentalModel
from laser.measles.compartmental.params import CompartmentalParams
from laser.measles.mixing.gravity import GravityMixing, GravityParams

# Scenario DataFrame — required columns: id, pop, lat, lon, mcv1
scenario = pl.DataFrame({
    "id":   ["patch_0", "patch_1", "patch_2", "patch_3"],
    "pop":  [100_000, 50_000, 20_000, 75_000],
    "lat":  [0.0, 0.5, 1.0, 0.5],
    "lon":  [0.0, 0.5, 0.0, 1.0],
    "mcv1": [0.85, 0.90, 0.80, 0.88],
})

# Create the model
params = CompartmentalParams(seed=42, num_ticks=365)
model = CompartmentalModel(scenario, params, name="4patch")

# Build the gravity mixer and inspect matrices
gravity_params = GravityParams(a=1.0, b=1.0, c=1.5, k=0.01)
mixer = GravityMixing(scenario=scenario, params=gravity_params)

# These are the computed matrices (cached on first access):
print(mixer.migration_matrix)  # raw gravity flows, diagonal = 0
print(mixer.mixing_matrix)     # row-normalized, diagonal = stay probability
```

The mixer is then passed to the infection component, which is how LASER connects the network to the model's transmission dynamics:

```python
from laser.measles.compartmental.components.process_infection import (
    InfectionProcess, InfectionParams,
)

infection_params = InfectionParams(beta=1.0, mixer=mixer)
# InfectionProcess uses mixer.mixing_matrix during force-of-infection calculation:
#   lambda_i = (beta * seasonal * prevalence) @ mixer.mixing_matrix
```

**Key points:**
- `gravity()` in `laser.core.migration` is the low-level function — it returns a raw NumPy array
- `GravityMixing` is the high-level wrapper — it auto-computes distances from lat/lon, normalizes so average outflow equals `k`, and caches the result
- The diagonal of the mixing matrix represents the fraction staying in each patch; off-diagonals are inter-patch travel probabilities
- The network enters the model through the `mixer` parameter on the infection component, not a standalone `model.network` attribute
