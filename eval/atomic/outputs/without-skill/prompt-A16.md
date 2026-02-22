Here's the complete code and result:

```python
from laser.core.migration import distance

# Lahore: 31.55°N, 74.35°E
# Karachi: 24.86°N, 67.01°E
km = distance(31.55, 74.35, 24.86, 67.01)
print(f"Haversine distance from Lahore to Karachi: {km:.2f} km")
```

**Result: 1034.15 km**

The `distance()` function in `laser.core.migration` implements the Haversine formula using an Earth radius of 6371.0 km. It accepts scalar lat/lon pairs — `distance(lat1, lon1, lat2, lon2)` — and returns the great-circle distance in kilometers. It also supports array inputs for computing pairwise distance matrices, which is how LASER uses it for migration network calculations.
