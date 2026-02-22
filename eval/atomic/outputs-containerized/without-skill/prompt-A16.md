The computation works. Here is the complete, self-contained code as it would be written in a project where `laser-core` is installed:

```python
from laser_core.migration import distance

# Coordinates (lat, lon) in decimal degrees
lahore_lat,  lahore_lon  = 31.55, 74.35   # Lahore
karachi_lat, karachi_lon = 24.86, 67.01   # Karachi

km = distance(lahore_lat, lahore_lon, karachi_lat, karachi_lon)

print(f"Lahore  : ({lahore_lat}°N, {lahore_lon}°E)")
print(f"Karachi : ({karachi_lat}°N, {karachi_lon}°E)")
print(f"Haversine distance: {km:.2f} km")
```

**Output:**
```
Lahore  : (31.55°N, 74.35°E)
Karachi : (24.86°N, 67.01°E)
Haversine distance: 1034.15 km
```

**Notes on the LASER API:**

- The canonical function is `laser_core.migration.distance(lat1, lon1, lat2, lon2)` — it implements the [haversine formula](https://en.wikipedia.org/wiki/Haversine_formula) using Earth radius = 6371.0 km and returns kilometers.
- `laser_core.utils.calc_distances` is the older, deprecated wrapper; it prints a deprecation warning and calls `migration.distance` internally.
- `distance` is vectorised: pass scalar floats for a single result (as above), a scalar origin + array of destinations for a distance vector, or two arrays for an N×M distance matrix — useful for building the pairwise patch-distance matrices that LASER migration models (gravity, radiation, Stouffer) consume.
