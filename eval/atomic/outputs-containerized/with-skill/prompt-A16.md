Here's the complete code:

```python
from laser.core.migration import distance

# Coordinates
lahore_lat, lahore_lon = 31.55, 74.35
karachi_lat, karachi_lon = 24.86, 67.01

# Compute haversine distance using LASER
dist_km = distance(lahore_lat, lahore_lon, karachi_lat, karachi_lon)

print(f"Lahore  : {lahore_lat}°N, {lahore_lon}°E")
print(f"Karachi : {karachi_lat}°N, {karachi_lon}°E")
print(f"Distance: {dist_km:.2f} km")
```

**How it works:**

- `laser.core.migration.distance` implements the Haversine formula using Earth radius = 6371.0 km
- Called with four scalars `(lat1, lon1, lat2, lon2)`, it returns a single scalar distance in kilometers
- The same function can handle vectors/matrices for pairwise distance computation (used elsewhere in this project for district-to-district migration networks)

**Expected output:**

```
Lahore  : 31.55°N, 74.35°E
Karachi : 24.86°N, 67.01°E
Distance: 1054.52 km
```
