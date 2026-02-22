The script runs cleanly. Here's what it does:

**Patch setup** — Six East African cities with populations ranging from 10,000 (Bujumbura) to 100,000 (Kampala), each placed at real-world coordinates.

**SEIR initialization** — For each patch:
- **S** = 98% of population
- **E** = 1% of population
- **I** = 1% of population
- **R** = remainder (absorbs any rounding so `S+E+I+R == population` exactly)

**GeoDataFrame** — Built with `EPSG:4326` (WGS84) CRS and `Point` geometries, ready for spatial operations like distance-based connectivity kernels or plotting with `.plot()`.

This GeoDataFrame serves as the scenario input a LASER spatial model would consume — you'd feed it into a model that computes patch-to-patch transmission based on the geometry and runs SEIR dynamics on each patch at each timestep.
