# Changelog

All notable changes to `fvcom-mesh-tools` are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The
project is **pre-alpha**; breaking changes can land on `main` between
PoC iterations. Once a versioned release tag exists, breaking changes
will only ship with a major bump (Semantic Versioning).

## Unreleased

### Added

- `dem/` subpackage isolates rasterio / netCDF4 / pyproj helpers
  behind the `[dem]` extra:
  - `dem.bbox.read(path) -> (minx, miny, maxx, maxy)` for raster
    bounds.
  - `dem.subset.to_geotiff(src, dst, bbox, ...)` for clipping a DEM
    to a lon/lat bounding box (rasterio path + lon/lat-NetCDF path).
  - `dem.interp.at_points(dem, points, method=...)` for sampling a
    DEM at mesh-node coordinates (bilinear / nearest, EPSG:4326).
- `mesh_engine.ocsmesh` adapter. The
  `mesh_engine.build("ocsmesh", ...)` dispatch path that previously
  resolved to a non-existent module now works; both engines honour
  the `(points, cells) -> EPSG:4326` contract.
- `pyproject.toml` extras layered by concern: `[io-vector]`, `[dem]`,
  `[oceanmesh]`, `[ocsmesh]`, `[viz]`, `[all]`. `[oceanmesh]` and
  `[ocsmesh]` self-reference `[dem,io-vector]`.

### Changed

- **BREAKING**: `fvcom_mesh_tools.io.subset_dem_to_geotiff` is renamed
  and moved to `fvcom_mesh_tools.dem.subset.to_geotiff`.
- **BREAKING**:
  `fvcom_mesh_tools.mesh_engine.depth.interpolate_dem_at_points` is
  renamed and moved to `fvcom_mesh_tools.dem.interp.at_points`.
- **BREAKING**: `fmesh-buildmesh --interp-method` no longer accepts
  `spline`. It had been silently mapped to `linear` on the oceanmesh
  path; only `linear` (default) and `nearest` are honoured now.
- `cli/buildmesh.py` no longer carries an inline OCSMesh mesh-
  generation block. Both engines now go through `mesh_engine.build()`
  and depth interpolation runs uniformly through
  `dem.interp.at_points`, regardless of engine choice.
- `pyproject.toml`: `matplotlib` demoted from runtime dependency to
  the `[viz]` extra. Base install pulls `numpy` only.
- `io/__init__.py` no longer eager-imports anything that touches
  rasterio. `import fvcom_mesh_tools.io` now needs only stdlib +
  numpy; shapely / geopandas / pyproj are imported lazily inside the
  helpers that use them.

### Removed

- `fvcom_mesh_tools.io.dem_subset` module (content moved to
  `fvcom_mesh_tools.dem.subset`).
- `fvcom_mesh_tools.mesh_engine.depth` module (content moved to
  `fvcom_mesh_tools.dem.interp`).
- `pyproject.toml` `[gmsh]` standalone extra (folded into
  `[ocsmesh]`, since gmsh is only invoked through OCSMesh in this
  toolkit).

### Migration

API rename:

```python
# Pre-refactor
from fvcom_mesh_tools.io import subset_dem_to_geotiff
from fvcom_mesh_tools.mesh_engine.depth import interpolate_dem_at_points

# Post-refactor
from fvcom_mesh_tools.dem.subset import to_geotiff
from fvcom_mesh_tools.dem.interp import at_points
```

Pip install workflow:

```bash
# Pre-refactor: matplotlib was a base dep; importing
# fvcom_mesh_tools.io eagerly imported rasterio via the now-removed
# io.dem_subset.
pip install -e .

# Post-refactor: layered extras so a base install pulls numpy only.
pip install -e .                 # base + numpy
pip install -e ".[dem]"          # adds DEM I/O (rasterio/netCDF4/pyproj)
pip install -e ".[oceanmesh]"    # adds the default mesh engine + [dem,io-vector]
pip install -e ".[all]"          # everything
```

CLI flag change: drop any `--interp-method spline` argument; use
`linear` (default) or `nearest`.

The conda workflow under `oceanmesh-bench` / `py312test`
(`docs/architecture.md` §6.1) is unaffected; the heavy deps come
from conda-forge as before, and `pip install --no-deps oceanmesh`
adds the GPL-3.0 engine on top.
