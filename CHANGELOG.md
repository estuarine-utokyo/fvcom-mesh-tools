# Changelog

All notable changes to `fvcom-mesh-tools` are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The
project is **pre-alpha**; breaking changes can land on `main` between
PoC iterations. Once a versioned release tag exists, breaking changes
will only ship with a major bump (Semantic Versioning).

## Unreleased

### Added

- `fvcom_mesh_tools.diagnostics` module + `fmesh-mesh-check` CLI for
  detection of inadequate FVCOM meshes, with no repair. Six detectors
  surface defects that are common in narrow water bodies (rivers,
  canals, harbours):
  1. disjoint dual-graph components (isolated wet pools);
  2. dead-end elements (degree-1 in the dual graph, no open-boundary
     edge);
  3. thin elements (all 3 vertices on a boundary);
  4. thin-chain elements (chain of ≥ `--min-thin-chain` adjacent thin
     elements — the 1-cell-wide-channel signature);
  5. over-connected nodes (valence > `--max-nbr-elem`, FVCOM
     `MAX_NBR_ELEM` cap);
  6. open-boundary unreachable elements.
  The CLI emits `<prefix>_summary.txt`, `<prefix>_diag.json` (per-element
  / per-node records with coordinates), and `<prefix>_map.png`. Exit
  code is non-zero if any detector fires, so the command works as a
  CI gate. Validated on the existing Tokyo Bay (PoC #19, #16) and Osaka
  Bay (PoC #20) fort.14 outputs in PoC #24; the over-connected-node
  finding for PoC #16 (440 nodes, max valence 26) was characterised in
  PoC #25 and traced to the OCSMesh+gmsh engine path (oceanmesh engine
  produces only 3 over-connected nodes for the same inputs). PoC #26
  ablation showed gmsh itself produces ~380 over-connected nodes
  (max v=18) on Tokyo Bay+rivers before any post-processing; turning
  off `--refine-min-angle` (longest-edge bisection) is the single
  largest improvement available within the OCSMesh path
  (440 → 313 over-connected, max v 26 → 21), but does not close the
  gap with the oceanmesh engine.
- `fvcom_mesh_tools.mesh_clean` module + `fmesh-mesh-clean` CLI for
  the safe-repair subset of the diagnostics surfaced by PoC #24-#26.
  Phase A keeps dual-graph connected components by size and / or
  open-boundary touch (default: only the largest). Phase B
  iteratively trims degree-1 elements that have no open-boundary edge
  ("spit" terminations of 1-cell channels). Phase C repairs
  1-cell-wide channels: by default (`--thin-chain-mode widen`) it
  inserts a centroid into every thin-chain element so each thin
  triangle becomes three sub-triangles fanning out from a strictly
  interior point, giving two cells across the channel; with
  `--thin-chain-mode delete` the chain is removed entirely. Phase A
  / B / C-delete deletions re-derive boundaries via DEM-bbox
  proximity (matching `fmesh-buildmesh`); Phase C-widen leaves
  boundary topology untouched because the centroid is interior.
  Validated on PoC #19's Tokyo Bay output: 5,496 disjoint elements
  + 628 dead-end elements removed (Phase A + B), 165 thin-chain
  elements widened with 165 interior nodes added (Phase C). The
  cleaned mesh has 1 connected component, 0 dead-ends, 0 unreachable
  elements, 0 thin chains, alpha mean 0.96, min-angle p50 51 deg,
  frac<20° = 0.16%, and a single open-boundary segment.
  Over-connected-node repair is **not** done here; that needs an
  additional structural policy that is still under design.
- `MESH_PNG_DPI = 600` shared default in `fvcom_mesh_tools.plotting`;
  every notebook that writes a mesh visualisation now passes
  `dpi=MESH_PNG_DPI` to `fig.savefig`. Histograms keep the matplotlib
  default. PoC #23's `outputs/23_overlap_mesh.png` is regenerated at
  600 dpi.
- PoC #23 (`notebooks/23_mesh_combine_overlap.py`) validates
  `fmesh-mesh-combine --strategy overlap` on real Tokyo Bay data:
  a coarse outer (hmin=1000 m, NP=4,224) and a fine northern-bay
  inner (hmin=200 m, NP=6,008) merge via
  `ocsmesh.ops.merge_overlapping_meshes` into a single fort.14
  with NP=8,227, NE=13,923, alpha 0.954, frac<20° 0.09 %, no
  flipped triangles. The CLI hooks for `overlap` and `neighbor`
  have existed since PoC #21 but had no end-to-end real-data
  exercise; `overlap` is now covered, `neighbor` still pending.
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
