# Architecture

This document explains how `fvcom-mesh-tools` is composed and gives a
concrete decision tree for which mesh engine to use for which task.
For the broader gap analysis vs. the legacy OceanMesh2D MATLAB
reference, see `python_pipeline_gap_analysis.md`.

## 1. Pipeline shape

`fmesh-buildmesh` is a single-shot pipeline whose stages, in order,
are:

```
DEM raster
  -> [engine] DEM + coastline -> (points, cells)         # mesher-specific
  -> depth interpolation  (DEM -> per-node depths)
  -> CW->CCW winding fix
  -> bbox-based open/land boundary classification
  -> river-mouth ibtype=21 segment injection
  -> [optional] quality pass (edge swap + Laplacian)
  -> [optional] longest-edge bisection refine
  -> first-ring perpfix on the open boundary
  -> fort.14 write
```

Only the first stage is engine-specific; everything from depth
interpolation onward is plain numpy operating on `(points, cells)`
plus a few index arrays. The engine adapter is responsible for
delivering CCW lon/lat coordinates and an integer triangle table.
Everything else is mesher-agnostic.

```
        +------------------+
        |  fmesh-buildmesh |  (CLI)
        +---------+--------+
                  |
                  v
   +------------------------------+
   |  mesh_engine.build(engine,..)|  (dispatch)
   +------------------------------+
        |               |
        v               v
   +---------+   +-----------+
   |oceanmesh|   |  ocsmesh  |
   +---------+   +-----------+
        |               |
        +-------+-------+
                |
                v   (points, cells, lon/lat)
       +-----------------+
       | depth interp    |  rasterio + scipy bilinear
       | CCW fix         |
       | bbox classify   |
       | river inflow    |
       | quality / refine|
       | perpfix         |
       | write_fort14    |
       +-----------------+
                |
                v
              fort.14
```

Multi-mesh composition (`fmesh-mesh-combine`) sits one step later: it
takes two or more `fort.14` files and produces a single `fort.14`,
either by pure concatenation (`disjoint`) or by delegating to OCSMesh's
`ops.combine_mesh` (`overlap`, `neighbor`).

## 2. Engine choice: oceanmesh vs. ocsmesh

`fmesh-buildmesh --engine` picks between two backends.

### When to use `--engine oceanmesh` (default)

- **Final / production meshes.**
- When near-equilateral element shape matters (e.g. CFL-sensitive
  FVCOM scenarios).
- When you have time for a 5-30 minute generation step.
- When you want bit-identical reproducibility — DistMesh is seeded
  via `--om-seed` (default 0).
- When the input coastline is high-quality (MLIT C23, OSM extracted
  detail). oceanmesh feeds the shoreline directly into its size
  function and benefits from detail.

### When to use `--engine ocsmesh`

- **Iteration / exploration phase.** You're tuning bbox, hmin/hmax,
  river-points, etc., and need fast feedback (~40 s vs ~25 min on
  Tokyo Bay).
- When the coastline is sparse and the bathymetric-gradient size
  function would dominate anyway.
- When OCSMesh's polygon area filters (`--min-polygon-area-m2`,
  `--min-island-area-m2`) are the right hammer for the geometry —
  they have no oceanmesh equivalent yet.

### Side-by-side numbers (Tokyo Bay, hmin=200/hmax=5000)

| metric            | `--engine ocsmesh` | `--engine oceanmesh` |
| ----------------- | ------------------ | -------------------- |
| NP                | 19,362             | 31,771               |
| NE                | 27,609             | 53,203               |
| alpha mean        | 0.847              | **0.959**            |
| frac<20°          | 1.13 %             | **0.10 %**           |
| min-angle p50     | 40.7°              | 51.05°               |
| flipped           | 0                  | 0                    |
| perpfix reverts   | 273                | 8                    |
| wall              | **40 s**           | 1717 s               |
| ibtype=21 rivers  | 5                  | 5                    |

Sources: PoCs #16 and #19. Osaka Bay (PoC #17 vs #20) shows the same
qualitative trend with a smaller wall-clock penalty (3.4× rather than
43×) because GSHHS-f L1 has fewer feature-size constraints than
MLIT C23.

### Engine-specific flags

`--engine oceanmesh` accepts:

- `--om-slope-parameter` (default 20) — `bathymetric_gradient_sizing_function`
  parameter; lower = more aggressive size variation across slope.
- `--om-gradation` (default 0.15) — `enforce_mesh_gradation` ratio.
- `--om-max-iter` (default 50) — DistMesh iteration cap.
- `--om-seed` (default 0) — PRNG seed; identical seed + identical
  inputs → bit-identical mesh.
- `--om-no-bathymetric-gradient` — fall back to feature-only sizing
  if the DEM is too coarse for slope-driven size to be meaningful.
- `--om-minimum-area-mult` (default 4.0) — drop inner-shoreline
  features (islands) below `mult * h0**2`. Raise to coalesce noisy
  islets; PoC #22 sweep on Tokyo Bay shows mult=25 brings the
  retained inner count from 39 down to 27 (matching the `--engine
  ocsmesh` path's `--min-island-area-m2 100000` behaviour).

`--engine ocsmesh` accepts the legacy OCSMesh-driven flags:

- `--coast-target-size`, `--coast-expansion-rate` —
  `Hfun.add_feature` controls.
- `--min-polygon-area-m2`, `--min-island-area-m2` — `Geom`
  filtering before meshing.

## 3. Multi-mesh composition: `fmesh-mesh-combine`

`fmesh-mesh-combine in1.14 in2.14 [in3.14 ...] out.14 --strategy ...`
fuses two or more `fort.14` meshes into one. Three strategies cover the
common cases.

### `--strategy disjoint` (default)

Plain concatenation with index renumbering. Boundaries (open, land,
ibtype=21 rivers) carry forward verbatim.

Use when:

- Inputs do not overlap and do not share any nodes (e.g. two different
  basins put into one model).
- You want to preserve the boundary structure of every input exactly.

Pure numpy; no OCSMesh required. PoC #21 used this to combine Tokyo
Bay (PoC #19, NP=31,771) and Osaka Bay (PoC #20, NP=13,811) into a
single 45,582-node fort.14 in 0.4 s.

### `--strategy overlap`

Wraps `ocsmesh.ops.merge_overlapping_meshes`. The first input is the
*background*; subsequent inputs are *foregrounds* that get carved into
the background, with a buffered seam re-meshed to bridge the
resolution gap.

Use when:

- You have a coarse outer mesh and one or more high-resolution inner
  meshes covering parts of the outer domain.
- You want a single fort.14 that respects the inner refinements.

Tunables: `--buffer-size`, `--buffer-domain`, `--min-int-ang`,
`--adjacent-layers`, `--no-clip-final`.

Boundaries are *not* preserved (OCSMesh's `MeshData` has no boundary
structure). Re-classify downstream with `fmesh-buildmesh`-style
post-processing or `omesh14-edit-bdy`.

### `--strategy neighbor`

Wraps `ocsmesh.ops.merge_neighboring_meshes`. KDTree-snaps shared
boundary vertices between two meshes whose edges already coincide
within ~1e-8 deg (~ 1 mm).

Use when:

- You generated two halves of a domain separately (e.g. east half +
  west half) with a shared seam line.
- Edge nodes on the seam already line up exactly between the two
  inputs.

Boundaries are dropped, like `overlap`.

## 4. Decision tree

```
                Need a single regional fort.14?
                          |
                +---------+---------+
                |                   |
            yes (one mesh)      no  (combining inputs)
                |                   |
        Iterating params?       Inputs overlap?
                |                   |
        +-------+-------+   +-------+-------+
        |               |   |               |
       yes             no  no overlap    yes overlap
        |               |   |               |
        v               v   v               v
   --engine ocsmesh   --engine     --strategy disjoint    Edges already
   (40 s)             oceanmesh                            coincide?
                      (default)                            |
                                                  +--------+--------+
                                                  |                 |
                                                 yes               no
                                                  |                 |
                                                  v                 v
                                           --strategy        --strategy
                                            neighbor          overlap
```

## 5. Module map

```
src/fvcom_mesh_tools/
+-- algorithms/         # mesher-agnostic: quality, edge_swap,
|                       # smoothing, refine, perpfix, boundary,
|                       # rivers, signed_areas, ...
+-- io/                 # fort.14 read/write, DEM subset,
|                       # coastline loader, geom filter,
|                       # river-points loader
+-- mesh_engine/        # generation backends
|   +-- __init__.py     # build(engine=...) dispatcher
|   +-- oceanmesh.py    # primary - DistMesh
|   +-- depth.py        # bilinear DEM sampling for non-OCSMesh path
|   `-- (ocsmesh path inline in cli/buildmesh.py for now)
+-- mesh_compose/       # multi-mesh stitching
|   +-- __init__.py     # combine(strategy=...) dispatcher
|   +-- disjoint.py     # numpy concat with index offsets
|   +-- convert.py      # Fort14Mesh <-> ocsmesh MeshData
|   `-- overlap.py      # OCSMesh-backed merge_overlapping/neighboring
+-- vis/, plot/         # visualization helpers
+-- scan/, mesh/        # legacy oceanmesh-tools features
`-- cli/                # console entry points
    +-- buildmesh.py        - fmesh-buildmesh
    +-- perpfix.py          - fmesh-perpfix
    +-- subset_dem.py       - fmesh-subset-dem
    `-- meshcombine.py      - fmesh-mesh-combine
```

## 6. Environment model

Two reproducible conda environments cover the toolkit:

* `py312test` -- the original env. Has OCSMesh + gmsh; lacks
  oceanmesh. Sufficient for `--engine ocsmesh` and unit tests.
* `oceanmesh-bench` -- the unified env. Has both oceanmesh and
  OCSMesh, plus all CLIs. Recommended for running PoCs #18-#21 and
  any production builds.

Setup scripts: `notebooks/18_setup_oceanmesh_env.pjsub` for initial
creation, `notebooks/18_setup_gdal_netcdf.pjsub` for the GDAL netCDF
driver, `notebooks/19_setup_env_extend.pjsub` for fvcom-mesh-tools
editable install, `notebooks/21_setup_ocsmesh_in_omsh.pjsub` for the
OCSMesh top-up.

## 7. Licensing

| Component                | License                     |
| ------------------------ | --------------------------- |
| `fvcom-mesh-tools` code  | Apache-2.0                  |
| `oceanmesh`              | GPL-3.0-or-later            |
| OCSMesh                  | CC0-1.0                     |
| gmsh (called by OCSMesh) | GPL-2.0-or-later (runtime)  |

Importing `oceanmesh` directly from this Apache-2.0 toolkit means the
**combined work** is subject to GPL-3.0 obligations when redistributed.
Source distributions of `fvcom-mesh-tools` itself remain Apache-2.0;
the GPL footprint kicks in only when you ship `oceanmesh` alongside.
The `--engine ocsmesh` path stays Apache-friendly (CC0 + GPL-only-runtime
through gmsh).

See `THIRD_PARTY_NOTICES.md` for the full attribution text.
