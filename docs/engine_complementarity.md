# oceanmesh vs. ocsmesh: capability map and complementarity

This note consolidates the empirical and source-level investigation of
the two mesh-generation libraries the toolkit depends on, so the
project can decide which pieces of each to keep, drop, or wrap.

It supersedes informal notes scattered across `docs/architecture.md`
§2, the PoC #16/#18/#19/#25/#26/#28/#30 reports, and the CHANGELOG.

> **Note on terminology.** Throughout this document, "oceanmesh"
> refers to the [`oceanmesh`](https://github.com/CHLNDDEV/oceanmesh)
> Python port of OceanMesh2D (DistMesh-based, GPL-3). "ocsmesh"
> refers to NOAA-OWP's [`ocsmesh`](https://github.com/noaa-ocs-modeling/OCSMesh)
> (CC0, wraps gmsh / Triangle). Both produce ADCIRC-flavoured fort.14.

## 1. Bottom-line summary

| Question | Answer |
| --- | --- |
| Which produces higher-quality output meshes today? | **oceanmesh** by a clear margin (alpha 0.96 vs. 0.85 on Tokyo Bay; max valence 9 vs. 26). PoC #16 vs. #19. |
| Can we drop ocsmesh entirely? | **No.** Several capabilities have no oceanmesh equivalent (next section). |
| Can we drop *gmsh* (the part of ocsmesh that misbehaves)? | **No, not without losing the buildmesh path.** ocsmesh's Triangle backend rejects varying sizing (PoC #30); gmsh is the only ocsmesh engine that consumes a `Hfun(raster)`. |
| Can we drop ocsmesh's *MeshDriver* and keep its libraries? | **Yes.** `ocsmesh.ops.combine_mesh`, `ocsmesh.utils`, and `ocsmesh.Raster` are gmsh-independent and individually useful. |

The strategic conclusion is therefore "**use oceanmesh for build,
ocsmesh as a library** — not as an engine." The next sections back up
that recommendation with the underlying inventory.

## 2. Empirical quality benchmark (Tokyo Bay, hmin=200 / hmax=5000)

Reproduced from `docs/architecture.md` §2.97-110 and the PoC log:

| metric            | `--engine ocsmesh` (gmsh) | `--engine oceanmesh` |
| ----------------- | ------------------------- | -------------------- |
| NP                | 19,362                    | 31,771               |
| NE                | 27,609                    | 53,203               |
| alpha mean        | 0.847                     | **0.959**            |
| frac<20°          | 1.13 %                    | **0.10 %**           |
| min-angle p50     | 40.7°                     | 51.05°               |
| max valence       | 26                        | **9**                |
| n over-connected  | 440                       | **3**                |
| flipped           | 0                         | 0                    |
| wall              | **40 s**                  | 1717 s               |

PoC #25/#26 traced the over-connected anomaly: ~380 of the 440
over-connected nodes are produced by gmsh itself (before any
post-processing). Switching the OCSMesh-side `--refine-min-angle` off
reduces the count to 313 but does not close the gap — gmsh's local
topology is the primary source.

PoC #30 (2026-05-10) confirmed that ocsmesh ships a Triangle backend
but it raises `NotImplementedError("Varying sizing is not supported
for Triangle engine!")` whenever `Hfun(raster)` is supplied. Triangle
is therefore *not* a drop-in replacement for gmsh in our build path.

## 3. Capability map: who has what

Public APIs split cleanly into three buckets.

### 3.1 oceanmesh-only

These are the load-bearing reasons to use `--engine oceanmesh` for
production output.

| Capability | API | Notes |
| --- | --- | --- |
| Force-based DistMesh (high-quality output) | `oceanmesh.generate_mesh` | Pure Python + C++ Delaunay extension; no gmsh. Default in `fmesh-buildmesh`. |
| Multiscale nested grids in one shot | `oceanmesh.generate_multiscale_mesh`, `oceanmesh.multiscale_signed_distance_function` | Validates nesting, blends gradation; **must be planned upfront** (cannot stitch independently-generated meshes). |
| Hamilton–Jacobi gradation | `oceanmesh.enforce_mesh_gradation` | C++ extension `_HamiltonJacobi`. ocsmesh has no analogue; its gradation is implicit in gmsh / Triangle. |
| Wavelength sizing (shallow water) | `oceanmesh.wavelength_sizing_function` | Useful for FVCOM tide-driven cases. |
| Bathymetric gradient with Rossby filter | `oceanmesh.bathymetric_gradient_sizing_function` | More physical than ocsmesh's plain gradient; uses Coriolis radius. |
| Distance-from-line / point sizing | `oceanmesh.distance_sizing_from_line_function`, `_from_point_function` | Cleaner river/jet refinement than ocsmesh's `add_feature`. |
| Determinism via PRNG seed | `--om-seed` (default 0) | ocsmesh's gmsh path is non-deterministic (`docs/python_pipeline_gap_analysis.md` §2.7). |
| Boundary section classification by depth | `oceanmesh.identify_ocean_boundary_sections` | Supplements our own bbox-tol classifier. |
| fort.14 writer | `oceanmesh.write_to_fort14` | Already used in `fmesh-buildmesh`. |

### 3.2 ocsmesh-only

These are the reasons not to drop ocsmesh as a library, even if we
stop using its `MeshDriver`.

#### 3.2.1 Sizing primitives that oceanmesh lacks

| Primitive | API | Why valuable |
| --- | --- | --- |
| Courant-number constraint | `Hfun(...).add_courant_num_constraint(target_courant, dt, ...)` | Directly enforces FVCOM CFL stability via `dt * sqrt(g h) / dx ≤ C`. oceanmesh has no Courant-aware sizing. |
| Subtidal flow limiter | `Hfun(...).add_subtidal_flow_limiter(...)` | Applies extra refinement in shoaling zones using the topographic gradient profile. |
| Elevation-bounded constraints | `Hfun(...).add_topo_bound_constraint(lower, upper, hmin, hmax)` | Force a fixed cell-size band to a bathymetric range (e.g. 0 m to −5 m gets `dx ≤ 100 m`). |
| Arbitrary depth-function constraint | `Hfun(...).add_topo_func_constraint(func, expansion_rate)` | Generic API; can express user-defined size laws. |
| Region constraint | `Hfun(...).add_region_constraint(gdf, hmin, hmax)` | Per-GeoDataFrame-feature size override (harbour basins etc.). |
| Refinement patch | `Hfun(...).add_patch(geom, target_size)` | Local fixed-size override with optional falloff. |
| Channel auto-detection | `Hfun(...).add_channel(...)` and `Raster.get_channels()` | Heuristic narrow-domain refinement; relates to detector 6 (PoC #28). |
| Multiple-Hfun composition | `HfunCollector` | Combine independently-built size functions with priority logic. |

> All of these registrations themselves are gmsh-free. **The cost is
> paid at apply time**: `Hfun(raster).meshdata()` builds a
> background-size mesh via `get_mesh_engine('gmsh')`, which is the
> bind we hit in PoC #30. Any port of these primitives onto
> oceanmesh's structured-grid sizing model would only need the
> registration logic, not the apply path.

#### 3.2.2 Mesh composition (multi-mesh stitching)

| Function | Capability | Notes |
| --- | --- | --- |
| `ocsmesh.ops.combine_mesh.merge_overlapping_meshes` | Carve foreground meshes into a background, re-mesh the seam buffer | **Internally uses Triangle**, not gmsh. PoC #23 validated overlap on Tokyo Bay. |
| `ocsmesh.ops.combine_mesh.merge_neighboring_meshes` | Stitch meshes whose seam edges already coincide within tolerance (KDTree snap) | Triangle-based. Used by `--strategy neighbor`. |
| `ocsmesh.ops.combine_mesh.clip_mesh_by_mesh` | Clip a mesh by another mesh's hull | Useful for nested model pre-processing. |
| `ocsmesh.ops.river_mesh.quadrangulate_rivermapper_arcs` | Generate quad-meshed river network from line arcs | gmsh-dependent; not currently used. |

oceanmesh has **no equivalent** to overlap/neighbor stitching: its
`generate_multiscale_mesh` requires all regions to be planned together
in a single call. For external models (regional ocean outputs that
need to be carved into a coastal mesh) ocsmesh.ops is the only path.

#### 3.2.3 Mesh post-processing utilities

| Function | Capability | gmsh-free? |
| --- | --- | --- |
| `ocsmesh.utils.cleanup_skewed_el` | Remove or improve high-skew triangles | Yes |
| `ocsmesh.utils.cleanup_isolates` | Drop unused nodes | Yes |
| `ocsmesh.utils.cleanup_duplicates` | Merge duplicate nodes/elements | Yes |
| `ocsmesh.utils.cleanup_pinched_nodes` | Remove zero-length-edge nodes | Yes |
| `ocsmesh.utils.cleanup_folded_bound_el` | Remove folded boundary triangles | Yes |
| `ocsmesh.utils.repartition_features` | Redistribute boundary nodes along feature lines | Yes |
| `ocsmesh.utils.interpolate_euclidean_mesh_to_euclidean_mesh` | Transfer nodal values between meshes | Yes |
| `ocsmesh.utils.get_polygon_channels` | Extract narrow waterway geometry | Yes |
| `ocsmesh.utils.sieve` / `finalize_mesh` | Drop small isolated regions | Yes |
| `ocsmesh.utils.calc_el_angles`, `calculate_edge_lengths` | Standard quality metrics | Yes |

oceanmesh's post-processing is narrower: `mesh_clean`, `fix_mesh`,
`laplacian2`, `delete_boundary_faces`, `delete_interior_faces`,
`delete_exterior_faces`, `delete_faces_connected_to_one_face`. No
equivalent of `cleanup_skewed_el` (sliver-aware), `repartition_features`
(boundary node redistribution), or `interpolate_*` (mesh-to-mesh
field transfer).

#### 3.2.4 Raster / DEM helpers (`ocsmesh.Raster`)

| Method | Capability | gmsh-free? |
| --- | --- | --- |
| `Raster.clip(polygon)` | Polygon clip | Yes |
| `Raster.warp(crs)`, `Raster.resample(dx)` | Re-project / re-grid | Yes |
| `Raster.fill_nodata()` | Interpolate NaN gaps | Yes |
| `Raster.gaussian_filter()`, `Raster.average_filter()` | DEM smoothing | Yes |
| `Raster.mask(polygon)` | Polygon mask | Yes |
| `Raster.get_contour(z)` | Extract contour lines at level | Yes |
| `Raster.get_channels(...)` | Auto-detect narrow regions | Yes |
| `Raster.sample(points)` | Query at arbitrary points | Yes |

oceanmesh's `DEM`/`Grid` classes do CRS transforms and point
interpolation but lack `clip`/`fill_nodata`/`gaussian_filter`/
`get_channels`. For DEM preprocessing (which `fmesh-subset-dem`
already partially uses), `ocsmesh.Raster` is the richer tool.

### 3.3 Both have it

Where the libraries overlap, the choice is a quality/usability trade,
not a feature gap.

| Capability | oceanmesh | ocsmesh |
| --- | --- | --- |
| `distance_sizing_function` | Yes | Yes (via `add_feature` + distance) |
| `feature_sizing_function` | Yes | Yes (via `add_contour` + `add_feature`) |
| Bathymetric gradient sizing | Yes (with Rossby filter) | Yes (`add_subtidal_flow_limiter`) |
| Shoreline-driven domain | `Shoreline` (.shp only) | `Geom(raster)` (any DEM) + `MultiPolygonGeom` |
| Island-area filter | `Shoreline.minimum_area_mult` | `Geom.get_multipolygon` + custom filter (we already implemented `filter_multipolygon_by_area` in `fvcom_mesh_tools.io`) |
| fort.14 / .grd / .2dm I/O | fort.14 only (read/write) | All three formats (`mesh.write(format=…)`) |
| Boundary classification | `identify_ocean_boundary_sections` (depth-based) | `get_boundary_segments`/`_edges` (topological only) |

## 4. Compatibility-table answers to common questions

> *"If we want X, which library?"*

| Goal | Pick |
| --- | --- |
| Final production fort.14 with high alpha | **oceanmesh** |
| Fast draft (~40 s vs. ~25 min) | ocsmesh+gmsh today; in future, oceanmesh with reduced `--om-max-iter` (untested but plausible) |
| Stitch independently-generated meshes (overlap or neighbour) | **ocsmesh** (`ocsmesh.ops.combine_mesh`); no oceanmesh equivalent |
| CFL-aware sizing | **ocsmesh** (`add_courant_num_constraint`); no oceanmesh equivalent |
| Sliver-aware post-processing | **ocsmesh** (`cleanup_skewed_el`); oceanmesh has only `delete_boundary_faces` |
| Mesh-to-mesh field transfer | **ocsmesh** (`utils.interpolate_*`); not in oceanmesh |
| DEM preprocessing (fill, smooth, clip, channel detect) | **ocsmesh.Raster** |
| Boundary classification by depth | **oceanmesh** (`identify_ocean_boundary_sections`) |
| Determinism (bit-identical reruns) | **oceanmesh** (`--om-seed`) |

## 5. Recommended division of labour

The project's actual workflow shape, based on this map:

```
┌──────────────────────────────────────────────────────────────────┐
│  fmesh-subset-dem        ── ocsmesh.Raster (clip / fill / smooth) │
│              │                                                    │
│              ▼                                                    │
│  fmesh-buildmesh         ── oceanmesh.generate_mesh   ◄ build     │
│              │             (default; high quality)                │
│              ▼                                                    │
│  fmesh-mesh-combine      ── ocsmesh.ops.combine_mesh  ◄ stitch    │
│   (overlap / neighbor)     (Triangle-based, gmsh-free)            │
│              │                                                    │
│              ▼                                                    │
│  fmesh-mesh-check        ── fvcom_mesh_tools.diagnostics ◄ detect │
│              │                                                    │
│              ▼                                                    │
│  fmesh-mesh-clean        ── fvcom_mesh_tools.mesh_clean (Phase A-E)│
│   Phase F (planned)      ── ocsmesh.utils.cleanup_skewed_el or   │
│                              oceanmesh.laplacian2  ◄ improve      │
└──────────────────────────────────────────────────────────────────┘
```

In other words:

* **`fmesh-buildmesh`**: oceanmesh becomes the only build path.
  `--engine ocsmesh` is deprecated and slated for removal in a
  future release. PoC #36 (``--om-max-iter`` sweep) closed the
  draft-turnaround case: `--om-max-iter 10` produces alpha 0.943
  on Tokyo Bay in ~7 min — still well above ocsmesh+gmsh's 0.847
  in ~40 s — so the quality gap dominates the wall-clock advantage
  and there is no remaining workflow that prefers the older path.
  PoC #30 had earlier ruled out a Triangle-without-gmsh escape
  hatch (`NotImplementedError("Varying sizing is not supported for
  Triangle engine!")`).
* **`fmesh-mesh-combine`**: stays on `ocsmesh.ops.combine_mesh`. This
  is gmsh-free and has no oceanmesh substitute.
* **`fmesh-mesh-clean`**: keeps its own pure-Python phases A–E, but
  Phase F (sliver clean) and a hypothetical Phase G (Laplacian
  smoothing) should evaluate `ocsmesh.utils.cleanup_skewed_el` and
  `oceanmesh.laplacian2` respectively rather than reinvent.
* **`fmesh-subset-dem`**: continues to use `ocsmesh.Raster` for the
  heavy lifting (already does so).
* **Sizing primitives** (`add_courant_num_constraint` etc.) are
  *not* automatically usable from the oceanmesh build because the
  apply path goes through gmsh. Porting them onto oceanmesh's
  structured-grid sizing model is a deferred task; the registration
  logic is gmsh-free and reusable but the application logic is
  not.

## 6. Decisions to confirm

The following are recommendations, not yet decided:

1. **Deprecate `--engine ocsmesh`** in `fmesh-buildmesh`. Keep the
   code path one release for migration, then remove. Driver: PoC #30
   showed the Triangle escape hatch does not exist in ocsmesh today.
2. **Keep `[ocsmesh]` extra in `pyproject.toml`** (renamed to e.g.
   `[ocsmesh-utils]` if we want to signal "library-only"). Drop
   `gmsh` from the extra once the build path is gone.
3. **Phase F = `ocsmesh.utils.cleanup_skewed_el` wrapper** instead of
   a from-scratch sliver detector. Saves implementation, inherits
   tested behaviour.
4. **Phase G = `oceanmesh.laplacian2` wrapper** with a per-step
   guard against boundary-node motion. Only enable on demand.
5. **Port the most valuable sizing primitive — `add_courant_num_constraint`**
   — onto oceanmesh's `Grid`/`compose_h0` structure. **Done** — see
   `fvcom_mesh_tools.mesh_engine.oceanmesh.courant_sizing_function`
   and the `--om-courant-sizing` CLI flag, validated by PoC #39.
   The remaining primitives (`add_topo_bound_constraint`,
   `add_subtidal_flow_limiter`) are still candidates for a
   follow-up batch.

## 7. Glossary of cross-references

| Reference | Where |
| --- | --- |
| Quality benchmarks | `docs/architecture.md` §2; `notebooks/16_*.py`, `notebooks/19_*.py` |
| Gmsh over-connected investigation | PoC #25 (`notebooks/25_*.py`), #26 (`notebooks/26_*.py`) |
| Thin-channel / under-resolved detector | PoC #28 (`notebooks/28_*.py`); productionised in `fvcom_mesh_tools.diagnostics` |
| Triangle engine probe | PoC #30 (`notebooks/30_triangle_engine_poc.py`) |
| Determinism gap | `docs/python_pipeline_gap_analysis.md` §2.7 |
| License policy | `THIRD_PARTY_NOTICES.md`; oceanmesh GPL → subprocess; ocsmesh CC0; gmsh GPL runtime only |
