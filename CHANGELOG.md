# Changelog

All notable changes to `fvcom-mesh-tools` are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The
project is **pre-alpha**; breaking changes can land on `main` between
PoC iterations. Once a versioned release tag exists, breaking changes
will only ship with a major bump (Semantic Versioning).

## Unreleased

### Deprecated

- **`fmesh-buildmesh --engine ocsmesh`** is now deprecated. The CLI
  emits a `DeprecationWarning` plus a visible stderr notice when
  the option is selected, and the `--help` text and CLI table label
  it accordingly. The code path is retained for one release for
  migration; production callers should switch to
  `--engine oceanmesh` (the default) immediately. Rationale is in
  `docs/engine_complementarity.md`: ocsmesh + gmsh produces
  alpha~0.85 / max valence 26 (vs 0.96 / 9 for oceanmesh) and
  ocsmesh's Triangle backend cannot consume a varying-size raster
  Hfun (PoC #30), so gmsh cannot be cheaply replaced. **Library
  use** of ocsmesh (`ops.combine_mesh` for
  `fmesh-mesh-combine --strategy {overlap,neighbor}`, `utils`,
  `Raster`) is unaffected — only the buildmesh engine path is
  going away.

### Added

- Phase G now repairs flipped triangles introduced by Laplacian
  smoothing. New helper
  `fvcom_mesh_tools.mesh_clean._repair_flipped_elements` detects
  any negative-signed-area triangle in the smoother's output and
  iteratively reverts the three nodes of every offending triangle
  to their pre-smoothing positions; the loop stops when no flips
  remain or `max_passes` (default 5) is reached. A full-rollback
  safety net guarantees the output is flip-free even when the
  iterative repair does not converge.
  `smooth_mesh_laplacian` now exposes `repair_flipped` (default
  True) and `max_repair_passes` parameters; the returned `info`
  dict gains `n_flipped_post_smooth`, `n_flipped_after_repair`
  (always 0 with `repair_flipped=True`), `n_nodes_rolled_back`,
  `n_rollback_passes`, and `full_rollback`. The `fmesh-mesh-clean`
  CLI gains `--smooth-no-repair-flipped` and
  `--smooth-max-repair-passes`. Discovered by
  `fmesh-mesh-quality` running over PoC #19's raw / cleaned /
  Phase-G output: the threshold gate `--max-flipped 0` flagged
  2 inverted triangles produced by `oceanmesh.laplacian2` on the
  cleaned mesh — a regression that the smoother's edge-length
  convergence metric does not catch on its own.
- `fmesh-mesh-quality` CLI + `fvcom_mesh_tools.quality` module: a
  unified quality-metrics dump that consolidates the per-mesh
  numbers PoCs computed ad-hoc (alpha mean / p05 / p50, min-angle
  p05 / p50, frac<20°, max valence, n_overconnected, n_flipped,
  n_components, n_disjoint_elems) into one place. `compute_metrics`
  takes a `Fort14Mesh` and returns a flat JSON-friendly dict with
  the keys listed in `quality.METRIC_KEYS`. The CLI accepts one or
  more fort.14 inputs: a single mesh prints metrics, two prints a
  side-by-side table with a `delta` column, three or more prints
  the matrix. Threshold flags (`--min-alpha`, `--max-frac-lt-20deg`,
  `--max-valence`, `--max-overconnected`, `--max-flipped`,
  `--max-disjoint-elems`) are evaluated against the LAST mesh and
  turn the command into a CI gate (exit 1 on any failure). All
  outputs land in a JSON summary alongside the last input. Validates
  the gate end-to-end: a smoke run against the PoC #19 raw / cleaned
  / Phase-G output flagged the 2 inverted triangles Phase G
  introduced — exactly the kind of regression the gate is for.
- `fmesh-mesh-clean` Phase G: Laplacian smoothing of interior nodes.
  New function `fvcom_mesh_tools.mesh_clean.smooth_mesh_laplacian`
  wraps `oceanmesh.laplacian2`, which derives the boundary set from
  the mesh topology and pins it automatically — connectivity, depth
  array, and open / land boundary lists are all preserved across
  the smoothing pass. Off by default. The CLI gains
  `--smooth-laplacian` plus `--smooth-laplacian-iters` (default 20)
  and `--smooth-laplacian-tol` (default 0.01) — both matching the
  oceanmesh defaults. Importing oceanmesh propagates GPL-3.0 into
  the redistributed combined work (already documented in
  `THIRD_PARTY_NOTICES.md`); callers who need a GPL-free path
  should leave Phase G off. PoC #32
  (`notebooks/32_phase_g_smooth_poc.py`) sweeps three iteration /
  tolerance presets on the PoC #19 cleaned Tokyo-Bay mesh:
  alpha-mean `0.9576 → 0.9590` (+0.0014), min-angle p05
  `39.90° → 40.27°` (+0.37°), bad-triangle fraction (frac<20°)
  `0.169 % → 0.120 %` (≈ 29 % relative reduction, 80 → 57 of
  47,409). Convergence is fast on cleaned input — the gentle
  preset (`max_iter=5`) yields the same numbers as the default
  (`max_iter=20`). Topology invariants (NP / NE / boundary
  counts) preserved across all presets.
- `fmesh-mesh-clean` Phase F: angle-based skewed-element removal. New
  function `fvcom_mesh_tools.mesh_clean.repair_skewed_elements` wraps
  `ocsmesh.utils.cleanup_skewed_el` (gmsh-free; ocsmesh used as a
  library only — see `docs/engine_complementarity.md` §3.2.3) and
  deletes triangles whose minimum interior angle is below
  `--repair-skewed-min-angle-deg` (default 1°) or whose maximum is
  at or above `--repair-skewed-max-angle-deg` (default 175°). The
  CLI gains `--repair-skewed-elements` (off by default) plus the
  two threshold flags. Boundaries are re-derived via DEM-bbox
  proximity after deletion (skipped when the run is a no-op).
  Validated by PoC #31 (`notebooks/31_phase_f_skewed_clean_poc.py`)
  on the PoC #19 cleaned Tokyo-Bay mesh across three threshold
  presets: at ocsmesh defaults `[1°, 175°]` the cleaned mesh has zero
  flagged elements (the prior 4-phase pipeline already left no
  degenerate slivers); at `[5°, 170°]` 3 of 47,409 elements (0.006 %)
  are removed; at `[10°, 160°]` 9 are removed (0.019 %). Phase F's
  real leverage is therefore on raw / unclean meshes where slivers
  survive — particularly the OCSMesh+gmsh path (alpha 0.85,
  frac<20°=1.13 %) — rather than on already-cleaned oceanmesh output.
- `docs/engine_complementarity.md` consolidates the empirical and
  source-level investigation of `oceanmesh` and `ocsmesh`: which
  capabilities each library has, which are exclusive (a long list on
  both sides), and which are gmsh-dependent. Headline findings:
  ocsmesh's Triangle engine cannot consume a varying-size `Hfun`
  (PoC #30 `NotImplementedError("Varying sizing is not supported for
  Triangle engine!")`), so it is not a drop-in replacement for gmsh
  in the build path. Several ocsmesh capabilities have no oceanmesh
  equivalent and are independently valuable to keep:
  `add_courant_num_constraint` and the other `Hfun.add_*` sizing
  primitives, `ops.combine_mesh.merge_overlapping_meshes` and
  `merge_neighboring_meshes` (Triangle-based, gmsh-free),
  `utils.cleanup_skewed_el` and `repartition_features`,
  `utils.interpolate_*` for mesh-to-mesh field transfer, and
  `Raster.{clip,fill_nodata,gaussian_filter,get_channels}` for DEM
  preprocessing. The recommended division of labour is "oceanmesh
  for build, ocsmesh as a library", with `--engine ocsmesh`
  deprecation planned (decision pending) and ocsmesh utility wrappers
  evaluated for Phase F sliver clean and Phase G Laplacian smoothing.
- PoC #30 (`notebooks/30_triangle_engine_poc.py`) tries swapping
  ocsmesh's `MeshDriver(engine_name="gmsh")` for `engine_name="triangle"`
  on the Tokyo Bay PoC #16 inputs. ocsmesh's Triangle wrapper raises
  `NotImplementedError` whenever a raster-driven Hfun is supplied —
  it only supports constant size — so Triangle is unusable for our
  build configuration without a from-scratch wrapper around Shewchuk
  Triangle. Documented in `docs/engine_complementarity.md` §2 as the
  reason gmsh cannot be cheaply removed.
- `.pre-commit-config.yaml` for contributor-side git hooks: `ruff
  check --fix`, trailing-whitespace + EOF-newline + mixed-line-ending
  fixers, YAML/TOML syntax check, merge-conflict marker check,
  large-file guard. CI installs `pre-commit` and runs
  `pre-commit run --all-files`, so local hooks and CI agree.
- CLI smoke step in GitHub Actions runs `--help` on every console
  script (`fmesh-buildmesh`, `fmesh-perpfix`, `fmesh-subset-dem`,
  `fmesh-mesh-combine`, `fmesh-mesh-check`, `fmesh-mesh-clean`) so
  any regression in the `pyproject.toml` entry-point wiring or
  argparse construction trips CI.

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
- `fmesh-mesh-clean` Phase D: over-connected node repair via
  valence-balancing edge swaps (graduated from PoC #27). The greedy
  Lawson-style flip is exposed as
  `fvcom_mesh_tools.algorithms.swap_edges_for_valence` and as
  `fvcom_mesh_tools.mesh_clean.repair_overconnected_nodes`. The CLI
  gains `--repair-overconnected-iters` (default 0 = OFF),
  `--max-nbr-elem` (default 8 = FVCOM legacy cap), and
  `--overconn-min-angle-floor` (default 0° — only triangle inversion
  forbidden, the value PoC #27 found practical on real meshes).
  Validated on PoC #19's cleaned Tokyo Bay mesh: enabling Phase D
  drives `n_overconnected: 3 → 0`, `max_valence: 9 → 8` after 20
  swaps in 2 iterations, with alpha mean 0.9577 → 0.9576, frac<20° =
  0.16 → 0.17 % — essentially zero quality cost. Severe gmsh-fan
  cases (PoC #16, max valence 26) are only partially fixable by edge
  swap alone; mitigation there is engine choice.
- `fmesh-mesh-clean` Phase E: widen or delete medial-axis-detected
  under-resolved channel elements (detector 6). Reuses the existing
  centroid-insertion mechanism from Phase C-widen — each flagged
  triangle becomes 3 sub-triangles fanning from a new interior
  centroid. Exposed as
  `fvcom_mesh_tools.mesh_clean.repair_under_resolved_channels`. The
  CLI gains `--under-resolved-mode {widen,delete,none}` (default
  `none`), `--under-resolved-min-w-h` (default 3.0), and the matching
  detector-6 parameters (`--under-resolved-sample-ds-m`,
  `--under-resolved-arc-separation-factor`,
  `--under-resolved-opposite-bank-cos-max`). Phase E is **off by
  default** because detector 6 typically flags thousands of elements
  on real meshes; enable deliberately when widening is the desired
  remediation. PoC #29 validates the widen path end-to-end on the
  PoC #19 cleaned Tokyo-Bay mesh: topology growth checks pass
  (NP +3,178, NE +6,356) and boundaries are preserved, but the
  detector-6 reduction is modest (3,178 → 3,032, 4.6 %). Reason:
  centroid insertion shrinks h_local by ~0.577× while the geometric
  channel width is unchanged, so w/h ≈ 1.73× the original ratio —
  only elements with originally-borderline ratios cross the
  threshold. Phase E should therefore be read as "lift local
  resolution one step", not "guarantee 3 cells across every narrow
  channel". The latter requires inserting nodes along the channel
  medial axis, a deeper remeshing operation outside `clean_mesh`'s
  scope.
- 7th detector `under_resolved_channels_flag` graduates from PoC #28
  into `fvcom_mesh_tools.diagnostics`. The metric is the smaller of
  two channel-width candidates divided by the median element edge
  length:
    1. **cross-polyline**: distance from the centroid to the two
       nearest distinct boundary polylines, summed (catches the
       "channel between mainland and an island" case);
    2. **same-polyline narrow inlet**: distance to the nearest
       sample on a polyline plus the distance to the nearest sample
       on that same polyline whose along-polyline arc separation is
       large (catches an inlet whose two banks lie on a single
       continuous coastline). The same-polyline candidate is
       accepted only when the two vectors (centroid → nearest
       sample, centroid → far-arc sample) point in roughly opposite
       directions (cos angle < `--channel-opposite-bank-cos-max`,
       default −0.8 = angle > 143°), which rejects coastal-corner
       false positives where the polyline wraps around a peninsula
       tip.
  Built per-polyline `cKDTree`s rather than a single combined tree,
  so cross-polyline distance is exact even on densely-sampled
  coasts. The CLI gains `--min-w-h` (default 3.0),
  `--channel-sample-ds-m` (default 50 m),
  `--channel-arc-separation-factor` (default 4.0), and
  `--channel-opposite-bank-cos-max` (default −0.8). Validated on the
  PoC #19 cleaned Tokyo Bay mesh: 3,178 elements (6.7 % of NE) are
  flagged, concentrated at the northern Tokyo Bay river mouths and
  along narrow coastal jetties / breakwaters. PoC #28's prototype
  underflagged (618) by missing same-polyline inlets; the
  productionised version's same-polyline + direction filter catches
  them while keeping the median ratio above 26 for the bay
  interior.
- PoC #28 (`notebooks/28_channel_width_poc.py`) prototypes a
  medial-axis-style channel-width / h ratio detector for
  under-resolved channels (2- to 3-cell wide) that the existing
  `thin_chain_elements_flag` (1-cell only) misses. Per element, the
  metric is `(d(centroid, polyline_A) + d(centroid, polyline_B)) /
  median_edge_length`, where `polyline_A`, `polyline_B` are the two
  closest *distinct* boundary polylines. Flag when the ratio is below
  3.

  Validated on the PoC #19 mesh (uncleaned and after Phase A+B+C+D
  cleaning):

    * On the raw mesh, the detector catches 1,145 elements at
      threshold 3, of which 1,009 are not flagged by `thin_chain` —
      these include the genuinely-under-resolved 2-cell channels
      around the Edogawa / Arakawa / Sumida river mouths in the
      northern bay, exactly the case the user originally raised.
    * On the cleaned mesh, 618 elements remain flagged — Phase C
      widened the 1-cell chains but the channels are still only
      ~2 cells wide.
    * Limitation: 691 thin-chain elements are *not* caught because
      the metric requires two distinct polylines. Where a continuous
      coastline carves a narrow inlet against itself (both banks on
      the same polyline), the second-nearest polyline is far away
      and the ratio inflates. Productionising this detector needs
      (i) splitting polylines at concave corners, (ii) Voronoi-based
      true medial axis, or (iii) a "K-nearest with along-polyline
      separation" filter.

  The PoC stays as a research notebook; productionising into a
  diagnostics-module flag is left for follow-up.
- PoC #27 (`notebooks/27_overconn_repair_poc.py`) explores
  edge-swap-based repair of over-connected nodes. A greedy
  Lawson-style flip scored by reduction of per-edge "valence excess"
  is shown to (i) eliminate the 3 over-connected nodes on the
  cleaned PoC #19 mesh (max v=9 → 8, +0.01% bad triangles), (ii)
  reduce 440 → 365 over-connected nodes on the PoC #16 OCSMesh+rivers
  mesh (max v=26 → 23, +0.25% bad triangles), but only when the
  min-angle floor is relaxed to 0°. The FVCOM-safe 20° floor rejects
  every flip on the real meshes; valence fixing is therefore
  fundamentally at odds with strict quality preservation when the
  over-connected node sits in a fan-like local topology. The
  algorithm is staged for productionising into a Phase D of
  `fmesh-mesh-clean`; the floor / threshold defaults are still under
  review.
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
