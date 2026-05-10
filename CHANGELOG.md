# Changelog

All notable changes to `fvcom-mesh-tools` are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The
project is **pre-alpha**; breaking changes can land on `main` between
PoC iterations. Once a versioned release tag exists, breaking changes
will only ship with a major bump (Semantic Versioning).

## Unreleased

### Highlights

A complete `clean → measure → loop` toolchain for FVCOM mesh quality
landed under this release. End-to-end:

```
fmesh-subset-dem → fmesh-buildmesh → fmesh-mesh-pipeline
                                             ↑
                                fmesh-mesh-check / -clean / -quality
```

Build, repair (7 phases A-G), unified metrics with threshold gate,
progressive 3-rung clean / quality / repeat loop with `--best-rung`
selection. ocsmesh's build engine is deprecated (PoC #30 / #36
together close the case); ocsmesh remains a library dependency for
`ops.combine_mesh` and `utils.cleanup_skewed_el` (no gmsh).

### Deprecated

- **`fmesh-buildmesh --engine ocsmesh`**. Selecting it now emits a
  `DeprecationWarning` plus a stderr notice. Quality (PoC #16:
  alpha 0.85, max valence 26) is far below the oceanmesh path
  (PoC #19: 0.96, max valence 9). PoC #30 ruled out a Triangle
  replacement (`NotImplementedError("Varying sizing is not supported
  for Triangle engine!")`); PoC #36 showed oceanmesh's `--om-max-iter
  10` matches the draft niche at 7 min wall with alpha 0.943
  (10× faster than default, still well above ocsmesh quality). The
  build path stays one release for migration. **Library** use of
  ocsmesh (`ops.combine_mesh`, `utils`, `Raster`) is unaffected.
  See `docs/engine_complementarity.md` for the full inventory.

### Added — CLIs

- **`fmesh-mesh-quality`** — unified quality metrics (alpha mean /
  p05 / p50, min-angle p05 / p50, frac<20°, max valence,
  n_overconnected, n_flipped, n_components, n_disjoint_elems) for
  one or more fort.14 inputs. Two inputs print a `delta` column;
  threshold flags (`--min-alpha`, `--max-frac-lt-20deg`,
  `--max-valence`, `--max-overconnected`, `--max-flipped`,
  `--max-disjoint-elems`) turn it into a CI gate evaluated against
  the LAST mesh (exit 1 on failure). Backed by
  `fvcom_mesh_tools.quality` with `compute_metrics` /
  `check_thresholds` / `format_comparison_table`.
- **`fmesh-mesh-pipeline`** — progressive `clean → quality → repeat`
  loop. Three cumulative rungs: rung 0 (A+B+C), rung 1 (+D+F+G),
  rung 2 (+E). Default early-stops at the first rung that satisfies
  the supplied threshold gate; `--best-rung` runs every rung up to
  `--max-iters` and writes the gate-passing rung with the highest
  `alpha_mean` (ties broken in favour of the lighter repair).
  Per-rung JSON history with `selection_reason`. PoC #33 validates
  end-to-end on the PoC #19 raw Tokyo Bay mesh under the
  FVCOM-friendly preset (`--min-alpha 0.95 --max-frac-lt-20deg
  0.005 --max-valence 8 --max-flipped 0 --max-disjoint-elems 0`):
  rung 1 (`+D+F+G`) passes, rung 2 not needed.

### Added — `fmesh-mesh-clean` phases

- **Phase A** `keep_components` — drop disjoint dual-graph
  components; default keeps only the largest.
- **Phase B** `trim_dead_ends` — iterative degree-1 trim.
- **Phase C** `repair_thin_chains` — widen 1-cell channels by
  centroid insertion (default), or delete the chain.
- **Phase D** `repair_overconnected_nodes` (off by default) —
  greedy Lawson edge swap that drives valence ≤ `--max-nbr-elem`.
  Graduated from PoC #27.
- **Phase E** `repair_under_resolved_channels` (off by default) —
  three modes: `widen` / `delete` / `medial`. The new
  `--under-resolved-min-channel-elements N` filter (default 1 = no
  filter) drops flagged elements whose face-face-connected
  component is smaller than N — see PoC #35 / #37 motivation
  below.
    * `widen` (centroid insert) lifts h_local by ~0.577× without
      changing the geometric channel width, so the post-widen w/h
      ratio is ~1.73× the original — borderline-flagged elements
      cross the threshold but very narrow channels stay flagged.
      PoC #29 validated 4.6 % reduction on PoC #19; the upper bound
      is characterised by PoC #35.
    * `medial` (Stage 2) replaces each face-face-connected channel
      of >= `min_channel_elements` flagged members with a Delaunay
      triangulation of (rim polygon ∪ centroid-spine sampled at
      `h_local_median` spacing). Skips components whose rim is
      branching or pathologically non-convex, leaving their original
      triangulation untouched. PoC #38 on the cleaned PoC #19 mesh
      with `min_channel_elements=10`: 40/51 components replaced,
      NP +206 / NE +274, alpha 0.9576 → 0.9551, frac<20° 0.17 % →
      0.33 % — 5-8× fewer new nodes / elements and 8× less alpha
      damage than `widen` at the same filter, with only modest
      growth in frac<20°. Public alias
      `fvcom_mesh_tools.mesh_clean.repair_under_resolved_channels`
      with `mode='medial'`.
- **Phase F** `repair_skewed_elements` (off by default) — wraps
  `ocsmesh.utils.cleanup_skewed_el` (gmsh-free). Deletes triangles
  whose interior angles fall outside
  `[--repair-skewed-min-angle-deg, --repair-skewed-max-angle-deg]`
  (default `[1°, 175°]`). PoC #31 sweep on the cleaned mesh: 0/3/9
  removed at default / conservative / aggressive thresholds —
  near-zero impact on already-clean output, real leverage on raw
  ocsmesh meshes.
- **Phase G** `smooth_mesh_laplacian` (off by default) — wraps
  `oceanmesh.laplacian2`. Connectivity, depths, and boundary lists
  preserved. Includes `repair_flipped=True` (default) safety net
  that reverts inverted triangles produced by the smoother — caught
  by `fmesh-mesh-quality --max-flipped 0` over PoC #19's
  raw / cleaned / Phase-G output. The same `repair_flipped_elements`
  helper now wraps the build-time `om.laplacian2` call in
  `mesh_engine/oceanmesh.py` (fixes the regression PoC #34
  surfaced). PoC #32 numbers on PoC #19's cleaned mesh:
  alpha 0.9576 → 0.9590, frac<20° drops 29 % relative.

### Added — `fmesh-buildmesh` (oceanmesh engine)

- **`--om-wavelength-sizing`** (off by default) — adds
  `oceanmesh.wavelength_sizing_function` (`dx ∝ T·√(g·h)/wl`) to
  the size composition alongside `feature_sizing_function` and
  `bathymetric_gradient_sizing_function`. The three are merged via
  `om.compute_minimum`. Tunables: `--om-wavelength-period` (default
  44712 s ≈ M2) and `--om-wavelength-grid-spacing` (default 100,
  implies dt ≈ T/wl ≈ 7.5 min). PoC #34 on Tokyo Bay: shoaling
  cells (≤ 5 m) refined +1.8 %, alpha +0.0012, frac<20° -22 %
  relative; **min CFL-feasible dt only +2 %** because Tokyo Bay's
  worst-case dt is set by coastline `feature_sizing`, not depth.
  Off-by-default is the right posture; turn on for basins with
  shoaling regions away from coastline detail.
- **Build-time `om.laplacian2` flip-rollback** — same safety net
  used by Phase G now wraps the build cleanup chain. Eliminates
  the 1 inverted triangle PoC #34 surfaced. Public alias
  `fvcom_mesh_tools.mesh_clean.repair_flipped_elements`.
- **Initial oceanmesh adapter** (`mesh_engine.oceanmesh`) — the
  default engine; pure-Python DistMesh + post-processing chain.
  Validated by PoCs #18-#22 (Tokyo Bay alpha 0.96, frac<20° 0.10 %
  vs ocsmesh+gmsh's 0.85 / 1.13 %).

### Added — diagnostics (`fmesh-mesh-check`)

- **Seven detectors** in `fvcom_mesh_tools.diagnostics`:
  `disjoint_components_flag`, `dead_end_elements_flag`,
  `thin_elements_flag`, `thin_chain_elements_flag`,
  `overconnected_nodes_flag`, `unreachable_elements_flag`,
  `under_resolved_channels_flag`. The 7th (medial-axis channel
  width) graduated from PoC #28: per-polyline cKDTrees + arc-
  separation filter + opposite-bank direction filter
  (`cos < --channel-opposite-bank-cos-max`, default −0.8).
- **`--min-channel-elements N`** filter on detector 6 (default 1 =
  no filter). Drops flagged elements whose face-face-connected
  component has fewer than N members. PoC #35 found that on real
  meshes the 3,178 default-flag elements split into 1,010
  components with mean ~3 elements / channel — mostly small
  isolated clusters, not the long ribbon-like inlets Phase E
  targets. Plumbed through `run_diagnostics`, `fmesh-mesh-check`,
  `repair_under_resolved_channels`, `clean_mesh`,
  `fmesh-mesh-clean`, and `fmesh-mesh-pipeline`.
- **`analyze_under_resolved_channels`** in `mesh_clean` — Stage 1
  measurement for the deferred "true medial-axis Phase E" project.
  Splits flagged elements into channels, reports per-channel
  `n_elements`, `h_local_median_m`, `long_axis_m`, and the
  centroid-widen vs medial-axis-to-N-cells new-node estimates.
  No re-meshing.

### Added — documentation

- `docs/engine_complementarity.md` — capability map
  oceanmesh ↔ ocsmesh; recommended division of labour ("oceanmesh
  for build, ocsmesh as a library").
- `docs/detector_repair_matrix.md` — single lookup table mapping
  each detector to the phase that fixes it, the metric that
  measures it, and the pipeline rung that turns it on.
  Phase-ordering rationale, side-effect summary, threshold-gate
  heuristics ("FVCOM-friendly preset"), recommended workflow,
  and a "where to add a new detector / phase" appendix.

### Added — infrastructure

- `.pre-commit-config.yaml` — ruff + standard hygiene hooks. CI
  installs pre-commit and runs `pre-commit run --all-files`.
- GitHub Actions CLI smoke step — `--help` invocation on every
  console script (`fmesh-buildmesh`, `fmesh-perpfix`,
  `fmesh-subset-dem`, `fmesh-mesh-combine`, `fmesh-mesh-check`,
  `fmesh-mesh-clean`, `fmesh-mesh-quality`, `fmesh-mesh-pipeline`).
- `MESH_PNG_DPI = 600` shared default in
  `fvcom_mesh_tools.plotting`.
- `dem/` subpackage isolating rasterio / netCDF4 / pyproj behind
  the `[dem]` extra (`dem.bbox.read`, `dem.subset.to_geotiff`,
  `dem.interp.at_points`).
- `pyproject.toml` extras layered by concern: `[io-vector]`,
  `[dem]`, `[oceanmesh]`, `[ocsmesh]`, `[viz]`, `[all]`.
- `mesh_engine.ocsmesh` adapter — keeps the deprecated dispatch
  path working until removal.

### PoC notes (research findings, deferred decisions)

These are the empirical results that drove the choices above.
Each links to a notebook in `notebooks/`.

- **PoC #23** — `fmesh-mesh-combine --strategy overlap` validated on
  real Tokyo Bay data (4,224-node coarse + 6,008-node inner →
  8,227 NP / 13,923 NE / alpha 0.954). `neighbor` strategy still
  pending end-to-end exercise.
- **PoC #25 / #26** — gmsh's over-connected anomaly on Tokyo Bay:
  ~380 of the 440 over-connected nodes come from gmsh itself
  before any post-processing; turning off `--refine-min-angle`
  drops 440 → 313 but cannot close the gap with the oceanmesh
  engine.
- **PoC #27** — Phase D feasibility: with the FVCOM-safe 20°
  min-angle floor every flip is rejected on real meshes, so Phase
  D defaults to floor=0° (only inversion forbidden).
- **PoC #28** — first medial-axis channel-width detector
  (cross-polyline only). Caught 1,009 elements the 1-cell
  thin-chain detector misses, but missed 691 same-polyline narrow
  inlets — productionised version added the same-polyline + cosine
  filter.
- **PoC #30** — ocsmesh's Triangle backend rejects raster-driven
  varying sizing. Drove the `--engine ocsmesh` deprecation.
- **PoC #34** — `--om-wavelength-sizing` A/B on Tokyo Bay. Refines
  shoaling cells +1.8 % and slightly improves quality, but min CFL
  dt only +2 % on this basin (coastline-pinned). Stays off by
  default.
- **PoC #35** — Stage 1 of "true medial-axis Phase E". Cleaned
  PoC #19 mesh: 3,178 flagged → 1,010 components, mean ~3 elements
  / channel; medial-axis estimate +48 % nodes vs centroid widen.
  Headline at the time: Stage 2 deferred because most channels
  looked like small isolated clusters where centroid widen is the
  right fix. PoC #37 below revisits this with the new filter and
  inverts the conclusion.
- **PoC #37** — `--min-channel-elements` sweep applied to the same
  cleaned PoC #19 mesh as Stage 1:

      min_n  flagged  comps  mean_la/h  centroid_n  medial_n  medial/centroid
      -----  -------  -----  ---------  ----------  --------  ----------------
          1    3,178  1,010      1.82       3,178     4,714        1.48
          3    2,239    301      3.06       2,239     2,112        0.94
          5    1,756    157      4.01       1,756     1,412        0.80
         10    1,070     51      6.51       1,070       710        0.66
         20      597     14     11.14         597       322        0.54

  **Stage 2 GO.** Once the small-cluster noise is filtered out
  (`min_n >= 3`) the medial-axis estimate becomes cheaper than
  centroid widen *and* the surviving channels are genuinely
  ribbon-like (mean ``long_axis_m / h_local_median`` rising from
  3.06 to 11.14). Production sweet spot: `--min-channel-elements
  10` (51 components / 1,070 elements). PoC #35's "Stage 2 not
  justified" headline came from the unfiltered noise dominating
  the average — the actual long ribbon-like channels (where
  medial-axis insertion is also topologically more correct than
  centroid widen, since centroid only lifts ``w/h`` to ~1.73× and
  cannot guarantee ``target_cells_across`` cells) are real and
  worth re-meshing.
- **PoC #38** — Stage 2 implementation validated end-to-end on the
  same cleaned PoC #19 input at `min_channel_elements=10`:

      mode      ΔNP    ΔNE   alpha     frac<20°  comps replaced
      ------    ----   ----  --------  --------  -----------------
      widen    +1070  +2140  0.9576→0.9332  0.17→1.42 %    51 / 51
      medial    +206   +274  0.9576→0.9551  0.17→0.33 %    40 / 51

  Stage 2 (``mode='medial'``) uses 5-8× fewer new nodes and
  elements, damages alpha 8× less, and limits frac<20° growth to
  +0.16 pp instead of +1.25 pp. Conservative skip rate on the 11
  rejected components: 8 non-convex rim, 3 branching rim — those
  patches keep their original triangulation as a safe fallback.
  Wall-clock identical (~2.6 s) on the 47 k-element mesh.
- **PoC #36** — `--om-max-iter` sweep on Tokyo Bay (50 → 25 → 10 → 5):

      iters   wall    alpha   frac<20°   max_v   n_overconn
      ------  ------  ------  ---------  ------  ----------
      50      26.0 m  0.9593  0.116 %    9       2
      25      14.1 m  0.9545  0.082 %    9       4
      10       6.8 m  0.9430  0.159 %   10      40
       5       4.6 m  0.9290  0.267 %   10      97

  At iters=10 oceanmesh produces alpha 0.943 — well above
  ocsmesh+gmsh's 0.847 — in 7 min. The draft niche `--engine
  ocsmesh` filled is now better served by `--om-max-iter 10` (or
  25 for "fast production"); ocsmesh's only remaining advantage
  was the ~40 s wall-clock, but the quality gap (0.943 vs 0.847,
  max_v 10 vs 26, n_overconn 40 vs 440) is so large that few real
  workflows would prefer the older path. The deprecation case is
  closed.

### Earlier groundwork (pre-CLI)

- `fvcom_mesh_tools.mesh_clean` module (initial 3-phase A+B+C
  before D/E/F/G landed) — graduated from PoC #24-#26.
- `fvcom_mesh_tools.diagnostics` module (initial 6 detectors
  before detector 7 landed).
- `fmesh-mesh-clean`, `fmesh-mesh-check` CLIs.

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
