# fvcom-mesh-tools

> ⚠️ **Pre-alpha** — under active development. APIs and CLIs are unstable.

Python toolkit for FVCOM unstructured mesh generation, repair, and visualization.

`fvcom-mesh-tools` provides a unified Python interface for building
high-quality FVCOM-ready unstructured meshes (`fort.14`), with a focus on:

- **Open-boundary edge orthogonalization** — enforce edges perpendicular to open boundaries
- **Mesh defect detection** — `fmesh-mesh-check` flags disjoint wet pools, dead-ends, 1-cell-wide channels (`thin_chain`), under-resolved 2-3 cell channels (medial-axis-style `w/h` ratio), over-connected nodes, and open-boundary-unreachable elements (no repair, JSON + map output)
- **Safe automated repair** — `fmesh-mesh-clean` prunes disjoint pools, trims dead-end "spits", and widens 1-cell channels to 2 cells via centroid insertion
- **Mesh quality inspection and visualization** — element quality, boundary classification, fort.14 plots

The package wraps several mature mesh tools behind a common backend interface
rather than reimplementing meshing algorithms from scratch.

## Backend strategy

| Role | Backend | License | How used |
|------|---------|---------|----------|
| Mesh generation (default) | [oceanmesh](https://github.com/CHLNDDEV/oceanmesh) (DistMesh) | GPL-3.0-or-later | imported |
| Mesh generation (alt / draft, **deprecated**) | [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) + [gmsh](https://gmsh.info/) | CC0-1.0 + GPL-2.0+ runtime | imported (OCSMesh); gmsh called via OCSMesh |
| Multi-mesh stitching | [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) (`ops.combine_mesh`) | CC0-1.0 | imported |
| Auxiliary mesh utilities (legacy) | [MeshKernelPy](https://github.com/Deltares/MeshKernelPy), [stompy](https://github.com/rustychris/stompy) | MIT | optional, imported |

`oceanmesh` is GPL-3.0; the combined work, when redistributed together
with this toolkit, must respect GPL-3.0. The `--engine ocsmesh` path is
GPL-only-runtime (gmsh is invoked as an external tool through OCSMesh,
not linked). See `THIRD_PARTY_NOTICES.md` and `docs/architecture.md`
for the full rationale.

## Installation

Python ≥3.12. Two install modes, depending on whether you want the
compiled scientific stack from conda-forge or a pip-only setup.

### Conda (recommended for full functionality)

```bash
mamba env create -f environment.yml
mamba activate fvcom-mesh
pip install --no-deps oceanmesh        # GPL-3.0-or-later
pip install --no-deps -e .             # this package, editable
```

Two pre-built env scripts under `notebooks/` reproduce the exact GENKAI
setup used to validate PoCs #18-#23:

| Script | Purpose |
|--------|---------|
| `notebooks/18_setup_oceanmesh_env.pjsub` | initial env (`oceanmesh-bench`) creation |
| `notebooks/18_setup_gdal_netcdf.pjsub`  | adds `libgdal-netcdf` so rasterio reads CF NetCDF |
| `notebooks/19_setup_env_extend.pjsub`   | editable install of this package |
| `notebooks/21_setup_ocsmesh_in_omsh.pjsub` | adds OCSMesh on top for `fmesh-mesh-combine` |

If you only need the OCSMesh path, the lighter `py312test` env (already
on GENKAI) works - oceanmesh is then unavailable and the default engine
must be flipped via `--engine ocsmesh`.

### Pip (modular extras)

Dependencies are layered by concern so pure-fort.14 callers do not pull
the compiled raster stack:

| Install | Pulls | Enables |
|---------|-------|---------|
| `pip install -e .` | numpy only | `algorithms.*`, `io.fort14`, `mesh_compose.disjoint`, `fmesh-perpfix`, `fmesh-mesh-combine --strategy disjoint` |
| `pip install -e ".[io-vector]"` | + shapely / geopandas / fiona | coastline / river-point / multipolygon-area helpers |
| `pip install -e ".[dem]"` | + rasterio / netCDF4 / pyproj | `dem.subset` / `dem.interp` / `dem.bbox`, `fmesh-subset-dem` |
| `pip install -e ".[oceanmesh]"` | + oceanmesh and the above | `fmesh-buildmesh --engine oceanmesh` (default) |
| `pip install -e ".[ocsmesh]"` | + ocsmesh, gmsh, and the above | `fmesh-mesh-combine --strategy {overlap,neighbor}` (Triangle-based, gmsh-free at runtime); `fmesh-buildmesh --engine ocsmesh` (**deprecated**, slated for removal) |
| `pip install -e ".[all]"` | superset of the above plus `viz` | every CLI and helper |

Note: `oceanmesh` on PyPI lists deps that conflict with conda-forge
versions; under `[oceanmesh]` pip will pull the PyPI variant. Under
the conda workflow we install `oceanmesh --no-deps` and rely on the
conda-forge stack instead.

## Quick start

Mesh Tokyo Bay end-to-end with the default oceanmesh engine:

```bash
# 1. (Optional) Clip a global DEM to a regional bounding box
fmesh-subset-dem data/SRTM15+.nc /tmp/tb.tif \
    --bbox 139.5 35.1 140.2 35.9 \
    --src-var z

# 2. Generate the mesh (DistMesh + post-processing chain)
fmesh-buildmesh /tmp/tb.tif /tmp/tokyo_bay.14 \
    --engine oceanmesh \
    --hmin 200 --hmax 5000 \
    --coastline data/coastline/MLIT_C23/C23-06_TOKYOBAY.shp \
    --river-inflow-points data/rivers/tokyo_bay/tokyo_bay_rivers.csv \
    --om-seed 0 \
    --perpfix-iters 1
```

> **`--engine ocsmesh` is deprecated.** It is kept one release for
> migration but emits a `DeprecationWarning` and a stderr notice when
> selected. Production meshes should use `--engine oceanmesh` (the
> default). The OCSMesh + gmsh path produces alpha~0.85 / max valence
> 26 vs. 0.96 / 9 for oceanmesh, and PoC #30 confirmed ocsmesh's
> Triangle backend cannot replace gmsh under varying sizing. See
> `docs/engine_complementarity.md` for the rationale. **Library** use
> of ocsmesh (`ops.combine_mesh`, `utils`, `Raster`) is unaffected.

To stitch independently generated meshes:

```bash
fmesh-mesh-combine tokyo_bay.14 osaka_bay.14 kanto_kansai.14 \
    --strategy disjoint
```

### Diagnose and repair an existing mesh

`fmesh-mesh-check` runs seven detectors over a `fort.14` and emits a
summary, a JSON dump of every flagged element / node (with
coordinates), and a map PNG. Exit code is non-zero when anything is
flagged so it works as a CI gate:

```bash
fmesh-mesh-check tokyo_bay.14 \
    --max-nbr-elem 8 \           # FVCOM MAX_NBR_ELEM cap (set to your build)
    --min-w-h 3.0                # min cells across an under-resolved channel
```

The detectors are: (1) disjoint dual-graph components, (2) dead-end
elements, (3) thin elements, (3b) thin chains (1-cell channels), (4)
over-connected nodes, (5) open-boundary-unreachable elements, and
(6) under-resolved channels — elements whose local channel width
(sum of distances to the two nearest non-adjacent boundary samples,
detected per polyline with an along-arc separation filter and an
opposite-bank direction check) divided by the median edge length is
below `--min-w-h`. Detector 6 catches 2- and 3-cell channels that
detector 3b misses.

`fmesh-mesh-clean` repairs the safe-to-fix subset that
`fmesh-mesh-check` surfaces — disjoint pools, dead-end "spits",
1-cell-wide channels, over-connected nodes, and (optionally)
under-resolved 2-cell channels. Pass the original DEM bbox (and the
same `--open-merge-coast-gap` you used at build time) so the
re-derived boundaries match `fmesh-buildmesh`:

```bash
fmesh-mesh-clean tokyo_bay.14 tokyo_bay_clean.14 \
    --bbox 139.565 35.102 140.171 35.856 \
    --open-merge-coast-gap 50 \
    --thin-chain-mode widen        # default; alternative: 'delete' or 'none'
```

Phase A keeps only the largest dual-graph component (or, with
`--require-open-boundary` / `--min-component-elements N`, a
configurable subset). Phase B trims degree-1 elements that have no
open-boundary edge, iterating until convergence. Phase C
(`--thin-chain-mode widen`, default) inserts a centroid into every
thin-chain element so each 1-cell channel becomes 2-cell;
`--thin-chain-mode delete` removes the chain instead, and
`--thin-chain-mode none` skips Phase C. Phase D
(`--repair-overconnected-iters N`, off by default) runs a greedy
Lawson edge-swap that drives every node valence to at most
`--max-nbr-elem` (default 8 = FVCOM legacy cap); with the default
`--overconn-min-angle-floor 0` it eliminates mild over-connection
(max v=9) at near-zero quality cost on real meshes. Severe gmsh-fan
cases need engine-level fixes instead. Phase E
(`--under-resolved-mode {widen,delete,none}`, default `none`)
widens (or deletes) elements flagged by detector 6 — the
medial-axis-style channel-width metric — using the same centroid
insertion as Phase C-widen, lifting 2-cell channels to 3-cell.
Detector 6 typically flags thousands of elements on real meshes, so
enable Phase E deliberately. Phase F
(`--repair-skewed-elements`, off by default) deletes triangles
whose minimum interior angle is below
`--repair-skewed-min-angle-deg` (default 1°) or whose maximum is
at or above `--repair-skewed-max-angle-deg` (default 175°). Wraps
`ocsmesh.utils.cleanup_skewed_el` (gmsh-free; ocsmesh used as a
library only — see `docs/engine_complementarity.md` §3.2.3).
Phase G (`--smooth-laplacian`, off by default) runs Laplacian
smoothing of all interior nodes via `oceanmesh.laplacian2`.
Boundary nodes are auto-pinned by oceanmesh; connectivity, depths,
and boundary lists are preserved. Tunable with
`--smooth-laplacian-iters` (default 20) and
`--smooth-laplacian-tol` (default 0.01). The smoother converges on
edge-length stability but does not check signed area, so by
default Phase G also runs an iterative rollback that detects any
flipped triangle in oceanmesh's output and reverts the offending
nodes back to their pre-smoothing positions; pass
`--smooth-no-repair-flipped` to surface the raw output instead.

`fmesh-mesh-quality` is the unified metrics + threshold-gate
companion to `fmesh-mesh-check` and `fmesh-mesh-clean`. It computes
`alpha_mean / alpha_p05 / min_angle_p05_deg / frac_lt_20deg /
max_valence / n_overconnected / n_flipped / n_components /
n_disjoint_elems` for one or more fort.14 files, prints a
side-by-side comparison (with a `delta` column when two meshes are
passed), and turns into a CI gate when threshold flags are
supplied. Exit 1 on any failure:

```bash
# Single mesh
fmesh-mesh-quality tokyo_clean.14

# Before/after comparison with a delta column
fmesh-mesh-quality tokyo.14 tokyo_clean.14 --labels before after

# CI gate
fmesh-mesh-quality tokyo_clean.14 \
    --min-alpha 0.95 --max-frac-lt-20deg 0.005 \
    --max-valence 8 --max-flipped 0
```

`fmesh-mesh-pipeline` chains the clean / quality steps into a
progressive `clean → quality → repeat` loop. It applies three
cumulative *rungs* of `fmesh-mesh-clean` phases, evaluating quality
thresholds after each, and stops at the first rung that passes —
or exits 1 if no rung does. The rungs are:

* **rung 0** — A + B + C (the conservative `fmesh-mesh-clean`
  default: drop disjoint pools, trim dead-end spits, widen 1-cell
  channels).
* **rung 1** — rung 0 + D + F + G: Lawson edge swap to balance
  over-connected nodes, skewed-element removal, Laplacian
  smoothing (with the flipped-triangle safety net).
* **rung 2** — rung 1 + E: widen detector-6 under-resolved
  channels; the most destructive rung (~3× new elements per
  flagged element).

```bash
fmesh-mesh-pipeline tokyo_raw.14 tokyo_passing.14 \
    --bbox 139.46 34.99 140.10 35.74 \
    --open-merge-coast-gap 50 \
    --min-alpha 0.95 --max-frac-lt-20deg 0.005 \
    --max-valence 8 --max-flipped 0 --max-disjoint-elems 0
```

The pipeline writes a JSON history that records each rung's metrics
and threshold-check results so the caller can audit which rung met
the gate. By default the loop early-stops at the first rung that
passes; pass `--best-rung` to disable the early-stop, run every
rung up to `--max-iters`, and pick the gate-passing rung that
maximises `alpha_mean` (ties broken in favour of the lighter
repair). Useful when one wants the maximum quality the pipeline can
produce, not just the first acceptable mesh.

`docs/architecture.md` is the full decision tree for engine choice and
combine strategy; `docs/python_pipeline_gap_analysis.md` has the
quality / runtime numbers vs. the OceanMesh2D MATLAB reference;
`docs/detector_repair_matrix.md` maps each `fmesh-mesh-check` detector
to the `fmesh-mesh-clean` Phase that fixes it, the
`fmesh-mesh-quality` metric that measures it, and the
`fmesh-mesh-pipeline` rung that turns it on automatically.

## Development

```bash
make install   # editable install with dev deps
make test      # pytest -q
make lint      # ruff check
```

A `.pre-commit-config.yaml` is checked in for contributor convenience.
After installing dev deps, run:

```bash
pre-commit install        # one-time: register the git hook
pre-commit run --all-files  # equivalent of what CI runs
```

The hook runs the same `ruff check --fix` and standard hygiene checks
(trailing whitespace, EOF newline, no merge conflict markers, no
oversized files) that GitHub Actions enforces, so local commits stay
in sync with CI expectations.

## Proof-of-concept notebooks

End-to-end smoke tests under `notebooks/` (each ships with a matching
`.pjsub` wrapper for GENKAI). Outputs land in `outputs/`.

| # | Notebook | What it does |
|---|----------|--------------|
| 01 | `01_read_fort14_summary.py` | Parse `tb_futtsu20220311.14`, plot the mesh, dump shape stats. |
| 02 | `02_open_boundary_orthogonality.py` | Measure FVCOM open-boundary edge perpendicularity on the reference mesh. |
| 03 | `03_mk_orthogonalize_effect.py` | Run MeshKernel global orthogonalize and confirm it does **not** improve open-boundary perpendicularity. |
| 04 | `04_boundary_perp_fix.py` | Custom first-ring perpendicularity fix; ships as the `fmesh-perpfix` CLI. |
| 05 | `05_ocsmesh_minimal.py` | Minimal OCSMesh pipeline: DEM → Geom/Hfun → gmsh → fort.14. |
| 06 | `06_parity_compare.py` | Edge-length / triangle-quality / boundary parity vs the reference mesh. |
| 07 | `07_buildmesh_e2e.py`   | End-to-end `fmesh-buildmesh` validation: depth interp + boundary classification + perpfix. |
| 08 | `08_quality_pass.py`    | Damped Laplacian smoothing alone — finding: insufficient by itself for slivers. |
| 09 | `09_edge_swap.py`       | Lawson / min-angle edge swap — monotonically improves quality. |
| 10 | `10_swap_smooth_combo.py` | Alternating swap + smooth; plateau at ~17 % bad triangles, set by initial sizing. |
| 11 | `11_buildmesh_quality_pass.py` | `fmesh-buildmesh --quality-pass 6` end-to-end smoke test. |
| 12 | `12_coastline_aware.py`  | Coastline-aware sizing (`Hfun.add_feature`); cuts bad-triangle fraction to **2.7 %**. |
| 13 | `13_island_filter.py`    | Island / polygon area filters; cuts 166 land segments → 23. |
| 14 | `14_boundary_cleanup.py` | Open-segment merge collapses 3 open segments → 1 contiguous arc. |
| 15 | `15_refine_bad.py`       | Longest-edge bisection refinement; cuts bad-triangle fraction to **1.82 %**. |
| 16 | `16_river_inflow.py`     | River-mouth ibtype=21 segments via `--river-inflow-points`; matches reference for Tokyo Bay (5 rivers). |
| 17 | `17_osaka_bay_validation.py` | Second-basin sanity check: same flag set on Osaka Bay (SRTM15+ subset + GSHHS-f L1). Three independent open arcs detected, alpha 0.90, frac<20° 0.18 %. |
| 18 | `18_oceanmesh_benchmark.py` | Head-to-head OCSMesh vs `oceanmesh` (DistMesh) on identical Tokyo Bay inputs. `oceanmesh` produces alpha 0.961 / frac<20° 0.034 % (vs 0.847 / 1.13 %); 39× slower. |
| 19 | `19_oceanmesh_full_pipeline.py` | Full `fmesh-buildmesh --engine oceanmesh` end-to-end on Tokyo Bay: alpha 0.959, frac<20° 0.10 %, 5 ibtype=21 river segments, perpfix reverts 8 (vs 273). |
| 20 | `20_osaka_bay_oceanmesh.py` | Second-basin de-risk of `--engine oceanmesh`: Osaka Bay end-to-end → alpha 0.966 / frac<20° 0.08 %, 3 open arcs, 4 ibtype=21 segments, no parameter changes from Tokyo Bay. |
| 21 | `21_mesh_combine_kanto_kansai.py` | `fmesh-mesh-combine --strategy disjoint`: stitch Tokyo Bay (PoC #19) + Osaka Bay (PoC #20) into one fort.14 with all boundaries (4 open + 138 land + 9 ibtype=21 river) preserved. |
| 22 | `22_minimum_area_mult_sweep.py` | Sweep `--om-minimum-area-mult` 1.0→2000.0 on Tokyo Bay; `om.Shoreline` inner-polygon count drops 53→39→27→19→5→0, confirming the new flag governs island filtering at the source. |
| 23 | `23_mesh_combine_overlap.py` | `fmesh-mesh-combine --strategy overlap` real-data validation: stitch a coarse Tokyo Bay outer (hmin=1000 m, NP=4,224) with a fine northern-bay inner (hmin=200 m, bbox 139.78–140.0 × 35.55–35.75, NP=6,008) via `ocsmesh.ops.merge_overlapping_meshes`. Combined NP=8,227 NE=13,923, alpha 0.954, frac<20° 0.09 %, no flipped triangles; edge length p50/p95 grades from ~393/1187 m (outer) through ~175/986 m (combined) to ~66/447 m (inner). |
| 24 | `24_mesh_check_poc.py` | Initial validation of the six `fmesh-mesh-check` detectors on the existing PoC #19 / #16 / #20 outputs: 144 disjoint components and 5,496 unreachable elements on the Tokyo Bay oceanmesh mesh, 440 over-connected nodes (max valence 26) on the OCSMesh+rivers mesh. |
| 25 | `25_overconnected_investigate.py` | Characterises the 440 over-connected nodes from PoC #24: 94.5 % are interior, 0 % are river-segment nodes, the worst clusters sit on the south open boundary or in tortuous coastal features. |
| 26 | `26_ocsmesh_overconn_ablation.py` | Post-processing ablation on Tokyo Bay+rivers (OCSMesh engine fixed): turning off `--refine-min-angle` drops over-connected nodes 440 → 313 (max valence 26 → 21). gmsh itself accounts for ~380 of them; the gap with the oceanmesh engine (3 over-connected, max v=9) cannot be closed by post-processing alone. |
| 27 | `27_overconn_repair_poc.py` | Greedy edge-swap repair of over-connected nodes scored by per-edge "valence excess" reduction. With `--overconn-min-angle-floor 0` it reduces 3 over-connected nodes (max v=9) to 0 on the cleaned PoC #19 mesh at +0.01% bad triangles, and 440 → 365 (max v 26 → 23) on the OCSMesh+rivers PoC #16 mesh at +0.25% bad triangles. With the FVCOM-safe 20° floor no swap is accepted on either real mesh, confirming that valence fixing requires accepting some sliver creation. |
| 28 | `28_channel_width_poc.py` | Medial-axis-style detector for under-resolved channels via `(d_to_polyline_A + d_to_polyline_B) / median_edge_length`. On the raw PoC #19 mesh it flags 1,145 elements at threshold 3, of which 1,009 are missed by the existing 1-cell `thin_chain` detector — including the under-resolved 2-cell channels around the northern Tokyo Bay river mouths. On the A+B+C+D-cleaned mesh, 618 elements remain flagged because Phase C only widens 1-cell chains. Known limitation: same-polyline narrow inlets (both banks on a continuous coastline) are missed; productionising needs polyline splitting or true Voronoi medial-axis. |
| 29 | `29_phase_e_widen_poc.py` | End-to-end validation of `fmesh-mesh-clean` Phase E (`--under-resolved-mode widen`) on the PoC #19 cleaned Tokyo-Bay mesh. Topology checks pass (NP +3,178, NE +6,356; boundaries preserved), but detector-6 flagged count drops only modestly (3,178 → 3,032; 4.6 %). Reason: centroid insertion shrinks h_local by ~0.577× while the geometric channel width is unchanged, so post-widen w/h ≈ 1.73× the original — only borderline-flagged elements cross the threshold. Phase E is "lift local resolution one step", not "guarantee 3 cells across every narrow channel". |
| 30 | `30_triangle_engine_poc.py` | Probe ocsmesh's Triangle backend on Tokyo Bay. Result: `NotImplementedError("Varying sizing is not supported for Triangle engine!")`. ocsmesh's Triangle wrapper accepts only constant size, so it cannot replace gmsh in our build path. Drove the deprecation of `--engine ocsmesh` (see `docs/engine_complementarity.md` §2). |
| 31 | `31_phase_f_skewed_clean_poc.py` | Phase F (`--repair-skewed-elements`) sweep on the PoC #19 cleaned mesh across three threshold presets. ocsmesh defaults `[1°, 175°]` flag 0; conservative `[5°, 170°]` flag 3 of 47,409 (0.006 %); aggressive `[10°, 160°]` flag 9 (0.019 %). Phase F has near-zero impact on already-clean oceanmesh output; its leverage is on raw or OCSMesh+gmsh meshes where slivers survive. |
| 32 | `32_phase_g_smooth_poc.py` | Phase G (`--smooth-laplacian`) sweep on the PoC #19 cleaned mesh across three iter/tol presets (default 20/1e-2, gentle 5/1e-2, deep 50/1e-4). Records max/mean node displacement, alpha-mean shift, frac<20° shift, and topology-preservation invariants (NP / NE / boundary counts unchanged). |
| 33 | `33_pipeline_poc.py` | End-to-end validation of `fmesh-mesh-pipeline` on the PoC #19 raw Tokyo-Bay mesh (144 components, 5,496 disjoint elements, max valence 9). Threshold preset `--min-alpha 0.95 --max-frac-lt-20deg 0.005 --max-valence 8 --max-flipped 0 --max-disjoint-elems 0`. Rung 0 (A+B+C) drops disjoint elements but leaves valence=9 — gate FAILS. Rung 1 (+D+F+G) brings valence to 8 and alpha to 0.9588 — gate PASSES. Rung 2 not needed. |
| 34 | `34_wavelength_sizing_poc.py` | A/B comparison of `--om-wavelength-sizing` on Tokyo Bay vs. PoC #19's gradient-only baseline. The wavelength contribution does refine shoaling regions as designed (n_shallow≤5m +1.8 %) and slightly improves quality (alpha 0.9594 → 0.9606; frac<20° -22 % rel). But min CFL-feasible dt at C=0.7 only goes 1.80 → 1.84 s (+2 %): Tokyo Bay's worst-case dt is set by `feature_sizing_function` along the coast, not by the gradient or wavelength fields. Off-by-default is the right posture; turn on when the basin's smallest cells sit in shoaling regions away from coastline detail. Side finding: 1 flipped triangle introduced by oceanmesh's build-time `laplacian2` (no flip-safety net at that call site — tracked as a follow-up). |
| 35 | `35_phase_e_potential_poc.py` | Stage 1 of the "true medial-axis Phase E" project. Splits detector-6-flagged elements into face-face-adjacent connected channel components and reports new-node cost: existing centroid-widen vs medial-axis-to-3-cell-width estimate. On the PoC #19 cleaned mesh: 3,178 flagged → 1,010 components (mean ~3 elements/component); centroid 3,178 vs medial-axis 4,714 nodes (+48 %). After centroid widen: 3,032 → 807 components, +27 %. Conclusion: the channels are predominantly small isolated clusters at river-mouth corners / jetty tips, not long ribbon-like narrow inlets — Stage 2 (real CDT re-meshing) is deferred until detector 6 gets a min-component-size filter that isolates the cases where medial-axis insertion actually wins. |
| 36 | `36_om_max_iter_sweep_poc.py` | `--om-max-iter` sweep on Tokyo Bay (50 / 25 / 10 / 5) under PoC #19 settings. iters=50 → alpha 0.9593 / frac<20° 0.116 % / max_v 9 in 26 min; iters=25 → 0.9545 / 0.082 % / 9 in 14 min; iters=10 → 0.9430 / 0.159 % / 10 in 7 min (still well above ocsmesh+gmsh's 0.847 / 1.13 % / 26 in ~40 s); iters=5 → 0.9290 / 0.267 % / 10 in 4.6 min. Closes the `--engine ocsmesh` deprecation case: the draft niche it filled is now better served by `--om-max-iter 10` (or 25 for "fast production") because the quality gap is far larger than the wall-clock advantage. |
| 37 | `37_phase_e_filter_sweep_poc.py` | `--min-channel-elements` sweep on the same two meshes PoC #35 used (cleaned PoC #19 + PoC #29 centroid-widened). On the cleaned mesh: at min_n=1 the medial-axis estimate is 1.48× the centroid-widen cost (mean long_axis/h_local 1.82); at min_n=3 the ratio falls to 0.94 (mean ratio 3.06); at min_n=10 it is 0.66 (mean ratio 6.51); at min_n=20 it is 0.54 (mean ratio 11.14). **Decision: Stage 2 GO** — once the small-cluster noise is filtered out (min_n ≥ 3), the surviving channels are genuinely ribbon-like and medial-axis insertion is both cheaper *and* topologically more correct than centroid widen. Production sweet spot is `--min-channel-elements 10` (51 components / 1,070 elements / la/h 6.5). PoC #35's "Stage 2 not justified" headline came from the unfiltered noise dominating the average. |
| 38 | `38_phase_e_medial_validation_poc.py` | End-to-end validation of `repair_under_resolved_channels(mode='medial')` (Stage 2) on the cleaned PoC #19 mesh against `widen` at the same `min_channel_elements=10` filter. medial: NP +206 / NE +274 / α 0.9576 → 0.9551 / frac<20° 0.17 % → 0.33 %; widen: NP +1,070 / NE +2,140 / α 0.9576 → 0.9332 / frac<20° 0.17 % → 1.42 %. Stage 2 uses 5–8× fewer new nodes/elements and damages quality 8× less for the same filter. 40/51 components replaced; 11 conservative skips (8 non-convex rim, 3 branching rim) keep the original triangulation in those patches as a safe fallback. The detector-6 residual gap (medial 947 vs widen 862) is by design — medial preserves channel width so the detector keeps reporting "narrow" even after the channel has 2 cells across; alpha / frac<20° are the right quality signals. |
| 39 | `39_courant_sizing_poc.py` | End-to-end validation of `--om-courant-sizing` on Tokyo Bay (PoC #19 baseline vs Courant `dt=10 s` / `C=0.7` / `nu=2 m`). NP +58 % / NE +69 %; alpha mean 0.9586 → 0.9731 (+1.5 %); alpha p05 0.8706 → 0.9171 (+5 %); min angle p05 39.8° → 43.6° (+3.7°); frac<20° 0.101 % → 0.023 % (-77 % rel); n_overconnected 3 → 0; max valence 9 → 8. Worst-case CFL dt p01 @C=0.7 only goes 1.44 → 1.62 s (+12 %) because Tokyo Bay's bottleneck cells sit in coastline-pinned shallow channels where the Courant envelope clamps to `--hmin`; **median dt p50 jumps 4.92 → 7.95 s (+62 %)**. Closes `docs/engine_complementarity.md` decision #5 (`add_topo_bound_constraint` / `add_subtidal_flow_limiter` remain on the follow-up list). |
| 40 | `40_phase_h_dry_run_poc.py` | Stage 1 of the planned per-element greedy optimizer ("Phase H"). On the pipeline-rung-1 output of PoC #19 (47,426 elements), 12,440 (26.2 %) fail the strict per-element gate `alpha >= 0.95` ∧ `min_angle >= 20°`. Dry-run on each fail element of (a) `smooth_node` to the 1-ring centroid and (b) `edge_swap` of each internal edge — re-run with the `_ring_by_node` bug fix: **9,297 (74.7 %) fixable by smooth, 2 by swap, 3,141 (25.3 %) unfixable by either**. Of the unfixable, 2,647 (84 %) touch a boundary, so a boundary-tangent smoother and coastline-aware `edge_split` are the highest-leverage v2 operators. PoC drives the operator-inventory decision for Phase H. (An earlier draft used `np.tile` instead of `np.repeat` to build the node→element ring map, which produced bogus rings and pushed the "smooth_node fixable" count to 0; the headline above is the corrected number.) |
| 41 | `41_phase_h_v1_validation_poc.py` | End-to-end validation of `phase_h_optimize` (Phase H v1) on the pipeline-rung-1 output. The driver alternates Pass A (batch Gauss-Seidel smooth, all interior nodes per sweep, accept iff strict per-1-ring penalty drop) with Pass B (per-element greedy on fail elements with `edge_swap` / `edge_split_interior` / `vertex_remove`); aux dicts are built once per smooth phase and rebuilt per topology accept. Result on the cleaned PoC #19 mesh (4 outer rounds, 61 smooth sweeps, 595 s wall): NP +65 / NE -926; alpha mean 0.9588 → 0.9655 (+0.0067); alpha p05 0.8758 → 0.9001 (+0.0243); min angle p05 40.2° → 41.9° (+1.7°); frac<20° 0.131 % → 0.013 % (-90 % rel). Operators applied: smooth_node 56,047 / vertex_remove 528 / edge_split_interior 65 / edge_swap 2. Fail count drops 12,440 → 11,182 (-10 %); the residual is dominated by boundary-touching elements that v1's conservative boundary handling cannot edit. v2 will add the boundary-tangent smooth + coastline-projecting boundary edge_split that PoC #40 identified as the unfixable-fixer. |
| 42 | `42_phase_h_v2_validation_poc.py` | Phase H v2 (boundary-aware) on the same input. Pass A's batch Gauss-Seidel smooth now also moves segment-interior boundary nodes along the prev-next tangent line (1-ring centroid projected to the line, clamped to (5 %, 95 %) of the segment to avoid collapse onto a neighbour). Pass B gains an `edge_split_boundary` operator that inserts a midpoint on a topological boundary edge, splits the single incident triangle into two, and threads the new node into the relevant open / land segment. Segment endpoints / corners stay pinned. Result on the same pipeline-rung-1 input (4 outer rounds, 72 smooth sweeps, 748 s wall): NP +205 / NE -676; alpha mean 0.9588 → 0.9675 (+0.0087); alpha p05 0.8758 → 0.9084 (+0.0326); min angle p05 40.2° → 42.7° (+2.5°); **frac<20° 0.131 % → 0.000 % — every element now meets the 20° gate**. Fail count drops 12,440 → 10,302 (-17 %, ~2× v1's reduction). Operators applied: smooth_node 57,585 / edge_split_boundary 108 / edge_split_interior 97 / vertex_remove 489 / edge_swap 2. The remaining 10,302 are all alpha < 0.95 cases that single-vertex moves cannot lift without dragging a neighbour below threshold — they are fundamental constraints of greedy local edits, not v2 gaps. v3 will add coastline-shapefile projection for new boundary nodes (currently lands at the straight midpoint). |
| 43 | `43_phase_h_v3_validation_poc.py` | Phase H v3 (coastline-projecting boundary ops) on the same input + the MLIT C23 Tokyo Bay coastline. `build_coastline_projector` loads every polyline into a `shapely.STRtree` and exposes a snap callable (`max_snap_m=500`, ~2.5 × hmin). `_batch_smooth_sweep` and `_apply_edge_split_boundary` route each proposed boundary position through the projector; when within snap range it lands on the nearest polyline, else falls through to the v2 chord midpoint. Result (4 outer rounds, 60 smooth sweeps, 1735 s wall — 2.3× v2 due to ~70k STRtree probes): NP +164 / NE -761; alpha mean 0.9588 → 0.9663 (+0.0076); alpha p05 0.8758 → 0.9040 (+0.0282); min angle p05 40.2° → 42.2° (+2.0°); **frac<20° 0.131 % → 0.000 %** (preserved). Fail 12,440 → 11,015 (-11 %). Operators: smooth_node 56,142 / edge_split_boundary 71 (v2: 108) / edge_split_interior 93 / vertex_remove 509 / edge_swap 2. The lower split_boundary count (71 vs v2's 108) and modest quality drop vs v2 reflect that projecting the midpoint off the chord can leave the two sub-triangles geometrically less symmetric — but every new boundary node now lies on the *actual* coastline, matching the SMS manual-edit standard for coastline fidelity. The user's choice is "v2 for absolute α / fail metrics; v3 for coastline-faithful boundaries". |

`docs/python_pipeline_gap_analysis.md` summarises what the Python
pipeline still has to gain to match the OceanMesh2D reference output.

`docs/architecture.md` is the user-facing decision tree: when to use
`--engine oceanmesh` vs. `--engine ocsmesh`, when to use which
`fmesh-mesh-combine --strategy`, and where the modules live.

## Command-line tools

Installed when `pip install -e .` is run.

| CLI | Purpose |
|-----|---------|
| `fmesh-buildmesh DEM out.14 [--engine oceanmesh\|ocsmesh]` | Single-shot DEM → fort.14. Default engine `oceanmesh` (OceanMesh2D Python port; alpha~0.96). `--engine ocsmesh` is **deprecated** (alpha~0.85, max valence 26; PoC #30 ruled out a Triangle replacement) and slated for removal — see `docs/engine_complementarity.md`. Shared post-processing: depth interp, bbox-based open/land split, river inflow, perpfix. |
| `fmesh-perpfix in.14 out.14`  | Stand-alone open-boundary first-ring perpendicularity correction. |
| `fmesh-subset-dem SRC OUT --bbox MINLON MINLAT MAXLON MAXLAT [--src-var z]` | Clip a global DEM (SRTM15+, GEBCO, GeoTIFF, ...) to a lon/lat bbox and emit a CF-tagged GeoTIFF for `fmesh-buildmesh`. Two read paths: rasterio (CRS-tagged inputs) and netCDF4 (lon/lat NetCDF without CRS, selected by `--src-var`). |
| `fmesh-mesh-combine in1.14 in2.14 [...] out.14 --strategy {disjoint,overlap,neighbor}` | Combine multiple fort.14 meshes. `disjoint` is pure-numpy concat with full boundary preservation (best for non-overlapping basins). `overlap` and `neighbor` wrap `ocsmesh.ops.combine_mesh` for nested-resolution and edge-snap scenarios respectively. |
| `fmesh-mesh-check fort.14 [--max-nbr-elem N] [--min-thin-chain N] [--min-w-h F] [--min-channel-elements N]` | Detect inadequate FVCOM meshes via seven detectors: disjoint wet-domain components, dead-end elements, thin / thin-chain (1-cell channel) elements, over-connected nodes, open-boundary-unreachable elements, and medial-axis-style under-resolved channels (`width / h < --min-w-h`, default 3); `--min-channel-elements N` (default 1) drops detector-6 flags whose connected component has fewer than N elements. Emits `*_summary.txt`, `*_diag.json` (per-id records with coordinates), and `*_map.png`. No repair. Exit code is non-zero when anything is flagged so the command is usable as a CI gate. |
| `fmesh-mesh-clean in.14 out.14 [--bbox] [--open-merge-coast-gap N] [--thin-chain-mode {widen,delete,none}] [--repair-overconnected-iters N] [--under-resolved-mode {widen,delete,none}] [--repair-skewed-elements] [--smooth-laplacian]` | Repair the safe-to-fix subset of the `fmesh-mesh-check` flags. Phase A prunes disjoint dual-graph components; Phase B iteratively trims degree-1 elements with no open-boundary edge; Phase C (default `widen`) inserts a centroid into every thin-chain element so 1-cell channels become 2-cell, or removes the chain entirely with `--thin-chain-mode delete`; Phase D (off by default) runs valence-balancing edge swaps that drive every node valence to at most `--max-nbr-elem`; Phase E (off by default) widens or deletes detector-6 under-resolved channel elements via the same centroid-insertion mechanism; Phase F (off by default) deletes triangles whose interior angles fall outside `[--repair-skewed-min-angle-deg, --repair-skewed-max-angle-deg]` via `ocsmesh.utils.cleanup_skewed_el`; Phase G (off by default) Laplacian-smooths interior nodes via `oceanmesh.laplacian2`. Boundaries are re-derived via DEM-bbox proximity, matching `fmesh-buildmesh`. |
| `fmesh-mesh-quality in.14 [in2.14 ...] [--labels ...] [--min-alpha F] [--max-frac-lt-20deg F] [--max-valence N] [--max-overconnected N] [--max-flipped N] [--max-disjoint-elems N]` | Compute unified mesh-quality metrics (`alpha_mean`, `alpha_p05/p50`, `min_angle_p05/p50_deg`, `frac_lt_20deg`, `max_valence`, `n_overconnected`, `n_flipped`, `n_components`, `n_disjoint_elems`) for one or more fort.14 files. Two inputs print a `delta` column. Threshold flags are evaluated against the LAST input and turn the command into a CI gate (exit 1 on failure). |
| `fmesh-mesh-pipeline in.14 out.14 [--bbox] [--max-iters N] [--best-rung] [--min-alpha F] [--max-frac-lt-20deg F] [--max-valence N] [--max-flipped N] ...` | Progressive `clean → quality → repeat` loop. Applies three cumulative rungs of `fmesh-mesh-clean` phases — rung 0 (A+B+C), rung 1 (+D+F+G), rung 2 (+E) — evaluating `fmesh-mesh-quality` thresholds after each. Stops at the first passing rung by default; with `--best-rung`, runs every rung up to `--max-iters` and outputs the gate-passing rung with the highest `alpha_mean` (ties broken in favour of the lighter repair). Exits 1 if no rung satisfies the gate when thresholds are supplied. JSON history records per-rung metrics and threshold-check results. |

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md) for breaking changes and migration
notes between revisions.

## Status and scope

This project covers the FVCOM-specific mesh tooling that was originally
considered for [`oceanmesh-tools`](https://github.com/estuarine-utokyo/oceanmesh-tools).
`oceanmesh-tools` continues to focus on OceanMesh2D (MATLAB) input scanning
and `fort.14` post-processing.

## License

Apache License 2.0. See `LICENSE` and `NOTICE`. Third-party backend licenses
are documented in `THIRD_PARTY_NOTICES.md`.
