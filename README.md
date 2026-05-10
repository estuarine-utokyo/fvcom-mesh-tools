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
| Mesh generation (alt / draft) | [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) + [gmsh](https://gmsh.info/) | CC0-1.0 + GPL-2.0+ runtime | imported (OCSMesh); gmsh called via OCSMesh |
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
| `pip install -e ".[ocsmesh]"` | + ocsmesh, gmsh, and the above | `fmesh-buildmesh --engine ocsmesh`, `fmesh-mesh-combine --strategy {overlap,neighbor}` |
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

To iterate fast, swap to the OCSMesh + gmsh backend (~40x faster, lower
quality). For draft work the same flag set works:

```bash
fmesh-buildmesh /tmp/tb.tif /tmp/tokyo_draft.14 --engine ocsmesh \
    --hmin 200 --hmax 5000 \
    --coastline data/coastline/MLIT_C23/C23-06_TOKYOBAY.shp \
    --coast-target-size 200 --coast-expansion-rate 0.005 \
    --min-island-area-m2 100000 \
    --quality-pass 6 --refine-min-angle 20
```

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
enable Phase E deliberately.

`docs/architecture.md` is the full decision tree for engine choice and
combine strategy; `docs/python_pipeline_gap_analysis.md` has the
quality / runtime numbers vs. the OceanMesh2D MATLAB reference.

## Development

```bash
make install   # editable install with dev deps
make test      # pytest -q
make lint      # ruff check
```

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

`docs/python_pipeline_gap_analysis.md` summarises what the Python
pipeline still has to gain to match the OceanMesh2D reference output.

`docs/architecture.md` is the user-facing decision tree: when to use
`--engine oceanmesh` vs. `--engine ocsmesh`, when to use which
`fmesh-mesh-combine --strategy`, and where the modules live.

## Command-line tools

Installed when `pip install -e .` is run.

| CLI | Purpose |
|-----|---------|
| `fmesh-buildmesh DEM out.14 [--engine oceanmesh\|ocsmesh]` | Single-shot DEM → fort.14. Default engine `oceanmesh` (OceanMesh2D Python port; alpha~0.96, slow). Alternative engine `ocsmesh` (OCSMesh+gmsh; alpha~0.85, ~40× faster) for draft / iteration. Shared post-processing: depth interp, bbox-based open/land split, river inflow, perpfix. |
| `fmesh-perpfix in.14 out.14`  | Stand-alone open-boundary first-ring perpendicularity correction. |
| `fmesh-subset-dem SRC OUT --bbox MINLON MINLAT MAXLON MAXLAT [--src-var z]` | Clip a global DEM (SRTM15+, GEBCO, GeoTIFF, ...) to a lon/lat bbox and emit a CF-tagged GeoTIFF for `fmesh-buildmesh`. Two read paths: rasterio (CRS-tagged inputs) and netCDF4 (lon/lat NetCDF without CRS, selected by `--src-var`). |
| `fmesh-mesh-combine in1.14 in2.14 [...] out.14 --strategy {disjoint,overlap,neighbor}` | Combine multiple fort.14 meshes. `disjoint` is pure-numpy concat with full boundary preservation (best for non-overlapping basins). `overlap` and `neighbor` wrap `ocsmesh.ops.combine_mesh` for nested-resolution and edge-snap scenarios respectively. |
| `fmesh-mesh-check fort.14 [--max-nbr-elem N] [--min-thin-chain N] [--min-w-h F]` | Detect inadequate FVCOM meshes via seven detectors: disjoint wet-domain components, dead-end elements, thin / thin-chain (1-cell channel) elements, over-connected nodes, open-boundary-unreachable elements, and medial-axis-style under-resolved channels (`width / h < --min-w-h`, default 3). Emits `*_summary.txt`, `*_diag.json` (per-id records with coordinates), and `*_map.png`. No repair. Exit code is non-zero when anything is flagged so the command is usable as a CI gate. |
| `fmesh-mesh-clean in.14 out.14 [--bbox] [--open-merge-coast-gap N] [--thin-chain-mode {widen,delete,none}] [--repair-overconnected-iters N] [--under-resolved-mode {widen,delete,none}]` | Repair the safe-to-fix subset of the `fmesh-mesh-check` flags. Phase A prunes disjoint dual-graph components; Phase B iteratively trims degree-1 elements with no open-boundary edge; Phase C (default `widen`) inserts a centroid into every thin-chain element so 1-cell channels become 2-cell, or removes the chain entirely with `--thin-chain-mode delete`; Phase D (off by default) runs valence-balancing edge swaps that drive every node valence to at most `--max-nbr-elem`; Phase E (off by default) widens or deletes detector-6 under-resolved channel elements via the same centroid-insertion mechanism. Boundaries are re-derived via DEM-bbox proximity, matching `fmesh-buildmesh`. |

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
