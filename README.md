# fvcom-mesh-tools

> ⚠️ **Pre-alpha** — under active development. APIs and CLIs are unstable.

Python toolkit for FVCOM unstructured mesh generation, repair, and visualization.

`fvcom-mesh-tools` provides a unified Python interface for building
high-quality FVCOM-ready unstructured meshes (`fort.14`), with a focus on:

- **Open-boundary edge orthogonalization** — enforce edges perpendicular to open boundaries
- **River-channel connectivity repair** — fix narrow channels where flow does not connect
- **Mesh quality inspection and visualization** — element quality, boundary classification, fort.14 plots

The package wraps several mature mesh tools behind a common backend interface
rather than reimplementing meshing algorithms from scratch.

## Backend strategy

| Role | Backend | License | How used |
|------|---------|---------|----------|
| Mesh generation (preferred) | [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 | imported |
| Orthogonalization / smoothing | [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | imported |
| Grid utilities | [stompy](https://github.com/rustychris/stompy) | MIT | imported |
| Geometry-aware sizing (optional) | [JIGSAW](https://github.com/dengwirda/jigsaw-python) | custom (non-OSI) | imported, optional extra |
| External mesher | [gmsh](https://gmsh.info/) | GPL-2.0+ | invoked as subprocess |

GPL-licensed backends are deliberately invoked as subprocesses to keep this
project's distribution under Apache-2.0. See `THIRD_PARTY_NOTICES.md`.

## Installation

Python ≥3.12. conda-forge is the recommended channel for the scientific stack:

```bash
mamba env create -f environment.yml
mamba activate fvcom-mesh
pip install -e .
```

Optional backends are pulled in via extras:

```bash
pip install -e ".[ocsmesh,meshkernel,test]"
```

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

`docs/python_pipeline_gap_analysis.md` summarises what the Python
pipeline still has to gain to match the OceanMesh2D reference output.

## Command-line tools

Installed when `pip install -e .` is run.

| CLI | Purpose |
|-----|---------|
| `fmesh-buildmesh DEM out.14 [--engine oceanmesh\|ocsmesh]` | Single-shot DEM → fort.14. Default engine `oceanmesh` (OceanMesh2D Python port; alpha~0.96, slow). Alternative engine `ocsmesh` (OCSMesh+gmsh; alpha~0.85, ~40× faster) for draft / iteration. Shared post-processing: depth interp, bbox-based open/land split, river inflow, perpfix. |
| `fmesh-perpfix in.14 out.14`  | Stand-alone open-boundary first-ring perpendicularity correction. |
| `fmesh-subset-dem SRC OUT --bbox MINLON MINLAT MAXLON MAXLAT [--src-var z]` | Clip a global DEM (SRTM15+, GEBCO, GeoTIFF, ...) to a lon/lat bbox and emit a CF-tagged GeoTIFF for `fmesh-buildmesh`. Two read paths: rasterio (CRS-tagged inputs) and netCDF4 (lon/lat NetCDF without CRS, selected by `--src-var`). |
| `fmesh-mesh-combine in1.14 in2.14 [...] out.14 --strategy {disjoint,overlap,neighbor}` | Combine multiple fort.14 meshes. `disjoint` is pure-numpy concat with full boundary preservation (best for non-overlapping basins). `overlap` and `neighbor` wrap `ocsmesh.ops.combine_mesh` for nested-resolution and edge-snap scenarios respectively. |

## Status and scope

This project covers the FVCOM-specific mesh tooling that was originally
considered for [`oceanmesh-tools`](https://github.com/estuarine-utokyo/oceanmesh-tools).
`oceanmesh-tools` continues to focus on OceanMesh2D (MATLAB) input scanning
and `fort.14` post-processing.

## License

Apache License 2.0. See `LICENSE` and `NOTICE`. Third-party backend licenses
are documented in `THIRD_PARTY_NOTICES.md`.
