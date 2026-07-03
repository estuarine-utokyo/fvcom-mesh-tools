# TOOLCHAIN_SURVEY — Phase 0 reconnaissance

Surveyed 2026-07-03 per `fvcom_mesh_kickoff.md` §6. Companion documents:
`docs/DATA_INVENTORY.md` (input datasets), `docs/fvcom_source_constraints.md`
(FVCOM-source-derived constraint set, §6 item 6), `docs/architecture.md`
(engine decision tree), `docs/python_pipeline_gap_analysis.md` (PoC history).

Verdict up front: **all mesh-generation and mesh-repair primitives the kickoff
needs already exist across the four toolkits; what is missing is the integration
layer** — recipe-driven orchestration, FVCOM native writers, one unified QA gate,
and the reference-grid figure system — plus four ports from OceanMesh2D
(CFL limiter, channel sizing, `interp`-equivalent, `make_bc`-equivalent).

---

## 1. Toolkit inventory

### 1.1 oceanmesh (Python, GPL-3.0) — the generator

Library only (no CLI). Modules: `edgefx.py` (sizing), `signed_distance_function.py`
(polygon/CSG domains), `geodata.py` (`DEM`, `Shoreline`), `clean.py` + `fix_mesh.py`
(cleaning), `mesh_generator.py` (`generate_mesh`, `generate_multiscale_mesh`,
`write_to_fort14`), `boundary.py` (`identify_ocean_boundary_sections`).

Sizing functions available (`edgefx.py`): `distance_sizing_function` (:356),
`distance_sizing_from_line_function` (:196), `distance_sizing_from_point_function`
(:280), `feature_sizing_function` (medial-axis, :752), `wavelength_sizing_function`
(:866), `bathymetric_gradient_sizing_function` (:429, with
`rossby_radius_filter`), `multiscale_sizing_function` (:947),
`enforce_mesh_gradation` (:91), `enforce_mesh_size_bounds_elevation` (:36).
Note: `__init__.py:110` exports a nonexistent `slope_sizing_function` (dead alias —
slope sizing is `bathymetric_gradient_sizing_function`).

Cleaning: `mesh_clean`, `make_mesh_boundaries_traversable` (the >2-boundary-edges
repair), `delete_{exterior,interior,boundary}_faces`,
`delete_faces_connected_to_one_face`, `laplacian2`, `fix_mesh`/`simp_qual`.

| Requirement | Status |
|---|---|
| distance / feature / wavelength / gradient sizing | exists |
| thalweg-polyline + min-size channel sizing | partial (generic line-distance sizer; not depth/width-aware) |
| FVCOM outputs | partial (ADCIRC fort.14 + t3s only) |
| recipe-driven config, QA gate, OBC tooling, figures | missing |

### 1.2 oceanmesh-tools — OM2D scanning + fort.14 post-processing

CLIs: `omt scan` (MATLAB-script input discovery → catalog JSON), `omt viz`
(fort.14 + DEM/coastline figures), `omesh14-view` (Plotly interactive),
`omesh14-edit-bdy` / `omesh14-fixstyle` (fort.14 boundary editing, byte-style
preserving). API: full fort.14 parser incl. boundary groups (`io/fort14.py`),
boundary-loop topology (`mesh/boundary.py`: `classify_open_boundary_edges`,
`walk_closed_loops`, …), matplotlib rendering (`plot/viz.py`).

| Requirement | Status |
|---|---|
| fort.14 parse + boundary segment tooling | exists |
| OBC classification | partial (labeling; no perpendicularity metric/fix) |
| land rendering / reference grid | partial (coastline overlay; lat/lon grid only, no 5 km metric grid) |
| FVCOM writers, QA gate, recipe | missing |

### 1.3 xcoast (MIT) — land/coast rendering + OSM

Library only. `CoastMask` + `load(name_or_bbox)` (cache-first; presets include
`tokyo_bay`, `tokyo_bay_inner`, `osaka_bay`, …); land = Geofabrik/government
polygons minus inland water, or OSM Overpass (`overpass.py`). Rendering:
`add_to_mpl` (Cartopy) and **`add_to_plain_mpl` / `add_land_to_plain_axes`
(plain matplotlib Axes, `zorder`-controllable)** — land behind a mesh triplot is
achievable today by passing a low zorder.

| Requirement | Status |
|---|---|
| land polygons + OSM/government coastline rendering | exists |
| land behind mesh in matplotlib | exists (zorder) |
| labeled 5 km reference grid (A1/B2/C4) | missing |

### 1.4 fvcom-mesh-tools (this repo, Apache-2.0) — FVCOM I/O + QA + repair

CLIs (all flag-driven): `fmesh-buildmesh` (DEM → engine → depth → boundaries →
rivers → quality → perpfix → fort.14, ~45 flags), `fmesh-perpfix` (first-ring
perpendicularity fix; a *fix*, not a placement tool — placement is
`algorithms/boundary.classify_boundaries_by_bbox`), `fmesh-subset-dem`,
`fmesh-mesh-combine`, `fmesh-mesh-check` (7 topology detectors incl. wet
connectivity + valence ≤ 8), `fmesh-mesh-clean` (phases A–G),
`fmesh-mesh-quality` (metrics + threshold gate), `fmesh-mesh-pipeline`
(progressive rungs 0–3), plus `mesh_clean_phase_h.py` (`phase_h_optimize`,
`phase_h_finish` — FVCOM strict gates C1/C2/C4/C5 as a repair loop).

FVCOM outputs today: **fort.14 only** (`io/fort14.py`, round-trip safe, rivers
ibtype 20/21). No `_grd.dat` / `_dep.dat` / `_obc.dat` / `_cor` / `_spg` / `.2dm`
writers anywhere. `plotting.py` is an empty stub (one DPI constant);
`docs/architecture.md` §5's `vis/`/`plot/`/`scan/` subpackages do not exist (stale
doc — fix when touched).

QA coverage vs the kickoff §9 list:

| §9 check | Status | Where |
|---|---|---|
| min angle ≥ 30° / max angle ≤ 130° / area-change ≤ 0.5 | partial | Phase H gates + notebook #48 only; not in a read-only QA CLI |
| valence ≤ 8 | exists | diagnostics + `fmesh-mesh-quality --max-valence` |
| Delaunay / overlap / flipped | partial | `n_flipped` gate; Delaunay only inside Phase H retriangulation |
| OBC perpendicularity | partial | metric + fix exist; no pass/fail gate |
| no dual-boundary-edge element (R4) | missing as explicit check | approximated by dead-end/thin detectors |
| wet connectivity | exists | `disjoint_components` / `unreachable_elements` |
| min-depth compliance | missing | — |
| implied Δt | missing (generation-side sizing only) | — |
| single command, JP pass/fail table | missing | split across 3 surfaces |

---

## 2. MATLAB reference workflow (OceanMesh2D/Tokyo_Bay)

Canonical script: `mesh_wide_futtsu_5r.m` (5-region multiscale; 5 reduced variants
share the same pipeline). Pipeline: per-region `safe_geodata('shp',…,'dem',…,
'h0',min_el,'bbox',…)` → `edgefx('fs',R,'wl',wl,'slp',slp,'fl',fl,'max_el',…,
'dt',dt,'g',grade)` → `meshgen('ef',{...},'bou',{...},'proj','trans','itmax',50)`
→ `m = interp(m, bou_list, 'mindepth', 0.01)` → `make_bc(m,'auto',gdat_01,'both')`
→ `write` (fort.14).

Region parameters (mesh_wide_futtsu_5r.m):

| Region | min_el | max_el | wl | slp | fl | grade | fs R | coastline | DEM |
|---|---|---|---|---|---|---|---|---|---|
| R1 Pacific | 10e3 | 50e3 | 30 | 50 | −50 | 0.25 | 3 | GSHHS_f_L1 | SRTM15+ slice |
| R2 entrance | 6e2 | 20e3 | 30 | 50 | −50 | 0.25 | 3 | C23-06_TOKYOBAY | SRTM15+ slice |
| R3 outer bay | 1e2 | 9e2 | 10 | 100 | −200 | 0.35 | 6 | C23 _OUTER | depth_0090 nc |
| R4 inner bay | 1e2 | 9e2 | 10 | 50 | −10 | 0.25 | 3 | C23 _INNER | depth_0090 nc |
| R5 Futtsu | 8 | 50 (ns 30) | 10 | 50 | −10 | 0.25 | 3 | Futtsu_coastline | depth_0090 nc |

Example patterns to reuse:

- **JBAY** (`Example_5_JBAY.m`): CFL limiting is just `edgefx(...,'dt',2,...)` —
  enlarges elements wherever a 2 s external timestep would violate CFL at local
  depth. Tokyo scripts use `dt=0` (auto from nearshore resolution). Post:
  `interp(...,'nan','fill','mindepth',1)`, `make_bc(...,'depth',5)` (OBC depth
  cutoff).
- **GBAY** (`Example_6_GBAY.m`): channel sizing is
  `edgefx(...,'Channels',pts2,'ch',0.1,...)` with thalweg polylines; `Example_9`
  shows `enforceMin` carrying the finest edgefunction (incl. channel sizing)
  inward through nested regions — the default multiscale behavior the Tokyo
  scripts rely on.

Caveat from the data audit: the `Futtsu_coastline.shp` / `coastline_2.shp` the
reference scripts consume have continental-scale extents (suspect stray
vertices) — verify/clean before reuse in the Python pipeline.

## 3. MATLAB → Python mapping (gaps = port candidates)

| OceanMesh2D feature | Python oceanmesh equivalent | Status |
|---|---|---|
| `edgefx 'fs'` feature sizing | `feature_sizing_function(r=…)` | OK |
| `edgefx 'wl'` wavelength | `wavelength_sizing_function` | OK |
| `edgefx 'slp'`/`'fl'` slope + Rossby filter | `bathymetric_gradient_sizing_function(slope_parameter, filter_quotient)` | OK-ish — verify slp/fl ↔ parameter semantics before relying on values |
| `edgefx 'g'` gradation | `enforce_mesh_gradation` | OK |
| **`edgefx 'dt'` CFL limiter** | none (repo has `courant_sizing_function` in `mesh_engine/oceanmesh.py` — reconcile/upstream) | **GAP #1** |
| **`edgefx 'Channels'/'ch'` thalweg sizing** | only generic `distance_sizing_from_line_function` | **GAP #2** |
| `interp(m, gdat, 'mindepth', 'nan fill')` | none at msh level (repo has `dem/interp.at_points` — partial) | **GAP #3** |
| `make_bc('auto'/'both'/'depth')` | `identify_ocean_boundary_sections` + repo bbox classifier — no depth-cutoff auto placement | **GAP #4** |
| `max_el_ns` nearshore cap | none | minor gap |
| multiscale nesting + `enforceMin` | `generate_multiscale_mesh` + `multiscale_sizing_function` | OK |
| `meshgen` proj='trans'/'utm' | crs/stereo model differs; **UTM 54N workflow not yet modeled** | gap (recipe layer) |
| `write` fort.14 | `write_to_fort14` / repo `write_fort14` | OK |

## 4. FVCOM source-derived constraints (§6 item 6) — summary

Full report: `docs/fvcom_source_constraints.md`. Headlines:

- Production = FVCOM 5.1, `-DGCN`, CARTESIAN, no PLBC (the R4 stop is active).
- **R4 confirmed and stronger than the kickoff prose**: `tge.F:558-581` PSTOPs any
  element whose three node flags sum > 4 — an open edge plus *any* third boundary
  node (solid or open) is fatal, not just an open+solid edge pair.
- `_obc.dat` order is not required to be spatially sequential, but every OBC node
  must be mesh-adjacent to another OBC node (connected chains; `PSTOP` otherwise),
  and junction nodes must be in the OBC list.
- `_grd.dat` must be CCW; FVCOM checks element #1 only and swaps columns to
  internal CW on read — QA must check all elements.
- FVCOM never checks: duplicate nodes, orphan nodes, flipped elements beyond #1,
  valence (silent scratch-array overrun above 50), or any angle/area-ratio
  quality metric. The C1–C5 gates are project-mandatory, not code-enforced.

## 5. Consolidated gap list → work items

Integration layer (Phase 1, priority order):

1. **Unified QA gate** (`fmesh-mesh-qa <mesh>`): §9 checklist = existing
   diagnostics + quality gates + Phase-H C1/C2/C4/C5 as read-only checks + new
   R4 / OBC-chain / CCW-all / duplicate-node / orphan-node / min-depth /
   implied-Δt checks from `fvcom_source_constraints.md`; single pass/fail table
   (Japanese output per kickoff §9), non-zero exit on fail. This unblocks the
   `/loop` validator.
2. **FVCOM native writers**: `casename_grd.dat` / `_dep.dat` / `_obc.dat`
   (+ `_cor`/`_spg` hooks) and SMS `.2dm` export, from `Fort14Mesh` — formats per
   `fvcom_source_constraints.md` §file-formats.
3. **recipe.yaml layer**: schema capturing CRS (UTM 54N), datum shifts per source
   (M7001 CD +1.13 m!), source precedence, region polygons, sizing stack,
   boundary treatment, quality targets, seed; a runner mapping recipe →
   existing buildmesh/pipeline calls; provenance record emission.
4. **Figure system**: xcoast land (plain-Axes, low zorder) + 5 km UTM reference
   grid overlay (A/B/C… × 1/2/3…, named-region aliases) + per-cell zoom panels;
   region-reference parser (cell/alias → UTM polygon) for prose addressing.

Ports from OceanMesh2D (Phase 1–2):

5. CFL/`dt` edge limiter (GAP #1 — reconcile with the repo's existing
   `courant_sizing_function`).
6. Depth/width-aware channel/thalweg sizing (GAP #2; GBAY pattern).
7. msh-level bathy `interp` with per-source datum shift + `mindepth` + nan-fill
   (GAP #3) — this is where the DATA_INVENTORY precedence/datum rules execute.
8. `make_bc`-style auto OBC placement with depth cutoff (GAP #4) + OBC
   perpendicularity as a *gate* (fix already exists).

Minor: `max_el_ns` nearshore cap; robust-shapefile loader (`safe_geodata`
equivalent — needed anyway for the broken Futtsu shapefiles); fix stale
`architecture.md` §5; remove dead `slope_sizing_function` export upstream.
