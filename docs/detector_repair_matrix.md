# Detector → repair → quality matrix

This is the canonical map from a problem flagged by
``fmesh-mesh-check`` to the ``fmesh-mesh-clean`` phase that fixes it,
the ``fmesh-mesh-quality`` metric that measures it, and the
``fmesh-mesh-pipeline`` rung that turns it on automatically.

It exists because the same concept (e.g. *over-connected node*)
appears under three names in three commands, and a new contributor
shouldn't have to grep the codebase to learn the mapping.

If you're looking for *"my mesh has X problem; what do I run?"*,
read [§ 1 — quick lookup table](#1-quick-lookup-table). For
*"how do these CLIs fit together?"*, read [§ 5 — recommended
workflow](#5-recommended-workflow).

## 1. Quick lookup table

| Problem (detector) | Detector function | Repair phase | `fmesh-mesh-clean` flag | Quality metric | Pipeline rung |
| --- | --- | --- | --- | --- | --- |
| Disjoint dual-graph components | `disjoint_components_flag` | **A** — `keep_components` | (default ON) `--no-remove-disjoint` to skip | `n_components`, `n_disjoint_elems` | rung 0 |
| Dead-end (degree-1) elements | `dead_end_elements_flag` | **B** — `trim_dead_ends` | `--trim-dead-ends-iters N` (default 10; `0` to skip) | (NE drop only) | rung 0 |
| Thin elements (3 vertices on boundary) | `thin_elements_flag` | (filter only — used by Phase C) | — | — | rung 0 |
| Thin chains (1-cell-wide channels) | `thin_chain_elements_flag` | **C** — `repair_thin_chains` (widen-at-centroid by default) | `--thin-chain-mode {widen,delete,none}` (default `widen`) | — | rung 0 |
| Over-connected nodes (valence > N) | `overconnected_nodes_flag` | **D** — `repair_overconnected_nodes` (Lawson edge swap) | `--repair-overconnected-iters N` (default `0` = OFF) | `max_valence`, `n_overconnected` | **rung 1** |
| Under-resolved channels (medial-axis `w/h < min_w_h`) | `under_resolved_channels_flag` | **E** — `repair_under_resolved_channels` (centroid widen, deletion, or medial-axis CDT re-mesh) | `--under-resolved-mode {widen,delete,medial,none}` (default `none`); `--under-resolved-min-channel-elements N` filters out small isolated clusters (default `1` = no filter; PoC #37 sweet spot for `medial` is `10`) | — | **rung 2** |
| Open-boundary-unreachable elements | `unreachable_elements_flag` | (no dedicated phase — overlaps with Phase A) | (use `--require-open-boundary` on Phase A) | `n_components` | rung 0 (incidental) |
| Skewed triangles (interior angle out of `[θmin, θmax]`) | (geometric, not a detector) | **F** — `repair_skewed_elements` (wraps `ocsmesh.utils.cleanup_skewed_el`) | `--repair-skewed-elements`, `--repair-skewed-min-angle-deg`, `--repair-skewed-max-angle-deg` | `min_angle_p05_deg`, `frac_lt_20deg` | **rung 1** |
| Locally-irregular interior nodes | (geometric, not a detector) | **G** — `smooth_mesh_laplacian` (wraps `oceanmesh.laplacian2`) | `--smooth-laplacian`, `--smooth-laplacian-iters`, `--smooth-laplacian-tol` | `alpha_mean`, `min_angle_p05_deg`, `frac_lt_20deg` | **rung 1** |

The detector list lives in
[`src/fvcom_mesh_tools/diagnostics.py`](../src/fvcom_mesh_tools/diagnostics.py);
the phase list lives in
[`src/fvcom_mesh_tools/mesh_clean.py`](../src/fvcom_mesh_tools/mesh_clean.py);
the metric list lives in
[`src/fvcom_mesh_tools/quality.py`](../src/fvcom_mesh_tools/quality.py)
(`METRIC_KEYS`).

## 2. Phase ordering and dependencies

Phases run in **A → B → C → D → E → F → G** when enabled. The
ordering matters:

| Order | Why |
| --- | --- |
| A → B | Trimming dead-end elements is meaningless until disjoint pools have been removed (an isolated triangle is structurally indistinguishable from a dead-end). |
| B → C | Thin chains often terminate in dead-ends; trimming first makes the chain detection cleaner. |
| C → D | The centroid insertions of Phase C-widen create new edges; running Phase D first would balance valences that are about to change. |
| D → E | Same logic — Phase E inserts more centroids and would re-introduce over-connection if D ran after E. |
| E → F | Phase E's new centroid sub-triangles are sometimes skewy; F mops them up. |
| F → G | Smoothing assumes a locally-valid topology — F removes degenerate triangles that would otherwise pull G's solution off. |

**A, B, C-delete, F all delete elements** and force a re-derivation
of the open / land boundary lists via the DEM-bbox classifier
(``classify_boundaries_by_bbox``). Phase C-widen, D, E-widen, and
G are topology- or coordinate-only changes and do not touch the
boundary edge set.

## 3. Side effects of each phase

| Phase | NP change | NE change | Boundary set | Invariants |
| --- | --- | --- | --- | --- |
| A `keep_components` | ↓ | ↓ | reset (re-derived) | Drops disjoint pools. Default keeps the largest component. |
| B `trim_dead_ends` | ↓ | ↓ | reset | Iteratively trims degree-1 elements with no open-boundary edge. |
| C-widen `repair_thin_chains` | **↑** (one centroid per flagged element) | **↑** (`+2` per flagged element: each becomes 3 sub-triangles) | preserved | Existing boundary node IDs stay valid. |
| C-delete `repair_thin_chains` | ↓ | ↓ | reset | Removes the entire thin chain. Aggressive; rare to want. |
| D `repair_overconnected_nodes` | unchanged | unchanged | preserved | Edge swap only; signed area never goes negative (`min_angle_floor_deg`). |
| E-widen `repair_under_resolved_channels` | **↑** (lots) | **↑** (lots — typically `+2 × O(thousands)`) | preserved | Detector 6 typically flags thousands of elements on real meshes; the output mesh size grows proportionally. |
| E-medial `repair_under_resolved_channels` | **↑** (one spine sample per `h_local_median` of channel length, summed across qualifying components) | **↑** (rim+spine Delaunay triangulation; per PoC #38: ~5-8× fewer elements than E-widen at the same filter) | reset only if `bbox`/`tol_deg` supplied; otherwise the existing boundary lists are preserved by node ID | Per face-face channel >= `min_channel_elements` flagged members. CCW-orients the patch rim, finds spine via two-BFS diameter, samples at `h_local_median` spacing, Delaunay-triangulates `(rim ∪ spine)` and prunes triangles outside the rim polygon. Components whose rim is branching or pathologically non-convex are skipped (their original triangulation is kept) — counted under `info["skip_reasons"]`. Output is flip-free by construction (degenerate / inverted triangles abort the per-component replacement). |
| E-delete | ↓ | ↓ | reset | Strips the flagged elements (rarely what you want). |
| F `repair_skewed_elements` | unchanged or ↓ | unchanged or ↓ | preserved if no-op, reset if any deletion | Wraps `ocsmesh.utils.cleanup_skewed_el`. ocsmesh is library-only here (no gmsh). |
| G `smooth_mesh_laplacian` | unchanged | unchanged | preserved | Boundary nodes auto-pinned. Built-in flipped-triangle safety net (`repair_flipped=True`) reverts any element that goes negative-area. |

## 4. Threshold gate heuristics

These are starting points for `fmesh-mesh-quality` and
`fmesh-mesh-pipeline` thresholds; tune for your specific FVCOM
build's CFL constraint.

| Threshold | Typical floor / ceiling | Source |
| --- | --- | --- |
| `--min-alpha 0.95` | oceanmesh production meshes (PoC #19) | Tokyo Bay PoC #19 reached 0.959 |
| `--max-frac-lt-20deg 0.005` | i.e. ≤ 0.5 % of triangles below 20° | `architecture.md` §2 |
| `--max-valence 8` | Conservative legacy FVCOM `MAX_NBR_ELEM` | match your build's compile-time cap |
| `--max-overconnected 0` | Hard zero is achievable on cleaned oceanmesh meshes (PoC #19, #33) | Phase D + edge swap |
| `--max-flipped 0` | Hard zero — Phase G's repair_flipped safety net guarantees this, and the same safety net now wraps the build-time `om.laplacian2` call inside `fmesh-buildmesh --engine oceanmesh` | added in commits `ba1d948` (Phase G) and the build-pipeline follow-up |
| `--max-disjoint-elems 0` | Hard zero — Phase A removes them | rung 0 always covers this |

A "FVCOM-friendly default preset":

```bash
--min-alpha 0.95 \
--max-frac-lt-20deg 0.005 \
--max-valence 8 \
--max-flipped 0 \
--max-disjoint-elems 0
```

`fmesh-mesh-pipeline` PoC #33 validates this preset end-to-end on the
PoC #19 raw Tokyo Bay mesh: rung 0 fails (max valence stays at 9),
rung 1 (`+D+F+G`) passes.

## 5. Recommended workflow

```
                   ┌─────────────────┐
DEM (NetCDF /      │ fmesh-subset-   │
 GeoTIFF)  ───────►│ dem             │ → DEM clipped to lon/lat bbox
                   └────────┬────────┘
                            ▼
                   ┌─────────────────┐
                   │ fmesh-buildmesh │ → fort.14 (oceanmesh DistMesh)
   coastline ─────►│ --engine        │
   river points    │ oceanmesh       │
                   └────────┬────────┘
                            ▼
                  ┌──────────────────┐
                  │ fmesh-mesh-      │ → flagged-element JSON,
                  │ check            │   summary, map.png
                  └────────┬─────────┘
                           ▼
                ┌────────────────────┐
                │ fmesh-mesh-pipeline│ ← preferred 1-shot loop
                │  (clean + quality) │
                └────────┬───────────┘
                         ▼
                  ┌─────────────┐
                  │ fmesh-mesh- │   ← optional re-check
                  │ quality     │
                  └─────────────┘
```

For a multi-basin model, also reach for ``fmesh-mesh-combine`` after
``fmesh-buildmesh`` and before ``fmesh-mesh-pipeline``.

`fmesh-mesh-pipeline` is preferred over hand-running
`fmesh-mesh-check` → `fmesh-mesh-clean` → `fmesh-mesh-quality`
because:

1. It encodes the "rung-0-then-rung-1-then-rung-2" escalation as a
   single decision tree.
2. It records per-rung metrics in JSON for audit (which rung
   actually fixed which threshold).
3. The threshold gates are evaluated against the **final** mesh, so
   the exit code is meaningful in CI.

## 6. Where to add a new detector / phase

If you find a new mesh pathology, the convention is:

1. Detection — add a `*_flag` function to `diagnostics.py` returning
   a boolean array. Wire it into `run_diagnostics` so
   `fmesh-mesh-check` surfaces it.
2. Repair — if there's a clean local fix, add a `repair_*` function
   to `mesh_clean.py`, plumb it through `clean_mesh` as a Phase H
   (or wherever it fits in the order), and expose a CLI flag.
3. Quality — if the new defect is metric-relevant, extend
   `quality.METRIC_KEYS` and `compute_metrics`.
4. Pipeline — choose a rung. New "always-safe" repairs go in rung 0;
   anything that adds elements goes no earlier than rung 2.
5. **Document here.** Add a row in
   [§ 1 quick lookup table](#1-quick-lookup-table) and a side-effect
   row in [§ 3](#3-side-effects-of-each-phase).

## 7. Cross-references

| Topic | Where |
| --- | --- |
| Engine choice (oceanmesh vs ocsmesh) | [`docs/engine_complementarity.md`](engine_complementarity.md) |
| Architecture decision tree | [`docs/architecture.md`](architecture.md) |
| Quality / runtime gap vs. OceanMesh2D MATLAB reference | [`docs/python_pipeline_gap_analysis.md`](python_pipeline_gap_analysis.md) |
| Per-PoC validation numbers | [README PoC notebook table](../README.md#proof-of-concept-notebooks) |
| Per-detector / per-phase implementation | source files linked from § 1 |
