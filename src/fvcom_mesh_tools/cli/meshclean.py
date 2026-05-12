"""``fmesh-mesh-clean`` CLI: clean an FVCOM mesh in seven composable phases.

Phase A removes dual-graph connected components by size and / or
open-boundary touch (default keeps only the largest component). Phase B
iteratively deletes degree-1 elements with no open-boundary edge so
"spit" terminations of 1-cell channels are removed. Phase C widens or
deletes 1-cell-wide channels (chains of "thin" triangles where all 3
vertices sit on a boundary). Phase D drives node valences down via
Lawson edge-flips. Phase E widens or deletes elements flagged as
under-resolved by the medial-axis channel-width detector
(``w/h < min_w_h``), catching 2- and 3-cell-wide channels Phase C does
not flag. Phase F deletes triangles whose minimum or maximum interior
angle exceeds the configured thresholds (wraps
``ocsmesh.utils.cleanup_skewed_el``; ocsmesh used as a library only,
no gmsh dependency). Phase G smooths interior nodes via Laplacian
relaxation (wraps ``oceanmesh.laplacian2``). Boundaries are re-derived
against a DEM-bbox classifier matching ``fmesh-buildmesh``. Phases D,
E, F, and G are off by default — enable deliberately.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.diagnostics import (
    DEFAULT_ARC_SEPARATION_FACTOR,
    DEFAULT_CHANNEL_SAMPLE_DS_M,
    DEFAULT_MIN_W_H,
    DEFAULT_OPPOSITE_BANK_COS_MAX,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import (
    DEFAULT_BBOX_TOL_M,
    DEFAULT_SKEWED_MAX_ANGLE_DEG,
    DEFAULT_SKEWED_MIN_ANGLE_DEG,
    DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    DEFAULT_SMOOTH_LAPLACIAN_TOL,
    DEFAULT_SMOOTH_REPAIR_PASSES,
    clean_mesh,
)


def _infer_bbox(mesh: Fort14Mesh) -> tuple[tuple[float, float, float, float], str]:
    """Pick a classification bbox from the input mesh.

    If the mesh already has at least one open-boundary segment, use the
    bbox of those OB node coordinates — this preserves the original
    fmesh-buildmesh classification (open boundary on one side of the
    DEM bbox) when re-deriving boundaries after deletion. Otherwise
    fall back to the mesh's node bbox, which assumes every side of
    the mesh is open.
    """
    if mesh.open_boundaries:
        idx = np.unique(np.concatenate([np.asarray(s) for s in mesh.open_boundaries]))
        pts = mesh.nodes[idx]
        return (
            (
                float(pts[:, 0].min()), float(pts[:, 1].min()),
                float(pts[:, 0].max()), float(pts[:, 1].max()),
            ),
            "input open-boundary nodes (auto)",
        )
    return tuple(mesh.bbox), "mesh.bbox (auto)"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-clean",
        description=(
            "Clean an FVCOM mesh in seven composable phases. Phase A: "
            "drop dual-graph components by size / open-boundary touch. "
            "Phase B: trim degree-1 dead-end elements iteratively. "
            "Phase C: widen or delete 1-cell-wide channels. Phase D: "
            "Lawson edge-flips to drive valence below MAX_NBR_ELEM. "
            "Phase E: widen or delete medial-axis-detected "
            "under-resolved channels (2- and 3-cell-wide). Phase F: "
            "delete skewed triangles by angle thresholds (wraps "
            "ocsmesh.utils.cleanup_skewed_el). Phase G: Laplacian "
            "smoothing of interior nodes (wraps oceanmesh.laplacian2). "
            "Boundaries are re-derived via DEM-bbox proximity. Phases "
            "D, E, F, and G are off by default."
        ),
    )
    p.add_argument("input", type=Path, help="Input fort.14.")
    p.add_argument("output", type=Path, help="Output fort.14.")
    p.add_argument(
        "--bbox", type=float, nargs=4, default=None,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help=(
            "Bounding box driving open / land classification of the "
            "surviving outer ring. Default: bbox of the input mesh's "
            "existing open-boundary nodes (preserves the original "
            "open-side classification); falls back to the full node "
            "bbox if the input has no open boundary. Pass the original "
            "DEM bbox explicitly for the most faithful reconstruction."
        ),
    )
    p.add_argument(
        "--bbox-tol-m", type=float, default=DEFAULT_BBOX_TOL_M,
        help=(
            "Tolerance in metres for 'on the bbox' open-boundary "
            f"classification. Default: {DEFAULT_BBOX_TOL_M:g} m "
            "(matches fmesh-buildmesh 0.75*hmin at hmin=200 m)."
        ),
    )
    p.add_argument(
        "--land-ibtype", type=int, default=20,
        help="ibtype written to re-derived land segments (default 20).",
    )
    p.add_argument(
        "--open-merge-coast-gap", type=int, default=0,
        help=(
            "Bridge a land run shorter than this (in nodes) sandwiched "
            "between two open runs into a single open segment. "
            "Default 0 (off)."
        ),
    )
    p.add_argument(
        "--no-remove-disjoint", dest="remove_disjoint", action="store_false",
        help="Skip Phase A (keep all dual-graph components).",
    )
    p.add_argument(
        "--min-component-elements", type=int, default=None,
        help=(
            "Phase A: keep components whose element count is >= this "
            "value. Default: keep only the single largest component."
        ),
    )
    p.add_argument(
        "--require-open-boundary", action="store_true",
        help=(
            "Phase A: keep only components that contain at least one "
            "open-boundary node."
        ),
    )
    p.add_argument(
        "--trim-dead-ends-iters", type=int, default=10,
        help=(
            "Phase B: maximum dead-end-trim iterations. 0 disables. "
            "Default 10."
        ),
    )
    p.add_argument(
        "--thin-chain-mode", choices=["widen", "delete", "none"],
        default="widen",
        help=(
            "Phase C policy for 1-cell-wide channels. 'widen' (default) "
            "inserts a centroid in every thin-chain element so each "
            "1-cell channel gains an interior node and becomes 2-cell. "
            "'delete' removes the chain entirely. 'none' skips Phase C."
        ),
    )
    p.add_argument(
        "--min-thin-chain", type=int, default=3,
        help=(
            "Phase C: minimum length of a connected thin-element run "
            "to treat as a 1-cell channel. Default 3 (matches "
            "fmesh-mesh-check)."
        ),
    )
    p.add_argument(
        "--repair-overconnected-iters", type=int, default=0,
        help=(
            "Phase D: maximum iterations of valence-balancing edge "
            "swaps. Default 0 (Phase D OFF). Set to e.g. 50 to enable. "
            "PoC #27 found this can eliminate mild over-connection "
            "(max v=9) in cleaned real meshes at a near-zero quality "
            "cost; severe gmsh-fan cases (max v >= 20) are only "
            "partially fixable by edge swap alone."
        ),
    )
    p.add_argument(
        "--max-nbr-elem", type=int, default=8,
        help=(
            "Phase D: FVCOM MAX_NBR_ELEM cap to drive every node "
            "valence to. Default 8 (matches fmesh-mesh-check). Raise "
            "if your FVCOM build is compiled with a larger cap."
        ),
    )
    p.add_argument(
        "--overconn-min-angle-floor", type=float, default=0.0,
        help=(
            "Phase D: minimum interior angle (degrees) a flip is "
            "allowed to produce. Default 0 — only triangle inversion "
            "forbidden, the value PoC #27 found practical on real "
            "meshes. Raise to 20 to forbid sliver creation; on "
            "fan-like local topology this typically rejects every "
            "candidate."
        ),
    )
    p.add_argument(
        "--under-resolved-mode",
        choices=["widen", "delete", "medial", "none"],
        default="none",
        help=(
            "Phase E policy for under-resolved channel elements "
            "(detector 6). 'widen' inserts a centroid in every flagged "
            "element so 2-cell channels become 3-cell. 'delete' removes "
            "the flagged elements. 'medial' replaces each face-face-"
            "connected channel of >= --under-resolved-min-channel-"
            "elements members with a Delaunay triangulation of (rim "
            "polygon ∪ centroid-spine sampled at h_local_median spacing) "
            "— the Stage 2 medial-axis CDT path; PoC #37 production "
            "sweet spot is min-channel-elements 10. 'none' (default) "
            "skips Phase E. Detector 6 typically flags thousands of "
            "elements on real meshes — enable deliberately."
        ),
    )
    p.add_argument(
        "--under-resolved-min-w-h", type=float, default=DEFAULT_MIN_W_H,
        help=(
            "Phase E threshold. An element is flagged when its local "
            "channel width divided by the median edge length is below "
            f"this value. Default {DEFAULT_MIN_W_H:g} (matches "
            "fmesh-mesh-check)."
        ),
    )
    p.add_argument(
        "--under-resolved-sample-ds-m", type=float,
        default=DEFAULT_CHANNEL_SAMPLE_DS_M,
        help=(
            "Phase E: boundary-sample spacing in metres for the "
            f"medial-axis detector. Default {DEFAULT_CHANNEL_SAMPLE_DS_M:g} m."
        ),
    )
    p.add_argument(
        "--under-resolved-arc-separation-factor", type=float,
        default=DEFAULT_ARC_SEPARATION_FACTOR,
        help=(
            "Phase E: arc-separation factor for same-polyline narrow "
            f"inlet detection. Default {DEFAULT_ARC_SEPARATION_FACTOR:g}."
        ),
    )
    p.add_argument(
        "--under-resolved-opposite-bank-cos-max", type=float,
        default=DEFAULT_OPPOSITE_BANK_COS_MAX,
        help=(
            "Phase E: maximum cosine of the angle between the two "
            "candidate-bank rays from the centroid. Default "
            f"{DEFAULT_OPPOSITE_BANK_COS_MAX:g} (matches fmesh-mesh-check)."
        ),
    )
    p.add_argument(
        "--under-resolved-min-channel-elements", type=int, default=1,
        help=(
            "Phase E: ignore detector-6 flags whose face-face-connected "
            "channel component has fewer than this many flagged "
            "elements. Default 1 (no filter). PoC #35 found that on "
            "real meshes most flagged clusters are tiny (~3 elements / "
            "channel); raise to e.g. 10 to limit Phase E to the long "
            "ribbon-like channels actually worth widening."
        ),
    )
    p.add_argument(
        "--repair-skewed-elements", action="store_true",
        help=(
            "Phase F switch (off by default). Delete triangles whose "
            "minimum interior angle is below "
            "--repair-skewed-min-angle-deg or whose maximum is at or "
            "above --repair-skewed-max-angle-deg. Wraps "
            "ocsmesh.utils.cleanup_skewed_el (gmsh-free)."
        ),
    )
    p.add_argument(
        "--repair-skewed-min-angle-deg", type=float,
        default=DEFAULT_SKEWED_MIN_ANGLE_DEG,
        help=(
            "Phase F: minimum interior angle a triangle is allowed to "
            f"have. Default {DEFAULT_SKEWED_MIN_ANGLE_DEG:g}° (matches "
            "ocsmesh)."
        ),
    )
    p.add_argument(
        "--repair-skewed-max-angle-deg", type=float,
        default=DEFAULT_SKEWED_MAX_ANGLE_DEG,
        help=(
            "Phase F: maximum interior angle a triangle is allowed to "
            f"have. Default {DEFAULT_SKEWED_MAX_ANGLE_DEG:g}° (matches "
            "ocsmesh)."
        ),
    )
    p.add_argument(
        "--smooth-laplacian", action="store_true",
        help=(
            "Phase G switch (off by default). Run Laplacian smoothing "
            "of all interior nodes via oceanmesh.laplacian2. Boundary "
            "nodes are auto-pinned. Connectivity, depths, and "
            "boundary lists are preserved. Note: oceanmesh is "
            "GPL-3.0-or-later; importing it propagates GPL into the "
            "redistributed combined work (see THIRD_PARTY_NOTICES.md)."
        ),
    )
    p.add_argument(
        "--smooth-laplacian-iters", type=int,
        default=DEFAULT_SMOOTH_LAPLACIAN_ITERS,
        help=(
            "Phase G: maximum smoothing iterations. Default "
            f"{DEFAULT_SMOOTH_LAPLACIAN_ITERS} (matches oceanmesh.laplacian2)."
        ),
    )
    p.add_argument(
        "--smooth-laplacian-tol", type=float,
        default=DEFAULT_SMOOTH_LAPLACIAN_TOL,
        help=(
            "Phase G: convergence tolerance on the max relative edge "
            "length change per iteration. Default "
            f"{DEFAULT_SMOOTH_LAPLACIAN_TOL:g} (matches oceanmesh)."
        ),
    )
    p.add_argument(
        "--smooth-no-repair-flipped", dest="smooth_repair_flipped",
        action="store_false",
        help=(
            "Phase G: surface raw oceanmesh.laplacian2 output without "
            "rolling back flipped triangles. Default behaviour repairs "
            "flips by reverting affected nodes' positions; pass this "
            "flag to diagnose whether a particular mesh causes "
            "flipping."
        ),
    )
    p.set_defaults(smooth_repair_flipped=True)
    p.add_argument(
        "--smooth-max-repair-passes", type=int,
        default=DEFAULT_SMOOTH_REPAIR_PASSES,
        help=(
            "Phase G: cap on the iterative-rollback loop used to "
            "repair flipped triangles. Default "
            f"{DEFAULT_SMOOTH_REPAIR_PASSES}; full rollback to the "
            "pre-smoothing positions fires if convergence is not "
            "reached within this many passes."
        ),
    )

    # Phase H: per-element greedy quality optimiser.
    p.add_argument(
        "--phase-h", action="store_true",
        help=(
            "Phase H: per-element greedy quality optimiser. Visits "
            "every element failing the per-element gate "
            "(--phase-h-alpha-target ∧ --phase-h-min-angle-target) "
            "and tries an operator inventory (Gauss-Seidel smooth, "
            "Lawson edge swap, interior + boundary edge split, "
            "vertex remove + 1-ring CDT) until the local 1-ring "
            "penalty strictly drops without flipping a triangle. "
            "Off by default — typically 10-30 min on a 47 k-element "
            "mesh and only worthwhile after Phase G has run. Phase "
            "H sits at the end of the A-G pipeline so it sees "
            "topologically-clean, Laplacian-smoothed input."
        ),
    )
    p.add_argument(
        "--phase-h-alpha-target", type=float, default=0.95,
        help=(
            "Phase H per-element alpha gate. Default 0.95; an element "
            "is 'fail' until its alpha >= this value AND its min "
            "interior angle >= --phase-h-min-angle-target."
        ),
    )
    p.add_argument(
        "--phase-h-min-angle-target", type=float, default=20.0,
        help=(
            "Phase H per-element minimum-interior-angle gate (degrees). "
            "Default 20."
        ),
    )
    p.add_argument(
        "--phase-h-max-outer-rounds", type=int, default=10,
        help=(
            "Phase H: cap on alternations of Pass A (batch smooth) "
            "and Pass B (topology operators). Default 10; convergence "
            "is typically reached in 3-5 rounds."
        ),
    )
    p.add_argument(
        "--phase-h-max-topology-per-round", type=int, default=10_000,
        help=(
            "Phase H: cap on topology accepts within a single Pass B. "
            "Default 10000."
        ),
    )
    p.add_argument(
        "--phase-h-max-smooth-sweeps", type=int, default=200,
        help=(
            "Phase H: cap on Gauss-Seidel sweeps within a single "
            "Pass A. Default 200; convergence is typically reached "
            "in 10-30 sweeps."
        ),
    )
    p.add_argument(
        "--phase-h-coastline", type=Path, action="append", default=[],
        help=(
            "Phase H (optional, repeatable): path to a coastline "
            "shapefile / GeoJSON / any GeoPandas-readable vector "
            "source. When supplied, new boundary nodes inserted by "
            "edge_split_boundary and moved by the boundary-tangent "
            "smooth are snapped onto the nearest coastline polyline "
            "within --phase-h-max-snap-m. Coordinates are auto-"
            "reprojected to EPSG:4326."
        ),
    )
    p.add_argument(
        "--phase-h-max-snap-m", type=float, default=500.0,
        help=(
            "Phase H: maximum snap distance (metres) for the "
            "coastline projector. Default 500 m (~ 2.5 × the "
            "Tokyo-Bay hmin). Proposals farther than this from any "
            "polyline fall through to the un-projected position."
        ),
    )
    p.add_argument(
        "--phase-h-lookahead", action="store_true",
        help=(
            "Phase H v4 (opt-in): enable Pass C 2-step lookahead "
            "after Pass B in every outer round. For each remaining "
            "fail element, op1 ∈ {smooth_node, vertex_remove} is "
            "applied with the penalty gate bypassed (force=True; "
            "validity unchanged), then op2 = smooth_node is searched "
            "on the elements overlapping op1's affected region. The "
            "pair is accepted iff the union penalty over the op1 ∪ "
            "op2 affected nodes strictly drops vs the round-start "
            "mesh. PoC #44 measured 61%% additional fixable on the "
            "Tokyo-Bay v3 residual (n=1000, 95%% CI ±3%%)."
        ),
    )
    p.add_argument(
        "--phase-h-max-lookahead-per-round", type=int, default=10_000,
        help=(
            "Phase H v4: cap on the number of (op1, op2) accepts in "
            "a single Pass C invocation. Default 10000. Each accept "
            "rebuilds the aux dicts on the new mesh so the cost is "
            "~250 ms / accept on a 47 k-element mesh."
        ),
    )

    p.add_argument(
        "--summary", type=Path, default=None,
        help=(
            "Optional path for the JSON summary. Default: "
            "<output_stem>_clean_summary.json next to the output fort.14."
        ),
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the before / after summary on stdout.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    if args.bbox_tol_m <= 0:
        print("--bbox-tol-m must be > 0.", file=sys.stderr)
        return 2
    if args.trim_dead_ends_iters < 0:
        print("--trim-dead-ends-iters must be >= 0.", file=sys.stderr)
        return 2
    if args.min_thin_chain < 1:
        print("--min-thin-chain must be >= 1.", file=sys.stderr)
        return 2
    if args.repair_overconnected_iters < 0:
        print("--repair-overconnected-iters must be >= 0.", file=sys.stderr)
        return 2
    if args.max_nbr_elem < 3:
        print("--max-nbr-elem must be >= 3.", file=sys.stderr)
        return 2
    if args.overconn_min_angle_floor < 0:
        print("--overconn-min-angle-floor must be >= 0.", file=sys.stderr)
        return 2
    if args.under_resolved_min_w_h <= 0:
        print("--under-resolved-min-w-h must be > 0.", file=sys.stderr)
        return 2
    if args.under_resolved_sample_ds_m <= 0:
        print("--under-resolved-sample-ds-m must be > 0.", file=sys.stderr)
        return 2
    if args.under_resolved_arc_separation_factor <= 0:
        print("--under-resolved-arc-separation-factor must be > 0.",
              file=sys.stderr)
        return 2
    if not (-1.0 <= args.under_resolved_opposite_bank_cos_max <= 1.0):
        print("--under-resolved-opposite-bank-cos-max must be in [-1, 1].",
              file=sys.stderr)
        return 2
    if args.under_resolved_min_channel_elements < 1:
        print("--under-resolved-min-channel-elements must be >= 1.",
              file=sys.stderr)
        return 2
    if args.repair_skewed_min_angle_deg < 0:
        print("--repair-skewed-min-angle-deg must be >= 0.", file=sys.stderr)
        return 2
    if args.repair_skewed_max_angle_deg > 180:
        print("--repair-skewed-max-angle-deg must be <= 180.", file=sys.stderr)
        return 2
    if args.repair_skewed_min_angle_deg >= args.repair_skewed_max_angle_deg:
        print("--repair-skewed-min-angle-deg must be < --repair-skewed-max-angle-deg.",
              file=sys.stderr)
        return 2
    if args.smooth_laplacian_iters < 1:
        print("--smooth-laplacian-iters must be >= 1.", file=sys.stderr)
        return 2
    if args.smooth_laplacian_tol <= 0:
        print("--smooth-laplacian-tol must be > 0.", file=sys.stderr)
        return 2
    if args.smooth_max_repair_passes < 0:
        print("--smooth-max-repair-passes must be >= 0.", file=sys.stderr)
        return 2
    if args.phase_h:
        if not 0.0 < args.phase_h_alpha_target <= 1.0:
            print(
                "--phase-h-alpha-target must be in (0, 1].",
                file=sys.stderr,
            )
            return 2
        if not 0.0 < args.phase_h_min_angle_target < 60.0:
            print(
                "--phase-h-min-angle-target must be in (0, 60).",
                file=sys.stderr,
            )
            return 2
        if args.phase_h_max_outer_rounds < 1:
            print(
                "--phase-h-max-outer-rounds must be >= 1.",
                file=sys.stderr,
            )
            return 2
        if args.phase_h_max_topology_per_round < 0:
            print(
                "--phase-h-max-topology-per-round must be >= 0.",
                file=sys.stderr,
            )
            return 2
        if args.phase_h_max_smooth_sweeps < 0:
            print(
                "--phase-h-max-smooth-sweeps must be >= 0.",
                file=sys.stderr,
            )
            return 2
        if args.phase_h_max_snap_m <= 0:
            print("--phase-h-max-snap-m must be > 0.", file=sys.stderr)
            return 2
        if args.phase_h_max_lookahead_per_round < 1:
            print(
                "--phase-h-max-lookahead-per-round must be >= 1.",
                file=sys.stderr,
            )
            return 2
        for coast_path in args.phase_h_coastline:
            if not coast_path.exists():
                print(
                    f"--phase-h-coastline file not found: {coast_path}",
                    file=sys.stderr,
                )
                return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(args.input)
    if args.bbox is None:
        bbox, bbox_source = _infer_bbox(mesh)
    else:
        bbox = tuple(args.bbox)
        bbox_source = "user-supplied"

    cleaned, info = clean_mesh(
        mesh,
        bbox=bbox,
        bbox_tol_m=args.bbox_tol_m,
        land_ibtype=args.land_ibtype,
        open_merge_coast_gap=args.open_merge_coast_gap,
        remove_disjoint=args.remove_disjoint,
        min_component_elements=args.min_component_elements,
        require_open_boundary=args.require_open_boundary,
        trim_dead_ends_iters=args.trim_dead_ends_iters,
        thin_chain_mode=args.thin_chain_mode,
        min_thin_chain=args.min_thin_chain,
        repair_overconnected_iters=args.repair_overconnected_iters,
        max_nbr_elem=args.max_nbr_elem,
        overconn_min_angle_floor_deg=args.overconn_min_angle_floor,
        under_resolved_mode=args.under_resolved_mode,
        under_resolved_min_w_h=args.under_resolved_min_w_h,
        under_resolved_sample_ds_m=args.under_resolved_sample_ds_m,
        under_resolved_arc_separation_factor=args.under_resolved_arc_separation_factor,
        under_resolved_opposite_bank_cos_max=args.under_resolved_opposite_bank_cos_max,
        under_resolved_min_channel_elements=args.under_resolved_min_channel_elements,
        repair_skewed=args.repair_skewed_elements,
        repair_skewed_min_angle_deg=args.repair_skewed_min_angle_deg,
        repair_skewed_max_angle_deg=args.repair_skewed_max_angle_deg,
        smooth_laplacian=args.smooth_laplacian,
        smooth_laplacian_iters=args.smooth_laplacian_iters,
        smooth_laplacian_tol=args.smooth_laplacian_tol,
        smooth_repair_flipped=args.smooth_repair_flipped,
        smooth_max_repair_passes=args.smooth_max_repair_passes,
        phase_h=args.phase_h,
        phase_h_alpha_target=args.phase_h_alpha_target,
        phase_h_min_angle_target=args.phase_h_min_angle_target,
        phase_h_max_outer_rounds=args.phase_h_max_outer_rounds,
        phase_h_max_topology_per_round=args.phase_h_max_topology_per_round,
        phase_h_max_smooth_sweeps=args.phase_h_max_smooth_sweeps,
        phase_h_coastline_paths=list(args.phase_h_coastline) or None,
        phase_h_max_snap_distance_m=args.phase_h_max_snap_m,
        phase_h_lookahead=args.phase_h_lookahead,
        phase_h_max_lookahead_per_round=args.phase_h_max_lookahead_per_round,
    )
    write_fort14(cleaned, args.output)

    summary_path = args.summary or args.output.with_name(
        args.output.stem + "_clean_summary.json"
    )
    payload = {
        "input_path": str(args.input.resolve()),
        "output_path": str(args.output.resolve()),
        "bbox_source": bbox_source,
        **info,
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if not args.quiet:
        print(f"input:  {args.input}")
        print(f"output: {args.output}")
        print(f"bbox:   {bbox}  ({bbox_source})")
        i_o = info["output"]
        i_i = info["input"]
        print(
            f"NP: {i_i['n_nodes']:,} -> {i_o['n_nodes']:,}    "
            f"NE: {i_i['n_elements']:,} -> {i_o['n_elements']:,}"
        )
        print(
            f"open: {i_i['n_open_boundaries']} -> {i_o['n_open_boundaries']}    "
            f"land: {i_i['n_land_boundaries']} -> {i_o['n_land_boundaries']}"
        )
        for ph in info["phases"]:
            print()
            print(f"phase: {ph['name']}")
            for k, v in ph.items():
                if k == "name":
                    continue
                if (
                    k == "all_component_sizes"
                    and isinstance(v, list)
                    and len(v) > 8
                ):
                    head = ", ".join(str(x) for x in v[:8])
                    print(f"  {k}: [{head}, ... ({len(v)} components)]")
                else:
                    print(f"  {k}: {v}")
        print(f"\nwrote {summary_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
