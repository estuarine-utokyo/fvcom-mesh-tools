"""``fmesh-mesh-clean`` CLI: prune disjoint pools and trim dead-end elements.

Phase A removes dual-graph connected components by size and / or
open-boundary touch (default keeps only the largest component). Phase B
iteratively deletes degree-1 elements with no open-boundary edge so
"spit" terminations of 1-cell channels are removed. Boundaries are
re-derived against a DEM-bbox classifier matching ``fmesh-buildmesh``.

This command does **not** repair thin elements, thin chains, or
over-connected nodes; those stages are deferred to follow-up work.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import DEFAULT_BBOX_TOL_M, clean_mesh


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
            "Clean an FVCOM mesh by pruning disjoint wet pools and "
            "trimming dead-end elements. Phase A keeps dual-graph "
            "components by size and / or open-boundary touch; Phase B "
            "iteratively trims degree-1 elements with no open-boundary "
            "edge. Boundaries are re-derived via DEM-bbox proximity, "
            "matching the fmesh-buildmesh convention. Thin / "
            "1-cell-channel / over-connected-node repair is NOT "
            "performed here."
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
