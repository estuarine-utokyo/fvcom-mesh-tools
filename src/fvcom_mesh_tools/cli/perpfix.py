"""``fmesh-perpfix`` CLI: align first-ring interior edges perpendicular to
the open-boundary tangent.

Reads a fort.14 file, runs
:func:`fvcom_mesh_tools.algorithms.align_open_boundary_first_ring`, and
writes the modified mesh out as another fort.14. Prints before/after
metrics and signed-area sanity-check to stdout.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14


def _stats(name: str, x: np.ndarray) -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  min={np.nanmin(x):.4f}  "
        f"p50={np.nanpercentile(x, 50):.4f}  p95={np.nanpercentile(x, 95):.4f}  "
        f"max={np.nanmax(x):.4f}  mean={np.nanmean(x):.4f}"
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-perpfix",
        description=(
            "Move the first-ring interior neighbours of every open-boundary "
            "node so each incident edge is as perpendicular to the local "
            "boundary tangent as possible while preserving original edge "
            "length. Boundary nodes (open + land) are kept fixed."
        ),
    )
    p.add_argument("input", type=Path, help="Input fort.14 file.")
    p.add_argument("output", type=Path, help="Output fort.14 file.")
    p.add_argument(
        "--alpha", type=float, default=1.0,
        help="Per-iteration damping factor in (0, 1] (default: 1.0).",
    )
    p.add_argument(
        "--iters", type=int, default=1, dest="n_iters",
        help="Number of perpendicular-projection iterations (default: 1).",
    )
    p.add_argument(
        "--smooth-iters", type=int, default=0,
        help="Optional Laplacian smoothing passes on second-ring nodes.",
    )
    p.add_argument(
        "--smooth-alpha", type=float, default=0.3,
        help="Damping factor for the Laplacian smoothing pass (default: 0.3).",
    )
    p.add_argument(
        "--segment-index", type=int, default=0,
        help="Open-boundary segment to operate on (default: 0).",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the before/after summary; just write the file.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    if not 0.0 < args.alpha <= 1.0:
        print("--alpha must lie in (0, 1].", file=sys.stderr)
        return 2
    if args.n_iters < 1:
        print("--iters must be >= 1.", file=sys.stderr)
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)

    before = read_fort14(args.input)
    perp_before = open_bdy_perpendicularity(before, segment_index=args.segment_index)
    n_flipped_before = int((signed_areas(before) <= 0).sum())

    after, info = align_open_boundary_first_ring(
        before,
        alpha=args.alpha,
        n_iters=args.n_iters,
        smooth_iters=args.smooth_iters,
        smooth_alpha=args.smooth_alpha,
        segment_index=args.segment_index,
    )
    perp_after = open_bdy_perpendicularity(after, segment_index=args.segment_index)
    n_flipped_after = int((signed_areas(after) <= 0).sum())

    write_fort14(after, args.output)

    if not args.quiet:
        print(f"input:                       {args.input}")
        print(f"output:                      {args.output}")
        print(f"NP={before.n_nodes:,}   NE={before.n_elements:,}")
        print(
            f"alpha={args.alpha}  iters={args.n_iters}  "
            f"smooth_iters={args.smooth_iters}  segment_index={args.segment_index}"
        )
        print(f"interior nodes moved:        {info['moved']:,}")
        print(
            f"movable first-ring by parent count: "
            f"{info['movable_first_ring_by_parent_count']}"
        )
        print()
        print("[FVCOM open-boundary perpendicularity (deg from 90)]")
        print(_stats("  before", perp_before))
        print(_stats("  after ", perp_after))
        print()
        print("[Element validity]")
        print(
            f"  flipped triangles before / after: "
            f"{n_flipped_before:,} / {n_flipped_after:,}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
