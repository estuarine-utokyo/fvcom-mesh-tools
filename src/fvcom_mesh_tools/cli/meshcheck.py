"""``fmesh-mesh-check`` CLI: detect inadequate FVCOM meshes.

Loads a fort.14 file, runs all six detectors in
:mod:`fvcom_mesh_tools.diagnostics`, and writes:

    * ``<prefix>_summary.txt`` — single-mesh summary table.
    * ``<prefix>_diag.json``   — structured per-element / per-node flag
      records with coordinates, ready for downstream repair pipelines.
    * ``<prefix>_map.png``     — mesh with flagged elements / nodes
      overlaid (skip with ``--no-plot``).

Exit code is 0 if every detector returns zero flags and 1 otherwise, so
the command is usable as a CI gate.
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
    DEFAULT_MAX_NBR_ELEM,
    DEFAULT_MIN_THIN_CHAIN,
    DEFAULT_MIN_W_H,
    DEFAULT_OPPOSITE_BANK_COS_MAX,
    plot_report,
    report_to_dict,
    report_to_summary_text,
    run_diagnostics,
)
from fvcom_mesh_tools.io import read_fort14


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-check",
        description=(
            "Detect inadequate FVCOM meshes. Six detectors flag elements "
            "and nodes likely to cause problems for FVCOM, especially in "
            "narrow water bodies (rivers, canals, harbours): disjoint "
            "wet-domain components, dead-end elements, thin elements, "
            "thin chains (1-cell channels), over-connected nodes, and "
            "open-boundary unreachable elements. No mesh repair is "
            "performed."
        ),
    )
    p.add_argument("input", type=Path, help="Input fort.14 file.")
    p.add_argument(
        "--out-prefix", type=Path, default=None,
        help=(
            "Output file prefix. Defaults to the input file's path with "
            "the '.14' suffix stripped, so 'foo.14' produces "
            "'foo_summary.txt', 'foo_diag.json', 'foo_map.png'."
        ),
    )
    p.add_argument(
        "--max-nbr-elem", type=int, default=DEFAULT_MAX_NBR_ELEM,
        help=(
            "Per-node element-neighbour cap. Nodes with valence > this "
            "value are flagged. The conservative legacy FVCOM value is "
            f"{DEFAULT_MAX_NBR_ELEM}; recent 4.x builds raise it. Set to "
            "match the cap of your FVCOM build."
        ),
    )
    p.add_argument(
        "--min-thin-chain", type=int, default=DEFAULT_MIN_THIN_CHAIN,
        help=(
            "Minimum length of a connected run of thin elements that we "
            f"flag as a 1-cell-wide channel. Default {DEFAULT_MIN_THIN_CHAIN}."
        ),
    )
    p.add_argument(
        "--min-w-h", type=float, default=DEFAULT_MIN_W_H,
        help=(
            "Minimum cells across a channel for the medial-axis "
            "detector. An element is flagged when its local channel "
            "width divided by the median edge length is below this "
            f"value. Default {DEFAULT_MIN_W_H:g}."
        ),
    )
    p.add_argument(
        "--channel-sample-ds-m", type=float,
        default=DEFAULT_CHANNEL_SAMPLE_DS_M,
        help=(
            "Boundary-sample spacing in metres used by the channel-width "
            f"detector. Default {DEFAULT_CHANNEL_SAMPLE_DS_M:g} m."
        ),
    )
    p.add_argument(
        "--channel-arc-separation-factor", type=float,
        default=DEFAULT_ARC_SEPARATION_FACTOR,
        help=(
            "Two boundary samples on the same polyline are treated as "
            "different banks only if their along-polyline arc "
            "separation exceeds this factor times the distance from "
            "the query point to the nearest sample. Default "
            f"{DEFAULT_ARC_SEPARATION_FACTOR:g}."
        ),
    )
    p.add_argument(
        "--channel-opposite-bank-cos-max", type=float,
        default=DEFAULT_OPPOSITE_BANK_COS_MAX,
        help=(
            "Maximum cosine of the angle between (centroid -> nearest "
            "sample) and (centroid -> far-arc sample) on the same "
            "polyline; only accept the far-arc sample as the 'other "
            "bank' when this cosine is below the threshold. Default "
            f"{DEFAULT_OPPOSITE_BANK_COS_MAX:g} (angle > "
            f"{180.0 / 3.14159 * np.arccos(DEFAULT_OPPOSITE_BANK_COS_MAX):.0f}°)."
        ),
    )
    p.add_argument(
        "--no-plot", action="store_true",
        help="Skip writing the *_map.png overlay.",
    )
    p.add_argument(
        "--no-json", action="store_true",
        help="Skip writing the *_diag.json structured report.",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the summary table on stdout (still writes the .txt).",
    )
    return p


def _resolve_prefix(input_path: Path, prefix: Path | None) -> Path:
    if prefix is not None:
        prefix.parent.mkdir(parents=True, exist_ok=True)
        return prefix
    if input_path.suffix == ".14":
        return input_path.with_suffix("")
    return input_path.with_name(input_path.name + "_check")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    if args.max_nbr_elem < 3:
        print("--max-nbr-elem must be >= 3 (no Delaunay mesh has fewer).",
              file=sys.stderr)
        return 2
    if args.min_thin_chain < 1:
        print("--min-thin-chain must be >= 1.", file=sys.stderr)
        return 2
    if args.min_w_h <= 0:
        print("--min-w-h must be > 0.", file=sys.stderr)
        return 2
    if args.channel_sample_ds_m <= 0:
        print("--channel-sample-ds-m must be > 0.", file=sys.stderr)
        return 2
    if args.channel_arc_separation_factor <= 0:
        print("--channel-arc-separation-factor must be > 0.", file=sys.stderr)
        return 2
    if not (-1.0 <= args.channel_opposite_bank_cos_max <= 1.0):
        print("--channel-opposite-bank-cos-max must be in [-1, 1].", file=sys.stderr)
        return 2

    prefix = _resolve_prefix(args.input, args.out_prefix)
    summary_path = prefix.with_name(prefix.name + "_summary.txt")
    json_path = prefix.with_name(prefix.name + "_diag.json")
    map_path = prefix.with_name(prefix.name + "_map.png")

    mesh = read_fort14(args.input)
    report = run_diagnostics(
        mesh,
        name=args.input.stem,
        path=args.input.resolve(),
        max_nbr_elem=args.max_nbr_elem,
        min_thin_chain=args.min_thin_chain,
        min_w_h=args.min_w_h,
        channel_sample_ds_m=args.channel_sample_ds_m,
        channel_arc_separation_factor=args.channel_arc_separation_factor,
        channel_opposite_bank_cos_max=args.channel_opposite_bank_cos_max,
    )

    summary = report_to_summary_text(report)
    summary_path.write_text(summary + "\n", encoding="utf-8")
    if not args.quiet:
        print(summary)
        print(f"\nwrote {summary_path}")

    if not args.no_json:
        json_path.write_text(
            json.dumps(report_to_dict(report), indent=2), encoding="utf-8",
        )
        if not args.quiet:
            print(f"wrote {json_path}")

    if not args.no_plot:
        plot_report(report, map_path)
        if not args.quiet:
            print(f"wrote {map_path}")

    return 0 if not report.any_flagged() else 1


if __name__ == "__main__":
    sys.exit(main())
