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

from fvcom_mesh_tools.diagnostics import (
    DEFAULT_MAX_NBR_ELEM,
    DEFAULT_MIN_THIN_CHAIN,
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
