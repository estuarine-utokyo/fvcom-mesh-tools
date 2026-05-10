"""``fmesh-mesh-quality`` CLI: unified quality metrics + threshold gate.

Three modes, picked from the number of input files supplied:

    * **Single mesh**: print one column of the standard metrics
      (``alpha_mean``, ``frac<20°``, ``max_valence``, etc.) and write
      a JSON dump.
    * **Two meshes (before / after)**: print a side-by-side table
      with a ``delta`` column. Useful as a one-shot
      ``fmesh-mesh-clean`` validation harness.
    * **Three or more meshes**: print all columns side-by-side. The
      first label is treated as "before" only for the threshold gate.

The threshold gate (``--min-alpha`` / ``--max-frac-lt-20deg`` /
``--max-valence`` / ...) is evaluated against the **last** mesh
supplied (i.e. "after" in a before/after pair). Any failure exits 1
so the CLI is usable as a CI gate, paralleling ``fmesh-mesh-check``.

JSON output (``--summary PATH``, default ``<last>_quality.json``):

    {
      "meshes": [
        {"label": "before", "path": ".../tokyo.14", "metrics": {...}},
        {"label": "after",  "path": ".../tokyo_clean.14", "metrics": {...}}
      ],
      "thresholds": {"min_alpha_mean": 0.95, ...},
      "checks": [{"metric": "alpha_mean", "op": "≥", ...}, ...],
      "passed": true
    }
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from fvcom_mesh_tools.diagnostics import DEFAULT_MAX_NBR_ELEM
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.quality import (
    check_thresholds,
    compute_metrics,
    format_comparison_table,
    format_threshold_table,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-quality",
        description=(
            "Compute unified quality metrics for one or more FVCOM "
            "fort.14 meshes. Single mesh prints metrics; two meshes "
            "compare with a delta column; thresholds turn the command "
            "into a CI gate (exit 1 on any threshold failure, "
            "evaluated against the last mesh)."
        ),
    )
    p.add_argument(
        "inputs", type=Path, nargs="+",
        help="One or more fort.14 files. Order matters: thresholds "
             "are evaluated against the LAST mesh.",
    )
    p.add_argument(
        "--labels", nargs="+", default=None,
        help="Display labels (one per input). Default: file stem.",
    )
    p.add_argument(
        "--max-nbr-elem", type=int, default=DEFAULT_MAX_NBR_ELEM,
        help=(
            "MAX_NBR_ELEM cap used to count over-connected nodes. "
            f"Default {DEFAULT_MAX_NBR_ELEM} (matches fmesh-mesh-check)."
        ),
    )

    g = p.add_argument_group("threshold gates (apply to the last input)")
    g.add_argument("--min-alpha", type=float, default=None,
                   help="Fail if alpha_mean < this (e.g. 0.95).")
    g.add_argument("--max-frac-lt-20deg", type=float, default=None,
                   help="Fail if frac_lt_20deg > this (fraction in [0, 1]).")
    g.add_argument("--max-valence", type=int, default=None,
                   help="Fail if max_valence > this (typical 8).")
    g.add_argument("--max-overconnected", type=int, default=None,
                   help="Fail if the number of over-connected nodes > this.")
    g.add_argument("--max-flipped", type=int, default=None,
                   help="Fail if n_flipped > this. Default usage: 0.")
    g.add_argument("--max-disjoint-elems", type=int, default=None,
                   help="Fail if NE - NE_largest_component > this.")

    p.add_argument(
        "--summary", type=Path, default=None,
        help="Path for the JSON summary. Default: "
             "<last>_quality.json next to the last input.",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the comparison / threshold tables on stdout.",
    )
    return p


def _resolve_labels(
    inputs: list[Path], labels: list[str] | None,
) -> list[str]:
    if labels is None:
        return [p.stem for p in inputs]
    if len(labels) != len(inputs):
        raise ValueError(
            f"--labels count ({len(labels)}) must match input count ({len(inputs)})"
        )
    return list(labels)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    for path in args.inputs:
        if not path.exists():
            print(f"input not found: {path}", file=sys.stderr)
            return 2

    if args.max_nbr_elem < 3:
        print("--max-nbr-elem must be >= 3.", file=sys.stderr)
        return 2

    try:
        labels = _resolve_labels(args.inputs, args.labels)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 2

    rows: list[tuple[str, dict]] = []
    for label, path in zip(labels, args.inputs):
        mesh = read_fort14(path)
        metrics = compute_metrics(mesh, max_nbr_elem=args.max_nbr_elem)
        rows.append((label, metrics))

    if not args.quiet:
        print(format_comparison_table(rows))

    last_metrics = rows[-1][1]
    thresholds = {
        k: v for k, v in {
            "min_alpha_mean": args.min_alpha,
            "max_frac_lt_20deg": args.max_frac_lt_20deg,
            "max_valence": args.max_valence,
            "max_overconnected": args.max_overconnected,
            "max_flipped": args.max_flipped,
            "max_disjoint_elems": args.max_disjoint_elems,
        }.items() if v is not None
    }
    passed, checks = check_thresholds(last_metrics, **thresholds)
    if checks and not args.quiet:
        print()
        print(f"thresholds (against '{rows[-1][0]}'):")
        print(format_threshold_table(checks))
        print()
        print(f"overall: {'PASS' if passed else 'FAIL'}")

    summary_path = args.summary or args.inputs[-1].with_name(
        args.inputs[-1].stem + "_quality.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meshes": [
            {
                "label": label,
                "path": str(path.resolve()),
                "metrics": metrics,
            }
            for (label, metrics), path in zip(rows, args.inputs)
        ],
        "thresholds": thresholds,
        "checks": [c.to_dict() for c in checks],
        "passed": bool(passed),
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if not args.quiet:
        print(f"\nwrote {summary_path}")

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
