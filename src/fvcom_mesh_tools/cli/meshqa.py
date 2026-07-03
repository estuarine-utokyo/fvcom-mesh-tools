"""``fmesh-mesh-qa`` CLI: the unified FVCOM acceptance gate.

One command, one pass/fail table (Japanese by default, ``--lang en``
for English), exit 0 only when every gated check passes. This is the
validator the kickoff §9 wires into ``/loop``.

Coverage (see :mod:`fvcom_mesh_tools.qa` and
``docs/fvcom_source_constraints.md``):

* FVCOM startup-fatal rules (CCW, isolated elements, R4 mixed-boundary
  elements, OBC chains, manifold topology) plus the silent hazards the
  model never checks (duplicates, orphans, tiny areas, fake ISBCE=2,
  OBC necks).
* Manual criteria C1 (min angle >= 30), C2 (max angle <= 130),
  C4 (area change <= 0.5), C5 (valence <= 8).
* OBC per-node best-edge perpendicularity and node-list ordering.
* Connectivity (single component, OBC-reachable) and min-depth clip.
* Informational: channel w/h, local Delaunay fraction, implied
  external-mode Δt (``--gate-channel-wh`` / ``--min-dt`` turn the
  first and last into gates).

Not covered yet: boundary conformity against an external coastline
polygon (needs the coastline input; planned with the recipe layer).

A JSON dump with every offender record is always written
(default ``<mesh>_qa.json``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.qa import format_report, run_qa


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-qa",
        description=(
            "Unified FVCOM mesh QA gate: FVCOM-source startup constraints, "
            "manual quality criteria C1/C2/C4/C5, open-boundary checks, "
            "connectivity, and min-depth in ONE pass/fail table. "
            "Exit 0 only when every gated check passes. Angles/areas are "
            "evaluated in a local metric projection for lon/lat meshes "
            "(FVCOM production builds are CARTESIAN)."
        ),
    )
    p.add_argument("input", type=Path, help="fort.14 mesh to check.")

    g = p.add_argument_group("gate thresholds")
    g.add_argument("--min-angle", type=float, default=30.0,
                   help="C1: minimum interior angle in degrees (default 30).")
    g.add_argument("--max-angle", type=float, default=130.0,
                   help="C2: maximum interior angle in degrees (default 130).")
    g.add_argument("--max-area-change", type=float, default=0.5,
                   help="C4: adjacent-element area change (default 0.5).")
    g.add_argument("--max-valence", type=int, default=8,
                   help="C5: elements per node (default 8).")
    g.add_argument("--min-depth", type=float, default=2.0,
                   help="Minimum-depth clip to verify, metres (default 2.0).")
    g.add_argument("--max-obc-perp-dev", type=float, default=20.0,
                   help="Max best-edge deviation from perpendicular per OBC "
                        "node, degrees (default 20).")
    g.add_argument("--gate-channel-wh", type=float, default=None,
                   help="Turn the advisory channel w/h check into a gate at "
                        "this ratio (e.g. 1.0 = at least one cell across).")
    g.add_argument("--min-dt", type=float, default=None,
                   help="Turn the advisory implied-dt check into a gate: fail "
                        "if any element implies an external dt below this (s).")

    p.add_argument("--coords", choices=("auto", "lonlat", "metric"),
                   default="auto",
                   help="Node-coordinate interpretation (default auto-detect).")
    p.add_argument("--no-channel", action="store_true",
                   help="Skip the (slow) channel-width metric entirely.")
    p.add_argument("--lang", choices=("ja", "en"), default="ja",
                   help="Report language (default ja, per the project spec).")
    p.add_argument("--max-offenders", type=int, default=5,
                   help="Worst offenders listed per failed check (default 5).")
    p.add_argument("--json", type=Path, default=None,
                   help="JSON report path (default <input>_qa.json).")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress the stdout table (JSON is still written).")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2

    mesh = read_fort14(args.input)
    report = run_qa(
        mesh,
        name=args.input.name,
        path=args.input.resolve(),
        min_angle_deg=args.min_angle,
        max_angle_deg=args.max_angle,
        max_area_change=args.max_area_change,
        max_valence=args.max_valence,
        min_depth_m=args.min_depth,
        max_obc_perp_dev_deg=args.max_obc_perp_dev,
        coords=args.coords,
        channel_check=not args.no_channel,
        min_channel_wh_gate=args.gate_channel_wh,
        min_dt_s=args.min_dt,
        max_offenders=args.max_offenders,
    )

    if not args.quiet:
        print(format_report(report, lang=args.lang))

    json_path = args.json or args.input.with_name(args.input.stem + "_qa.json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    if not args.quiet:
        print(f"\nwrote {json_path}")

    return 0 if report.passed else 1


if __name__ == "__main__":
    sys.exit(main())
