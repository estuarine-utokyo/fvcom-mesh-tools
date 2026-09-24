"""``fmesh-plot-views``: the finished mesh, whole and in close-ups.

Every solid boundary -- coastline, quay, and a wall represented as a line --
is drawn black, the open boundary red, interior edges thin grey
(``plotting.SOLID_BOUNDARY_COLOR`` etc.), so the edges water cannot cross
are recognisable in every figure.

    fmesh-plot-views outputs/refine_kimitsu_port_hires \\
        --view port:391900:394100:3908350:3910600 --view basin:392150:392850:3908650:3909350

The input is a fort.14, or a refinement output directory (its ``*.14``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import MESH_PNG_DPI, plot_mesh_views


def _view(spec: str):
    name, *vals = spec.split(":")
    if not name or len(vals) != 4:
        raise argparse.ArgumentTypeError(f"--view is NAME:x0:x1:y0:y1, got {spec!r}")
    return name, tuple(float(v) for v in vals)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="fmesh-plot-views",
                                description="Plot a mesh whole and in close-ups, "
                                            "with fixed boundary colours.")
    p.add_argument("source", type=Path, help="a fort.14, or a refinement output directory")
    p.add_argument("--view", type=_view, action="append", default=[],
                   help="a close-up, NAME:x0:x1:y0:y1 in the mesh's coordinates "
                        "(repeatable)")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="where to write (default: next to the mesh)")
    p.add_argument("--prefix", default="final_mesh", help="PNG name prefix")
    p.add_argument("--dpi", type=int, default=MESH_PNG_DPI)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    src = args.source.expanduser().resolve()
    f14 = next(src.glob("*.14"), None) if src.is_dir() else src
    if f14 is None or not f14.exists():
        print(f"fmesh-plot-views: no fort.14 at {src}", file=sys.stderr)
        return 2
    mesh = read_fort14(f14)
    out = plot_mesh_views(mesh, args.out_dir or f14.parent, dict(args.view),
                          prefix=args.prefix, dpi=args.dpi)
    for png in out:
        print(f"wrote {png}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
