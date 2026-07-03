"""``fmesh-plot-mesh`` CLI: kickoff §10 figures.

Whole-domain figure (xcoast land + labeled 5 km reference grid +
open-boundary highlight) plus optional per-cell zoom panels addressed
by grid references ("C4", "C4-D5") or named aliases.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import plot_mesh_overview


def _parse_alias(spec: str) -> tuple[str, tuple[float, float, float, float]]:
    name, _, rest = spec.partition("=")
    vals = [float(v) for v in rest.split(",")]
    if not name or len(vals) != 4:
        raise ValueError(f"alias must be NAME=xmin,ymin,xmax,ymax — got {spec!r}")
    return name, (vals[0], vals[1], vals[2], vals[3])


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-plot-mesh",
        description=(
            "Render the whole-domain mesh figure (land + labeled "
            "reference grid + open boundary) and optional zoom panels."
        ),
    )
    p.add_argument("input", type=Path, help="fort.14 mesh to plot.")
    p.add_argument("--out", type=Path, default=None,
                   help="Overview PNG path (default <input>_overview.png).")
    p.add_argument("--crs", default="EPSG:32654",
                   help="Projected CRS of the mesh nodes (default EPSG:32654).")
    p.add_argument("--cell-km", type=float, default=5.0,
                   help="Reference-grid cell size in km (default 5).")
    p.add_argument("--no-grid", action="store_true",
                   help="Disable the reference grid overlay.")
    p.add_argument("--coast-preset", default=None,
                   help="xcoast preset name for land rendering "
                        "(e.g. tokyo_bay).")
    p.add_argument("--coast-bbox", type=float, nargs=4, default=None,
                   metavar=("W", "S", "E", "N"),
                   help="xcoast lon/lat bbox for land rendering.")
    p.add_argument("--zoom", action="append", default=[],
                   metavar="REF",
                   help="Zoom panel reference (cell 'C4', range 'C4-D5', "
                        "or alias). Repeatable.")
    p.add_argument("--alias", action="append", default=[],
                   metavar="NAME=XMIN,YMIN,XMAX,YMAX",
                   help="Named region in CRS coordinates. Repeatable.")
    p.add_argument("--dpi", type=int, default=200,
                   help="Raster DPI (default 200; mesh detail needs 400+).")
    return p


def main(argv: list[str] | None = None) -> int:
    os.environ.setdefault("MPLBACKEND", "Agg")
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    try:
        aliases = dict(_parse_alias(s) for s in args.alias)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 2

    coast = None
    if args.coast_preset and args.coast_bbox:
        print("use either --coast-preset or --coast-bbox, not both",
              file=sys.stderr)
        return 2
    if args.coast_preset:
        coast = args.coast_preset
    elif args.coast_bbox:
        coast = tuple(args.coast_bbox)

    mesh = read_fort14(args.input)
    cell_m = None if args.no_grid else args.cell_km * 1000.0
    out = args.out or args.input.with_name(args.input.stem + "_overview.png")

    written = [plot_mesh_overview(
        mesh, out, crs=args.crs, cell_m=cell_m, aliases=aliases,
        coast=coast, dpi=args.dpi,
    )]
    for ref in args.zoom:
        safe = ref.replace("–", "-").replace("=", "_")
        zoom_png = out.with_name(f"{args.input.stem}_{safe}.png")
        try:
            written.append(plot_mesh_overview(
                mesh, zoom_png, crs=args.crs, cell_m=cell_m,
                aliases=aliases, coast=coast, zoom=ref, dpi=args.dpi,
            ))
        except ValueError as e:
            print(f"zoom {ref!r} skipped: {e}", file=sys.stderr)

    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
