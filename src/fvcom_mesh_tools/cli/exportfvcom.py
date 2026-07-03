"""``fmesh-export-fvcom`` CLI: fort.14 -> FVCOM native input set.

Writes ``<case>_grd.dat`` / ``<case>_dep.dat`` / ``<case>_obc.dat``
(+ optional ``_cor.dat`` / ``_spg.dat`` hooks and an SMS ``.2dm``)
per the formats in ``docs/fvcom_source_constraints.md``.

The Coriolis column is FVCOM's per-node latitude (CARTESIAN builds):
``--cor y`` takes it from the node y coordinate (lon/lat meshes),
``--cor crs --crs EPSG:32654`` inverse-projects projected nodes back
to latitude, ``--cor none`` (default) skips the file.

Writers refuse flipped/zero-area elements and unreferenced nodes;
run ``fmesh-mesh-qa`` first — export is meant for gate-passing meshes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-export-fvcom",
        description=(
            "Export a fort.14 mesh to FVCOM's native ASCII inputs "
            "(_grd/_dep/_obc[.../_cor/_spg] .dat + SMS .2dm)."
        ),
    )
    p.add_argument("input", type=Path, help="fort.14 mesh to export.")
    p.add_argument("--casename", default=None,
                   help="Case name prefix (default: input file stem).")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Output directory (default: alongside the input).")
    p.add_argument("--obc-type", type=int, default=1,
                   help="FVCOM OBC type 1-10 for every open node "
                        "(default 1 = tidal elevation, linear).")
    p.add_argument("--cor", choices=("none", "y", "crs"), default="none",
                   help="Coriolis latitude source: none (skip _cor.dat), "
                        "'y' (node y is latitude), or 'crs' "
                        "(inverse-project via --crs).")
    p.add_argument("--crs", default=None,
                   help="Projected CRS of the mesh (e.g. EPSG:32654); "
                        "required with --cor crs.")
    p.add_argument("--write-empty-spg", action="store_true",
                   help="Write a valid zero-node _spg.dat hook.")
    p.add_argument("--no-2dm", action="store_true",
                   help="Skip the SMS .2dm export.")
    p.add_argument("--z-convention", choices=("depth", "elevation"),
                   default="depth",
                   help="ND z column of the .2dm: fort.14 depth as-is "
                        "(default) or -depth.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    if args.cor == "crs" and not args.crs:
        print("--cor crs requires --crs (e.g. --crs EPSG:32654)", file=sys.stderr)
        return 2

    mesh = read_fort14(args.input)
    casename = args.casename or args.input.stem
    outdir = args.outdir or args.input.parent

    cor = None
    if args.cor == "y":
        cor = mesh.nodes[:, 1]
    elif args.cor == "crs":
        from pyproj import Transformer

        tr = Transformer.from_crs(args.crs, "EPSG:4326", always_xy=True)
        _lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        cor = lat

    try:
        written = export_fvcom_case(
            mesh, outdir, casename,
            obc_type=args.obc_type,
            cor=cor,
            write_empty_spg=args.write_empty_spg,
            twodm=not args.no_2dm,
            z_convention=args.z_convention,
        )
    except ValueError as e:
        print(f"export refused: {e}", file=sys.stderr)
        return 1

    for kind, path in written.items():
        print(f"{kind}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
