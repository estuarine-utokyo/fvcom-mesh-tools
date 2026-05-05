"""``fmesh-subset-dem`` CLI: clip a DEM raster to a lon/lat bbox and
write a CF-tagged GeoTIFF that ``fmesh-buildmesh`` can ingest.

Two source families are recognised automatically:

* GeoTIFF / CF-NetCDF with an embedded CRS (rasterio path), or
* lon/lat NetCDF without CRS metadata - SRTM15+, GEBCO 2024, ... -
  selected by passing ``--src-var <z|elevation|...>``.

Example::

    fmesh-subset-dem SRTM15+.nc osaka.tif \\
        --bbox 134.90 34.20 135.55 34.85 \\
        --src-var z
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fvcom_mesh_tools.io import subset_dem_to_geotiff


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-subset-dem",
        description=(
            "Clip a global DEM to a lon/lat bounding box and write a "
            "CF-tagged GeoTIFF for downstream meshing."
        ),
    )
    p.add_argument("src", type=Path, help="Source DEM (GeoTIFF / NetCDF).")
    p.add_argument("dst", type=Path, help="Output GeoTIFF path.")
    p.add_argument(
        "--bbox", type=float, nargs=4,
        metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"),
        required=True,
        help="Lon/lat bounding box (degrees).",
    )
    p.add_argument(
        "--src-var", type=str, default=None, metavar="NAME",
        help=(
            "Force the lon/lat-NetCDF read path with this variable "
            "name. Required for SRTM15+ / GEBCO style inputs that lack "
            "an embedded CRS."
        ),
    )
    p.add_argument(
        "--src-crs", type=str, default="EPSG:4326",
        help=(
            "CRS to tag the output with on the lon/lat-NetCDF path "
            "(default: EPSG:4326). Ignored on the rasterio path."
        ),
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.src.exists():
        print(f"source not found: {args.src}", file=sys.stderr)
        return 2
    minlon, minlat, maxlon, maxlat = args.bbox
    if not (minlon < maxlon and minlat < maxlat):
        print(
            "--bbox must satisfy minlon < maxlon and minlat < maxlat.",
            file=sys.stderr,
        )
        return 2

    info = subset_dem_to_geotiff(
        args.src,
        args.dst,
        (minlon, minlat, maxlon, maxlat),
        src_var=args.src_var,
        src_crs=args.src_crs,
    )

    print(f"src:    {args.src}")
    print(f"dst:    {args.dst}")
    print(f"path:   {info['path_dependent']}")
    print(f"shape:  {info['shape'][0]} rows x {info['shape'][1]} cols")
    print(f"crs:    {info['crs']}")
    print(
        f"bbox:   ({info['bbox'][0]:.6f}, {info['bbox'][1]:.6f}, "
        f"{info['bbox'][2]:.6f}, {info['bbox'][3]:.6f})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
