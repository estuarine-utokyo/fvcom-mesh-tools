"""fmesh-prep-shoreline: engineered-shoreline preprocessing CLI.

One command turns a bbox into the two v5 front-end artifacts:

* ``land_opened.shp``   — LAND-opened true-land polygons (thin
  artificial structures erased, water connectivity intact);
* ``skeleton_seeds.shp`` — water medial-axis seed lines for
  narrow-but-essential water (optional).

Every parameter of the PoC #90/#91 scripts is exposed; a
``provenance.json`` records them all for reproducibility.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="fmesh-prep-shoreline",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--bbox", type=float, nargs=4, required=True,
        metavar=("LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX"),
        help="Fetch/processing extent in EPSG:4326.",
    )
    p.add_argument(
        "--out-dir", type=Path, required=True,
        help="Output directory (created if missing).",
    )
    p.add_argument(
        "--r-open", type=float, default=150.0, metavar="METRES",
        help="Land-opening radius: land thinner than 2r is erased "
             "(default 150 -> erases <300 m structures; match hmin/2).",
    )
    p.add_argument(
        "--min-island", type=float, default=3.6e5, metavar="M2",
        help="Drop islands smaller than this after opening "
             "(default 0.36 km^2 = 4*(300 m)^2).",
    )
    p.add_argument(
        "--min-water-area", type=float, default=1e-5, metavar="DEG2",
        help="xcoast: minimum water-polygon area subtracted from land.",
    )
    p.add_argument(
        "--land-shp", type=Path, default=None,
        help="OSM land-polygons source shapefile (needed on a cache "
             "miss; default: $DATA_DIR/OSM/land-polygons-split-4326/"
             "land_polygons.shp when DATA_DIR is set).",
    )
    p.add_argument(
        "--cache-dir", type=Path, default=None,
        help="xcoast download cache (default: xcoast's own).",
    )
    p.add_argument(
        "--force-fetch", action="store_true",
        help="Re-download/rebuild the xcoast mask even if cached.",
    )
    p.add_argument(
        "--skeleton", action=argparse.BooleanOptionalAction, default=True,
        help="Also extract water medial-axis seed lines.",
    )
    p.add_argument(
        "--half-width", type=float, nargs=2, default=(140.0, 460.0),
        metavar=("LO", "HI"),
        help="Skeleton kept where local half-width is in [LO, HI] m "
             "(default 140 460 -> widths 280-920 m at hmin 300).",
    )
    p.add_argument(
        "--px", type=float, default=50.0, metavar="METRES",
        help="Skeleton raster cell size.",
    )
    p.add_argument(
        "--min-water-area", type=float, default=1e-5,
        help="min_water_area_deg2 for the xcoast true-land product. "
             "1e-5 (~1e5 m2) drops fragmented riverbank polygons and "
             "leaves rivers as LAND (Tamagawa/Arakawa gap, PoC #114); "
             "1e-7 opens them.",
    )
    p.add_argument(
        "--obc-line", type=float, nargs="+", default=None,
        metavar="LONLAT",
        help="Artificial open-boundary line as lon lat pairs "
             "(southern end first). Water seaward of it is walled "
             "off so the domain ends at the line.",
    )
    p.add_argument(
        "--obc-perpendicular-ends",
        action=argparse.BooleanOptionalAction, default=True,
        help="Replace each OBC-line end with a straight segment "
             "along the local coast normal (perpendicular junction; "
             "storm-surge stability practice).",
    )
    p.add_argument(
        "--utm-epsg", type=int, default=None,
        help="Metric CRS override (default: auto UTM zone).",
    )
    args = p.parse_args(argv)

    from fvcom_mesh_tools.prep import (
        fetch_true_land,
        open_land,
        water_skeleton_lines,
    )

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[prep] fetching OSM true-land for bbox {args.bbox} ...",
          flush=True)
    land = fetch_true_land(
        tuple(args.bbox),
        land_shp_path=args.land_shp,
        min_water_area_deg2=args.min_water_area,
        cache_dir=args.cache_dir,
        force=args.force_fetch,
    )
    print(f"[prep] true-land polygons: {len(land)}", flush=True)

    opened = open_land(
        land,
        r_open_m=args.r_open,
        min_island_area_m2=args.min_island,
        clip_bbox=tuple(args.bbox),
        utm_epsg=args.utm_epsg,
    )
    obc_pts = None
    if args.obc_line:
        from fvcom_mesh_tools.prep.shoreline import (
            cut_domain_at_obc_line,
            extend_obc_ends_perpendicular,
        )

        obc_pts = [(args.obc_line[i], args.obc_line[i + 1])
                   for i in range(0, len(args.obc_line), 2)]
        if args.obc_perpendicular_ends:
            obc_pts = extend_obc_ends_perpendicular(
                obc_pts, opened, utm_epsg=args.utm_epsg,
            )
            print("[prep] OBC ends extended perpendicular to the "
                  "coast", flush=True)
        opened = cut_domain_at_obc_line(opened, obc_pts,
                                        tuple(args.bbox))
        print(f"[prep] domain cut at OBC line ({len(obc_pts)} pts)",
              flush=True)
    land_path = out_dir / "land_opened.shp"
    opened.to_file(land_path)
    print(f"[prep] wrote {land_path} ({len(opened)} polygons, "
          f"r_open={args.r_open:g} m)", flush=True)

    skel_path = None
    if args.skeleton:
        seeds = water_skeleton_lines(
            opened,
            px_m=args.px,
            half_width_range_m=tuple(args.half_width),
            utm_epsg=args.utm_epsg,
        )
        skel_path = out_dir / "skeleton_seeds.shp"
        seeds.to_file(skel_path)
        print(f"[prep] wrote {skel_path} ({len(seeds)} lines, "
              f"half-width {args.half_width[0]:g}-"
              f"{args.half_width[1]:g} m)", flush=True)

    (out_dir / "provenance.json").write_text(json.dumps({
        "tool": "fmesh-prep-shoreline",
        "bbox": list(args.bbox),
        "r_open_m": args.r_open,
        "min_island_area_m2": args.min_island,
        "min_water_area_deg2": args.min_water_area,
        "skeleton": bool(args.skeleton),
        "half_width_m": list(args.half_width),
        "px_m": args.px,
        "obc_line": args.obc_line,
        "obc_line_effective": obc_pts,
        "obc_perpendicular_ends": bool(args.obc_perpendicular_ends),
        "utm_epsg": args.utm_epsg,
        "outputs": {
            "land_opened": str(land_path),
            "skeleton_seeds": str(skel_path) if skel_path else None,
        },
    }, indent=2), encoding="utf-8")
    print(f"[prep] provenance -> {out_dir / 'provenance.json'}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
