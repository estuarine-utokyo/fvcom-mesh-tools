# Did a change to the tool move a coastline it was not meant to move?
#
# QA passing is not enough: a rule added for one port can move another
# port's coastline and still pass every gate (the Odaiba step rule moved the
# Futtsu coast up to 47 m off OSM with 0 violations). This compares the solid
# boundary of two refinements of the same recipe, and each one's distance to
# its own filtered OSM coast where they differ.
#
#   python notebooks/435_compare_boundaries.py <before dir> <after dir>
#
# Exit status 0 when the boundaries are identical, 1 when they differ (read
# the report: a difference can be intended), 2 on bad input.
import glob
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import geopandas as gpd  # noqa: E402
import shapely  # noqa: E402

from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import boundary_segments  # noqa: E402


def solid(d: Path):
    f14 = glob.glob(str(d / "*.14"))
    if not f14:
        raise SystemExit(f"no mesh in {d}")
    m = read_fort14(f14[0])
    s, _ = boundary_segments(m.nodes, m.elements, m.open_boundaries)
    return shapely.MultiLineString(list(s)), m


def osm_distance(line, d: Path):
    shp = d / "shoreline_filtered.shp"
    if not shp.exists() or line.is_empty:
        return None
    coast = shapely.union_all(gpd.read_file(shp).geometry.values).boundary
    xy = shapely.get_coordinates(shapely.segmentize(line, 2.0))
    return shapely.distance(shapely.points(xy), coast)


def main(before: Path, after: Path) -> int:
    a, ma = solid(before)
    b, mb = solid(after)
    h = a.hausdorff_distance(b)
    for tag, d, m, line in (("before", before, ma, a), ("after ", after, mb, b)):
        print(f"{tag} {d.name}: NP={m.n_nodes:,} NE={m.n_elements:,}, "
              f"boundary {line.length:.1f} m")
    print(f"solid boundary: Hausdorff {h:.2f} m")
    if h < 1e-6:
        print("identical")
        return 0
    diff = shapely.symmetric_difference(a.buffer(0.01), b.buffer(0.01))
    pieces = sorted((g for g in getattr(diff, "geoms", [diff]) if g.area > 1.0),
                    key=lambda g: -g.area)
    # group pieces that belong to one change (within 50 m of each other)
    zones = shapely.union_all([g.buffer(50.0) for g in pieces])
    for z in getattr(zones, "geoms", [zones]):
        x0, y0, x1, y1 = z.bounds
        box = shapely.box(x0, y0, x1, y1)
        la, lb = shapely.intersection(a, box), shapely.intersection(b, box)
        da, db = osm_distance(la, before), osm_distance(lb, after)
        c = z.centroid
        print(f"  changed near ({c.x:.0f}, {c.y:.0f}), {max(x1 - x0, y1 - y0) - 100:.0f} m across, "
              f"moved up to {la.hausdorff_distance(lb):.1f} m")
        if da is not None and db is not None:
            print(f"    distance to the filtered OSM coast: before median {np.median(da):.1f} / "
                  f"max {da.max():.1f} m, after median {np.median(db):.1f} / max {db.max():.1f} m")
    return 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__ or "usage: 435_compare_boundaries.py <before dir> <after dir>")
        raise SystemExit(2)
    raise SystemExit(main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()))
