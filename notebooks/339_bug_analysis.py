# Reproduce the geometry-policy classification and dump the wide
# components + every corridor's contacts: WHY did the Hanami-gawa
# river chain get widened?
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union

poly = np.array(
    [[139.83, 34.973], [140.12, 34.973], [140.12, 35.75],
     [139.60, 35.75], [139.60, 35.20], [139.6642, 35.1546],
     [139.6713, 35.1396], [139.6737, 35.1288], [139.6772, 35.1168],
     [139.6816, 35.1031], [139.6871, 35.0877], [139.6946, 35.0705],
     [139.7000, 35.0576], [139.7069, 35.0445], [139.7134, 35.0327],
     [139.7216, 35.0184], [139.7289, 35.0047], [139.7373, 34.9916],
     [139.7497, 34.9750], [139.83, 34.973]])
H = 290.0 * 1.2
cosw = float(np.cos(np.deg2rad(35.35)))
scale = 0.5 * (111e3 * cosw + 111e3)
r_open = 0.5 * 1.2 * H / scale
a_cell = (np.sqrt(3) / 4) * (H / scale) ** 2

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
dom = Polygon(poly)
water = dom.difference(land)
wide = water.buffer(-r_open).buffer(r_open * 1.02, join_style="mitre",
                                    mitre_limit=1.2)
wide = wide.intersection(water)
parts = [g for g in getattr(wide, "geoms", [wide]) if not g.is_empty]
areas = np.array([g.area for g in parts]) / a_cell
order = np.argsort(-areas)
print("wide components (cells):")
for k in order[:8]:
    c = parts[k].representative_point()
    print(f"  #{k}: {areas[k]:8.1f} cells at ({c.x:.3f},{c.y:.3f})")
obc_pt = shapely.Point(139.7000, 35.0576)
main_i = int(np.argmin([obc_pt.distance(g) for g in parts]))
print(f"main = #{main_i} ({areas[main_i]:.1f} cells)")
# the Hanami-gawa corridor: pick narrow components near 140.10,35.62
narrow = water.difference(wide)
eps = 0.02 * H / scale
for N in getattr(narrow, "geoms", [narrow]):
    c = N.representative_point()
    if 140.05 < c.x < 140.13 and 35.55 < c.y < 35.70:
        nb = [k for k, g in enumerate(parts) if N.distance(g) < eps]
        print(f"river piece at ({c.x:.3f},{c.y:.3f}) "
              f"area={N.area/a_cell:.1f} cells, touches "
              f"{[(k, round(float(areas[k]),1)) for k in nb]}")
