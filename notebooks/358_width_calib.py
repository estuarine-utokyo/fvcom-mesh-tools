# One-shot calibration: mean corridor width vs sample coverage for
# every keep-eligible network (floor disabled), to choose the
# resolve-width floor from DATA (measure-before-recommend).
import os
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from pyproj import Transformer
from fvcom_mesh_tools.waterways import detect_waterways

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
OBC = [[139.6713, 35.1396], [139.6737, 35.1288],
       [139.6772, 35.1168], [139.6816, 35.1031],
       [139.6871, 35.0877], [139.6946, 35.0705],
       [139.7000, 35.0576], [139.7069, 35.0445],
       [139.7134, 35.0327], [139.7216, 35.0184],
       [139.7289, 35.0047], [139.7373, 34.9916],
       [139.7497, 34.9750]]
dom = Polygon([[139.83, 34.973], [140.12, 34.973],
               [140.12, 35.75], [139.60, 35.75], [139.60, 35.20],
               [139.6642, 35.1546]] + OBC + [[139.83, 34.973]])
cosw = float(np.cos(np.deg2rad(35.35)))
recs = detect_waterways(land, dom, h_mesh_m=1.2 * 290,
                        obc_point=(139.7, 35.05),
                        metric_scale=(111e3 * cosw, 111e3),
                        min_resolve_width_frac=0.0)
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326",
                          always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon, lat = tr.transform(Ps[:, 0], Ps[:, 1])
cent = np.column_stack([lon, lat])[Ts].mean(axis=1)
pts = shapely.points(cent[:, 0], cent[:, 1])
rows = []
for r in recs:
    if r["action"] != "keep" or r["extent_cells"] < 2:
        continue
    gbuf = r["geometry"].buffer(100 / 111e3)
    nin = int(np.sum(shapely.covers(gbuf, pts)))
    c = r["geometry"].representative_point()
    rows.append((r["mean_width_cells"], nin, r["kind"],
                 r["extent_cells"], float(c.x), float(c.y)))
print("w_cells  sample_cells  kind  ext  center")
for w, nin, kind, ext, x, y in sorted(rows):
    print(f"w={w:4.2f}  sample={nin:3d}  {kind:8s} ext={ext:5.1f}"
          f"  ({x:.4f},{y:.4f})", flush=True)
