# Diagnose the Keihin network's skeleton branches: per-branch
# natural width + location, to see (a) which thin branches were
# wrongly carved (G8-d4/e4 fabrication), (b) why the wide
# H8-b5 -> H8-c5 reach ended up closed.
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from fvcom_mesh_tools.channel_arcs import (
    skeleton_branches, snap_arc_to_channel)
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
MS = (111e3 * cosw, 111e3)
h = 1.2 * 290
recs = detect_waterways(land, dom, h_mesh_m=h,
                        obc_point=(139.7, 35.05),
                        metric_scale=MS)
# the big Keihin network = keep record whose geometry covers
# (139.7167, 35.4900)
import shapely as sh
p0 = sh.Point(139.7167, 35.4900)
rec = min((r for r in recs if r["action"] == "keep"),
          key=lambda r: r["geometry"].distance(p0))
print("network kind/ext/w:", rec["kind"], rec["extent_cells"],
      rec["mean_width_cells"], flush=True)
brs = skeleton_branches(rec["geometry"].buffer(0),
                        metric_scale=MS, density_m=0.15 * h,
                        prune_m=0.7 * h)
print(f"branches: {len(brs)}", flush=True)
scale = 0.5 * sum(MS)
for bi, br in enumerate(brs):
    ls = sh.LineString(br)
    L = ls.length * scale
    try:
        snap = snap_arc_to_channel(land, br, metric_scale=MS,
                                   step_m=0.35 * h,
                                   max_halfwidth_m=2.5 * h)
        wm = snap["width_m"]
        sat = wm >= 0.95 * 2 * 2.5 * h
        wn = wm[~sat] if (~sat).any() else wm
        w_med = float(np.median(wn))
        note = ""
    except Exception as e:
        w_med = -1
        note = f"SNAP FAIL: {str(e)[:60]}"
    c = ls.interpolate(0.5, normalized=True)
    print(f"  br{bi:02d} len={L:6.0f} m w_med={w_med:6.0f} m "
          f"({w_med/h:5.2f} h) at ({c.x:.4f},{c.y:.4f}) {note}",
          flush=True)
