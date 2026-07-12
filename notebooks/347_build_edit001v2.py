# Rebuild edit_001 (W10 Haneda D-runway channel) as a FAITHFUL
# arc-band edit (owner 2026-07-12): the v1 fixed 609 m corridor
# eroded the runway because the guide arc (smoothed sample
# centroids) ran south of the real channel. v2:
#   * snap the arc onto the real OSM water gap (cross-sections),
#   * carve with the MEASURED per-station width -> land is only
#     removed at the pier/bridge pinch (data correction), the true
#     banks stay untouched,
#   * constrain BOTH banks with pfix+egfix so the sub-cell-width
#     stretch is meshed as an explicit 1-row band instead of being
#     bridged by DistMesh.
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from fvcom_mesh_tools.channel_arcs import (
    bank_chains,
    snap_arc_to_channel,
)

H0 = 290.0
H_MESH = 1.2 * H0                  # ~350 m rendered edge
SPACING = round(0.9 * H_MESH)      # bank node spacing ~313 m
COSW = float(np.cos(np.deg2rad(35.35)))
SCALE = (111e3 * COSW, 111e3)

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))

# guide arc = v1 arc (rough, from sample centroids)
GUIDE = np.array([
    [139.79566, 35.52345], [139.80061, 35.52521],
    [139.80649, 35.52615], [139.8133, 35.52628]])

snap = snap_arc_to_channel(land, GUIDE, metric_scale=SCALE,
                           step_m=120.0)
w = snap["width_m"]
print("station widths (m):",
      np.round(w).astype(int).tolist(), flush=True)
print("snapped:", snap["snapped"].astype(int).tolist(), flush=True)
pfix, egfix = bank_chains(snap, spacing_m=SPACING,
                          metric_scale=SCALE)
print(f"bank constraints: {len(pfix)} nodes, {len(egfix)} edges",
      flush=True)

ed = {
    "id": "W10-edit001",
    "rev": 2,
    "note": "Haneda D-runway channel, FAITHFUL rebuild: arc "
            "snapped to the real OSM water gap; carve width = "
            "measured per-station width, so land is removed only "
            "at the pier stretch OSM draws as land (~2 % of v1's "
            "erosion); both banks pfix+egfix-constrained so the "
            "sub-cell-width stretch meshes as an explicit 1-row "
            "band. 1-wide cells there are EXPECTED and flagged by "
            "the one-wide reporter for the next manual-edit round.",
    "arc": [[round(float(x), 6), round(float(y), 6)]
            for x, y in snap["arc"]],
    "widths_m": [round(float(v), 1) for v in snap["width_carve_m"]],
    "min_gap_m": 150.0,
    "bank_pfix": [[round(float(x), 6), round(float(y), 6)]
                  for x, y in pfix],
    "bank_egfix": [[int(a), int(b)] for a, b in egfix],
}
out = Path("recipes/edits/sample_repro/"
           "edit_001_haneda_d_runway.json")
out.write_text(json.dumps(ed, indent=1))
print(f"wrote {out}", flush=True)

# faithfulness check: how much ORIGINAL land does v2 remove?
from fvcom_mesh_tools.channel_arcs import carve_channel_corridor
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
_, info = carve_channel_corridor(
    land, np.asarray(ed["arc"], float),
    np.asarray(ed["widths_m"], float),
    min_gap_m=150.0, metric_scale=SCALE, domain_poly=dom)
print(f"v2 carve: len={info['arc_length_m']:.0f} m "
      f"arc_on_land={info['arc_on_land_m']:.0f} m "
      f"land_removed={info['land_removed_m2']/1e4:.2f} ha "
      f"(v1 removed ~{876667 * 0.35 / 1e4:.0f}+ ha)", flush=True)
