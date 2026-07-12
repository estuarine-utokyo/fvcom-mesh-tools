# edit_001 rev 3 (owner 2026-07-12): THE SAMPLE IS RIGHT. The
# goto2023 sample meshes the whole Haneda D-runway pier region --
# the wedge OSM draws as land AND the sub-cell NW channel -- as a
# broad water band. rev 2 snapped to the wrong (southern) gap and
# left the runway-airport channel 1 cell wide. rev 3 replaces the
# carve with a WATER PATCH: the union of the sample's triangles
# around the site, minus nothing -- land inside that footprint is
# opened exactly as the sample designs it. A medial arc + measured
# widths are stored for the cross-section one-wide checker.
import json
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from pyproj import Transformer
from shapely.geometry import box, mapping
from shapely.ops import unary_union

from fvcom_mesh_tools.channel_arcs import (
    arc_from_points,
    snap_arc_to_channel,
)

COSW = float(np.cos(np.deg2rad(35.35)))
SCALE = (111e3 * COSW, 111e3)
# rev 4 (owner 2026-07-12): the rev-3 window's northern edge
# (35.538) clipped the sample footprint right at cell 3290's
# station, so the NW channel's NE half kept its narrow OSM banks
# and meshed 1-wide. Extend the window over the WHOLE runway and
# the channel's NE exit -- the sample has 2 rows there.
WIN = box(139.788, 35.516, 139.832, 35.552)

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326",
                          always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
P = np.column_stack([lon_s, lat_s])
cent = P[Ts].mean(axis=1)

sel = np.where(
    (cent[:, 0] > 139.788) & (cent[:, 0] < 139.832)
    & (cent[:, 1] > 35.516) & (cent[:, 1] < 35.552))[0]
tris = [shapely.Polygon(P[Ts[i]]) for i in sel]
patch = unary_union(tris).buffer(0).intersection(WIN)
land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
added = patch.intersection(land)
sc = 0.5 * (SCALE[0] + SCALE[1])
print(f"sample triangles in window: {len(sel)}; patch "
      f"{patch.area*sc*sc/1e4:.0f} ha; OSM land opened "
      f"{added.area*sc*sc/1e4:.1f} ha", flush=True)

# medial arc of the band (bank-hugging sample cells only) for the
# cross-section checker: diameter path of centroid kNN graph
d_land = np.array([shapely.Point(c).distance(land) * sc
                   for c in cent[sel]])
band = sel[d_land < 700.0]
arc = arc_from_points(
    np.column_stack([cent[band, 0] * COSW, cent[band, 1]]),
    smooth_passes=2)
arc = np.column_stack([arc[:, 0] / COSW, arc[:, 1]])
land_v3 = land.difference(patch)
snap = snap_arc_to_channel(land_v3, arc, metric_scale=SCALE,
                           step_m=120.0)
print("checker arc stations:", len(snap["arc"]),
      " widths(m):",
      np.round(snap["width_m"]).astype(int).tolist(), flush=True)

# second checker arc: the NW channel between the airport main
# island and the remaining runway body (the owner's 1-wide
# complaint site) -- the medial diameter path does not pass
# through this branch, so it needs its own cross-section arc
guide_nw = np.array([
    [139.7947, 35.5260], [139.7997, 35.5290],
    [139.8055, 35.5318], [139.8135, 35.5365],
    [139.8205, 35.5420]])
snap_nw = snap_arc_to_channel(land_v3, guide_nw,
                              metric_scale=SCALE, step_m=120.0)
print("NW-branch checker stations:", len(snap_nw["arc"]),
      " widths(m):",
      np.round(snap_nw["width_m"]).astype(int).tolist(),
      flush=True)

ed = {
    "id": "W10-edit001",
    "rev": 4,
    "type": "water_patch",
    "note": "Haneda D-runway area, rev 4: adopt the goto2023 "
            "sample's water footprint verbatim (owner: the sample "
            "is right -- the pier region must be meshed). The "
            "patch is the union of sample triangles in the window "
            "(139.788-139.832, 35.516-35.552; rev 4 extends over "
            "the whole runway + NE channel exit -- rev 3 clipped "
            "the footprint at cell 3290 and left it 1-wide); OSM "
            "land inside it "
            "(pier wedge + channel banks) becomes water. The "
            "stored arc/widths drive the one-wide cross-section "
            "checker, not a carve.",
    "geometry": mapping(patch),
    "arc": [[round(float(x), 6), round(float(y), 6)]
            for x, y in snap["arc"]],
    "widths_m": [round(float(v), 1) for v in snap["width_m"]],
    "check_arcs": [
        [[round(float(x), 6), round(float(y), 6)]
         for x, y in snap["arc"]],
        [[round(float(x), 6), round(float(y), 6)]
         for x, y in snap_nw["arc"]],
    ],
    "check_widths_m": [
        [round(float(v), 1) for v in snap["width_m"]],
        [round(float(v), 1) for v in snap_nw["width_m"]],
    ],
}
out = Path("recipes/edits/sample_repro/"
           "edit_001_haneda_d_runway.json")
out.write_text(json.dumps(ed, indent=1))
print(f"wrote {out}", flush=True)
