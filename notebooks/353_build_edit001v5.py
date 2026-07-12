# edit_001 rev 5 (owner 2026-07-12): SAMPLE-FREE manual edit.
# The general pipeline must not manufacture waterways from a
# reference mesh -- waterways come from the OSM-native detector.
# This edit remains for one reason only: OSM DATA IS WRONG at the
# Haneda D-runway (the pier-supported flow passage is drawn as
# land). Inputs: a hand-drawn guide arc through the passage +
# cross-section measurements against the OSM shoreline. Kept
# waterway rule applies: widen to two standard rows.
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.ops import unary_union

from fvcom_mesh_tools.channel_arcs import snap_arc_to_channel

H0 = 290.0
TWO_ROWS = round(0.875 * 2.0 * 1.2 * H0)      # ~609 m
COSW = float(np.cos(np.deg2rad(35.35)))
SCALE = (111e3 * COSW, 111e3)

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))

# hand-drawn guide through the pier passage (Tama mouth -> bay,
# under the D-runway pier section); NOT derived from any mesh
GUIDE = np.array([
    [139.79566, 35.52345], [139.80061, 35.52521],
    [139.80649, 35.52615], [139.8133, 35.52628]])

snap = snap_arc_to_channel(land, GUIDE, metric_scale=SCALE,
                           step_m=120.0)
w = np.maximum(snap["width_carve_m"], float(TWO_ROWS))
print("stations:", len(snap["arc"]),
      " widths(m):", np.round(w).astype(int).tolist(), flush=True)

from shapely.geometry import LineString
scale = 0.5 * (SCALE[0] + SCALE[1])
crossing = LineString(snap["arc"]).intersection(land)
cross_m = float(crossing.length) * scale
print(f"arc-on-land (pier stretch): {cross_m:.0f} m", flush=True)

ed = {
    "id": "W10-edit001",
    "rev": 5,
    "note": "Haneda D-runway pier passage, rev 5 (SAMPLE-FREE): "
            "OSM draws the pier-supported flow passage as land -- "
            "this manual edit corrects the DATA with a hand-drawn "
            "arc + OSM cross-section widths, widened to two "
            "standard rows like every kept waterway. "
            "arc_on_land_tol_m covers the measured pier stretch "
            "explicitly.",
    "arc": [[round(float(x), 6), round(float(y), 6)]
            for x, y in snap["arc"]],
    "widths_m": [round(float(v), 1) for v in w],
    "min_gap_m": 150.0,
    "arc_on_land_tol_m": round(cross_m + 200.0),
}
out = Path("recipes/edits/sample_repro/"
           "edit_001_haneda_d_runway.json")
out.write_text(json.dumps(ed, indent=1))
print(f"wrote {out} (rev 5, sample-free)", flush=True)
