"""OW05 axis-centering check (owner rule 2026-07-15): corridor-
centre offset from the REAL channel axis at fixed longitudes.
Acceptance metric for replacing edit_004/edit_005 with the
automatic normalization pass -- the manual-edit result measured
-3 m @139.830 and -34 m @139.833 (commit 3dfa621)."""

import json
import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from pyproj import Transformer
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.io import read_fort14

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
LONS = [139.827, 139.830, 139.833]
LAT0, LAT1, LATC = 35.620, 35.645, 35.6320

land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
# SR_AXIS_MESH: explicit override for A/B comparisons (loud)
_mp = os.environ.get("SR_AXIS_MESH")
_mesh_path = Path(_mp) if _mp else OUT / "sample_repro_final.14"
print(f"[372] mesh = {_mesh_path}", flush=True)
mesh = read_fort14(_mesh_path)
tr = Transformer.from_crs(32654, 4326, always_xy=True)
lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
P = np.column_stack([lon, lat])
T = mesh.elements
cent = P[T].mean(axis=1)
win = ((cent[:, 0] > 139.815) & (cent[:, 0] < 139.852)
       & (cent[:, 1] > 35.615) & (cent[:, 1] < 35.650))
mwater = unary_union([shapely.Polygon(P[t]) for t in T[win]])


def _intervals(geom_1d):
    out = []
    for s in getattr(geom_1d, "geoms", [geom_1d]):
        if s.is_empty or s.geom_type != "LineString":
            continue
        ys = [c[1] for c in s.coords]
        out.append((float(min(ys)), float(max(ys))))
    return sorted(out)


def _pick(ivs, lat_ref):
    """Interval containing lat_ref, else the nearest by centre."""
    for iv in ivs:
        if iv[0] <= lat_ref <= iv[1]:
            return iv, "contains"
    if not ivs:
        return None, "none"
    j = int(np.argmin([abs(0.5 * (a + b) - lat_ref)
                       for a, b in ivs]))
    return ivs[j], "nearest"


rows = []
print(f"{'lon':>9} {'axis_lat':>9} {'axis_w_m':>8} "
      f"{'mesh_lat':>9} {'mesh_w_m':>8} {'offset_m':>8}  flags",
      flush=True)
for L in LONS:
    line = shapely.LineString([(L, LAT0), (L, LAT1)])
    ax_iv, ax_how = _pick(_intervals(line.difference(land)), LATC)
    if ax_iv is None:
        print(f"{L:9.3f}  NO ORIGINAL WATER on transect", flush=True)
        continue
    axc = 0.5 * (ax_iv[0] + ax_iv[1])
    m_iv, m_how = _pick(_intervals(line.intersection(mwater)), axc)
    if m_iv is None:
        print(f"{L:9.3f} {axc:9.5f} "
              f"{(ax_iv[1] - ax_iv[0]) * 111e3:8.0f} "
              f"  NO MESH WATER on transect", flush=True)
        rows.append({"lon": L, "axis_lat": round(axc, 5),
                     "offset_m": None})
        continue
    mc = 0.5 * (m_iv[0] + m_iv[1])
    off = (mc - axc) * 111e3
    flags = ("" if ax_how == "contains" else "axis:nearest ") + \
            ("" if m_how == "contains" else "mesh:nearest")
    print(f"{L:9.3f} {axc:9.5f} "
          f"{(ax_iv[1] - ax_iv[0]) * 111e3:8.0f} {mc:9.5f} "
          f"{(m_iv[1] - m_iv[0]) * 111e3:8.0f} {off:8.0f}  "
          f"{flags}", flush=True)
    rows.append({
        "lon": L, "axis_lat": round(axc, 5),
        "axis_w_m": round((ax_iv[1] - ax_iv[0]) * 111e3),
        "mesh_lat": round(mc, 5),
        "mesh_w_m": round((m_iv[1] - m_iv[0]) * 111e3),
        "offset_m": round(off), "flags": flags.strip()})
(OUT / "axis_check.json").write_text(json.dumps(rows, indent=1))
print(f"[372] saved {OUT / 'axis_check.json'}", flush=True)
