# Design an open boundary from a design recipe and write its nodes.
#
#   python notebooks/444_design_obc.py recipes/extend/tokyo_bay_enshu_obc_design.yaml \
#       recipes/extend/tokyo_bay_enshu_obc.csv
#
# The nodes are the INPUT of mesh generation (owner, 2026-09-28: the open
# boundary is given, not generated). This script is how that input is made
# and checked: orthogonal to the coast at both ends, straight sides, rounded
# corners, and a node spacing no finer than the time-step floor. A report
# (<csv>.json) and a figure (<csv>.png) are written beside the CSV.
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml
from pyproj import Transformer
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from fvcom_mesh_tools.dem.sources import sample, source_files  # noqa: E402
from fvcom_mesh_tools.obc_design import (  # noqa: E402
    bearing_vector,
    coast_normal,
    fillet,
    ray_intersection,
    resample,
)

MESH_EPSG = 32654          # the base mesh's plane (UTM 54N)
design_path = Path(sys.argv[1]).resolve()
out_csv = Path(sys.argv[2]).resolve()
cfg = yaml.safe_load(design_path.read_text())
DATA = Path(os.environ["DATA_DIR"])

to_m = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)
to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)


def ll(p):
    return np.array(to_ll.transform(*np.asarray(p, float).T)).T


# ------------------------------------------------------------------ land
import geopandas as gpd  # noqa: E402
import shapely  # noqa: E402

bb = tuple(cfg["land_bbox"])
land_ll = gpd.read_file(DATA / "geodata/OSM/land-polygons-split-4326/land_polygons.shp",
                        bbox=bb).clip(bb)
land = shapely.ops.transform(lambda x, y, z=None: to_m.transform(x, y),
                             shapely.unary_union(land_ll.geometry.values))
chord = float(cfg.get("chord_m", 3000))
bs, S0 = coast_normal(land, *to_m.transform(*cfg["start"]), chord_m=chord)
be, E0 = coast_normal(land, *to_m.transform(*cfg["end"]), chord_m=chord)
print(f"[obc] start {ll(S0).round(5).tolist()} normal {bs:.1f} deg; "
      f"end {ll(E0).round(5).tolist()} normal {be:.1f} deg", flush=True)

# ------------------------------------------------------------------ sides
verts = [S0]
cur = S0
for k, leg in enumerate(cfg["legs"]):
    b = {"start_normal": bs}.get(leg["bearing"], leg["bearing"])
    u = bearing_vector(float(b))
    if leg.get("until") == "end_normal":
        nxt = ray_intersection(cur, u, E0, bearing_vector(be))
    else:
        (key, target), = [(kk, v) for kk, v in leg.items() if kk.startswith("until_")]
        axis = {"until_lon": 0, "until_lat": 1}[key]

        def miss(s, cur=cur, u=u, axis=axis, target=target):
            return ll(cur + s * u)[axis] - target

        s_hi = 50_000.0
        while miss(s_hi) * miss(0.0) > 0:
            s_hi *= 2
            if s_hi > 5e6:
                raise SystemExit(f"leg {k} never reaches {key} = {target}")
        nxt = cur + brentq(miss, 0.0, s_hi, xtol=0.01) * u
    verts.append(nxt)
    cur = nxt
    if leg.get("until") == "end_normal":
        break
verts.append(E0)
verts = np.array(verts)
line = fillet(verts, cfg["radii_m"])
print(f"[obc] corners {ll(verts[1:-1]).round(4).tolist()}", flush=True)

# ---------------------------------------------------------------- spacing
sp = cfg["spacing"]
names = sp["bathymetry"]


def spacing(p):
    lo, la = to_ll.transform(p[:, 0], p[:, 1])
    d, _ = sample(names, np.asarray(lo), np.asarray(la))
    if np.isnan(d).any():
        # uncovered water would get no time-step floor (review F6)
        raise SystemExit(f"{int(np.isnan(d).sum())} boundary point(s) outside every source {names}")
    floor = sp["cfl_dt_s"] * np.sqrt(9.81 * np.clip(d, 0, None)) / sp["cfl_cr"]
    return np.maximum(sp["min_m"], floor)


nodes = resample(line, spacing)
# what is checked is what is published: the CSV holds lon/lat to 1e-9 deg
# (about 0.1 mm), and the nodes are taken back from those rounded values
# (review round 2 F7)
lon, lat = (np.round(np.asarray(v), 9) for v in to_ll.transform(nodes[:, 0], nodes[:, 1]))
nodes = np.column_stack(to_m.transform(lon, lat))
depth, which = sample(names, np.asarray(lon), np.asarray(lat))

# ----------------------------------------------------------------- checks


def end_angle(first, second):
    p = shapely.Point(first)
    ring = min((g.exterior for g in getattr(land, "geoms", [land])), key=lambda r: r.distance(p))
    s = ring.project(p)
    c = (np.array(ring.interpolate((s + chord) % ring.length).coords[0])
         - np.array(ring.interpolate((s - chord) % ring.length).coords[0]))
    e = np.asarray(second) - np.asarray(first)
    return float(np.degrees(np.arccos(abs(c @ e) / np.linalg.norm(c) / np.linalg.norm(e))))


# the whole line, end edges included: the ends only touch the coast (review F7)
full = shapely.LineString(nodes)
steps = np.linalg.norm(np.diff(nodes, axis=0), axis=1)
# each edge (a chord) against the floor at both of its ends (review F5)
floor_n = spacing(nodes)
edge_floor = np.maximum(floor_n[:-1], floor_n[1:])
report = {
    "design": str(design_path), "csv": str(out_csv),
    "start_normal_deg": bs, "end_normal_deg": be,
    "corners_lonlat": ll(verts[1:-1]).tolist(),
    "n_nodes": int(len(nodes)), "length_km": float(steps.sum() / 1e3),
    "spacing_m": {"min": float(steps.min()), "median": float(np.median(steps)),
                  "max": float(steps.max())},
    "end_angle_to_coast_deg": [end_angle(nodes[0], nodes[1]), end_angle(nodes[-1], nodes[-2])],
    "crosses_land_m": float(land.intersection(full).length),
    "min_edge_over_floor": float((steps / edge_floor).min()),
    "min_distance_to_land_km": float(min(land.distance(shapely.Point(p))
                                         for p in nodes[2:-2]) / 1e3),
    "depth_m": {"min": float(np.nanmin(depth)), "max": float(np.nanmax(depth))},
    "bbox_lonlat": [float(min(lon)), float(min(lat)), float(max(lon)), float(max(lat))],
    "csv_decimals": 9,
    "bathymetry": names,
    "bathymetry_files": {k: [str(p) for p in v] for k, v in source_files(names).items()},
}
bad = []
if any(abs(a - 90) > 1.0 for a in report["end_angle_to_coast_deg"]):
    bad.append(f"not orthogonal to the coast: {report['end_angle_to_coast_deg']}")
# a crossing, not the round-off of an end point lying on the coast
# (review round 2 F6): 1 mm
if report["crosses_land_m"] > 1e-3:
    bad.append(f"crosses land over {report['crosses_land_m']:.0f} m")
if report["min_edge_over_floor"] < 0.995:      # a chord is a little shorter than its arc
    bad.append(f"an edge is below the spacing floor: {report['min_edge_over_floor']:.4f}")
report["problems"] = bad
print("[obc] " + json.dumps({k: report[k] for k in (
    "n_nodes", "length_km", "spacing_m", "end_angle_to_coast_deg", "crosses_land_m",
    "min_distance_to_land_km", "depth_m", "bbox_lonlat")}), flush=True)

out_csv.parent.mkdir(parents=True, exist_ok=True)
if bad:
    # a rejected design must not replace a usable boundary (review F8): its
    # report goes beside it under another name, the CSV is left alone
    rej = out_csv.with_suffix(".rejected.json")
    rej.write_text(json.dumps(report, indent=1, default=float))
    raise SystemExit(f"the design is rejected ({'; '.join(bad)}); report in {rej}, "
                     f"{out_csv.name} left as it was")
header = (f"# Open boundary nodes (lon,lat, EPSG:4326) designed by notebooks/444_design_obc.py\n"
          f"# from {design_path.name}; {len(nodes)} nodes, {report['length_km']:.1f} km.\n"
          "lon,lat\n")
csv_text = header + "".join(f"{a:.9f},{b:.9f}\n" for a, b in zip(lon, lat))
# the report names the CSV it belongs to (review round 2 F2)
report["csv_sha256"] = hashlib.sha256(csv_text.encode()).hexdigest()
# one publisher at a time, each with its own temporary files: a shared
# "<name>.tmp" let one writer rename another's payload (review round 2 F2)
lock = out_csv.with_name(out_csv.name + ".lock")
try:
    os.close(os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
except FileExistsError:
    raise SystemExit(f"{lock} exists: another design is being published") from None
try:
    for path, text in ((out_csv, csv_text),
                       (out_csv.with_suffix(".json"), json.dumps(report, indent=1, default=float))):
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)              # atomic on one file system
finally:
    lock.unlink()

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(figsize=(10, 9))
land_ll.plot(ax=ax, color="0.85", edgecolor="k", linewidth=0.4)
ax.plot(lon, lat, "r-", lw=2, label=f"open boundary ({len(nodes)} nodes)")
ax.plot(lon, lat, "r.", ms=3)
ax.set_xlim(bb[0], bb[2])
ax.set_ylim(bb[1], bb[3])
ax.set_aspect(1 / np.cos(np.radians(np.mean(lat))))
ax.legend(loc="lower right")
ax.set_title(f"{out_csv.name} (black = land, red = open boundary)")
fig.savefig(out_csv.with_suffix(".png"), dpi=110, bbox_inches="tight")
print(f"[obc] wrote {out_csv}", flush=True)
