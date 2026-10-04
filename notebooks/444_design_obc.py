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
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
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
# The products -- the CSV, its report, its figure and a rejection report --
# are four distinct files, none of them the design by name, link or file
# identity (review rounds 17 F3, 18 F2, F3). A .csv suffix keeps the CSV
# from being its own report or figure.
if out_csv.suffix != ".csv":
    raise SystemExit(f"the boundary output must end in .csv, not {out_csv.name}")
PRODUCTS = (out_csv, out_csv.with_suffix(".json"), out_csv.with_suffix(".png"),
            out_csv.with_suffix(".rejected.json"))
for prod in PRODUCTS:
    if (prod == design_path or prod.resolve() == design_path
            or (prod.exists() and prod.samefile(design_path))):
        raise SystemExit(f"{prod} would overwrite the design {design_path}")
# nor may two products be one file (round 19 F5)
_ids = [p.resolve() for p in PRODUCTS]
if len(set(_ids)) != len(_ids) or any(
        a.exists() and b.exists() and a.samefile(b)
        for i, a in enumerate(PRODUCTS) for b in PRODUCTS[i + 1:]):
    raise SystemExit(f"two of {[p.name for p in PRODUCTS]} are the same file")


def _umask() -> int:
    m = os.umask(0)
    os.umask(m)
    return m
from fvcom_mesh_tools.yaml_strict import load_unique  # noqa: E402

cfg = load_unique(design_path.read_text())      # no repeated keys (round 21 F7)
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
# the whole leg list first: `until: end_normal` ends the walk, so it may only be
# the final leg, and nothing after it may be silently dropped (review round 50 F2)
_legs = cfg["legs"]
if not isinstance(_legs, list) or not _legs:
    raise SystemExit("legs must be a non-empty list")
_ends = [i for i, g in enumerate(_legs) if g.get("until") == "end_normal"]
if _ends and _ends != [len(_legs) - 1]:
    raise SystemExit(f"'until: end_normal' must be on the final leg only; found on leg(s) "
                     f"{_ends} of {len(_legs)}")
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
# finite positive controls, or the floor silently vanishes (review round 3 F11)
for key in ("min_m", "cfl_dt_s", "cfl_cr"):
    v = sp.get(key)
    if not (isinstance(v, (int, float)) and not isinstance(v, bool) and np.isfinite(v) and v > 0):
        raise SystemExit(f"spacing.{key} must be a finite positive number, not {v!r}")


def spacing(p):
    lo, la = to_ll.transform(p[:, 0], p[:, 1])
    d, _ = sample(names, np.asarray(lo), np.asarray(la))
    if np.isnan(d).any():
        # uncovered water would get no time-step floor (review F6)
        raise SystemExit(f"{int(np.isnan(d).sum())} boundary point(s) outside every source {names}")
    floor = sp["cfl_dt_s"] * np.sqrt(9.81 * np.clip(d, 0, None)) / sp["cfl_cr"]
    return np.maximum(sp["min_m"], floor)


nodes = resample(line, spacing)
if len(nodes) < 6:
    # the generation puts a ladder inside the boundary, which needs six nodes,
    # and the report's interior statistics need some (review round 3 F12)
    raise SystemExit(f"the boundary has {len(nodes)} nodes; 6 or more are needed")
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
# on the published geometry: resample keeps chords at the floor, so only
# round-off is tolerated (review round 3 F2)
if report["min_edge_over_floor"] < 1 - 1e-6:
    bad.append(f"an edge is below the spacing floor: {report['min_edge_over_floor']:.4f}")
# a line that crosses or touches itself, or repeats a node, cannot bound the
# domain; generation would refuse it only later (review round 4 F10)
if not full.is_simple:
    bad.append("the boundary crosses or touches itself")
if (steps <= 0).any():
    bad.append(f"{int((steps <= 0).sum())} repeated node(s)")
# The land is read inside land_bbox and clipped to it, so the clip edges
# look like coast. The whole boundary, with its end chords and a margin,
# must lie inside the window, where every coast is real (review round 6 F2).
window = shapely.ops.transform(lambda x, y, z=None: to_m.transform(x, y),
                               shapely.box(*bb).segmentize(0.01))
reach = chord + 5000.0
if not full.buffer(reach).within(window):
    bad.append(f"the boundary comes within {reach / 1e3:.0f} km of the land window's edge "
               f"{list(bb)}; widen land_bbox")
report["problems"] = bad
print("[obc] " + json.dumps({k: report[k] for k in (
    "n_nodes", "length_km", "spacing_m", "end_angle_to_coast_deg", "crosses_land_m",
    "min_distance_to_land_km", "depth_m", "bbox_lonlat")}), flush=True)

out_csv.parent.mkdir(parents=True, exist_ok=True)
if bad:
    # a rejected design must not replace a usable boundary (review F8): its
    # report goes beside it under another name, the CSV is left alone
    # written to a temporary file and moved into place: a rejection name
    # that is a link to the published CSV must not write through it
    # (review round 19 F5)
    rej = out_csv.with_suffix(".rejected.json")
    fd, tmp = tempfile.mkstemp(dir=rej.parent, prefix=rej.name + ".", suffix=".tmp")
    with os.fdopen(fd, "w") as fh:
        fh.write(json.dumps(report, indent=1, default=float))
    os.chmod(tmp, 0o666 & ~_umask())
    os.replace(tmp, rej)
    raise SystemExit(f"the design is rejected ({'; '.join(bad)}); report in {rej}, "
                     f"{out_csv.name} left as it was")
header = (f"# Open boundary nodes (lon,lat, EPSG:4326) designed by notebooks/444_design_obc.py\n"
          f"# from {design_path.name}; {len(nodes)} nodes, {report['length_km']:.1f} km.\n"
          "lon,lat\n")
csv_text = header + "".join(f"{a:.9f},{b:.9f}\n" for a, b in zip(lon, lat))
# the report names the CSV it belongs to (review round 2 F2)
report["csv_sha256"] = hashlib.sha256(csv_text.encode()).hexdigest()
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

# One publisher at a time, each with its own temporary files: a shared
# "<name>.tmp" let one writer rename another's payload (review round 2 F2).
# All three products are written to temporaries first; only then is each
# published file kept as "<name>.prev" and replaced. A failure while
# replacing puts the previous set back, so a failed redesign leaves the
# previous boundary usable (round 4 F6). If the process dies half-way, the
# report's hash no longer matches the CSV and the recipe loader refuses it
# (round 3 F5); the ".prev" files hold the previous set.
lock = out_csv.with_name(out_csv.name + ".lock")
try:
    os.close(os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
except FileExistsError:
    raise SystemExit(f"{lock} exists: another design is being published") from None
staged, kept, done = [], {}, []
# Cleanup (backups, lock) happens only after a known outcome: nothing
# replaced yet, all published, or all put back. Any other exit -- a failed
# restore, an interrupt during one -- keeps the .prev backups and the lock
# (review rounds 5-7, 10 F5).
state = "staging"
UMASK = os.umask(0)
os.umask(UMASK)
try:
    if out_csv.with_name(out_csv.name + ".RECOVER").exists():
        raise SystemExit(f"{out_csv.name}.RECOVER: an earlier publication failed half-way; "
                         "restore from its .prev files and remove the marker first")
    for path, payload in ((out_csv.with_suffix(".json"),
                           json.dumps(report, indent=1, default=float)),
                          (out_csv.with_suffix(".png"), fig),
                          (out_csv, csv_text)):
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
        staged.append((Path(tmp), path))
        if isinstance(payload, str):
            with os.fdopen(fd, "w") as fh:
                fh.write(payload)
        else:
            os.close(fd)
            payload.savefig(tmp, dpi=110, bbox_inches="tight", format="png")
        os.chmod(tmp, 0o666 & ~UMASK)      # mkstemp makes 0600
    # an existing backup is an unresolved earlier failure: never overwrite
    # it (review round 7 F11)
    old_prev = [p.with_name(p.name + ".prev") for _, p in staged
                if os.path.lexists(p.with_name(p.name + ".prev"))]   # a dangling link too (r35 F2)
    if old_prev:
        raise SystemExit(f"backups of an earlier failed publication remain: "
                         f"{[q.name for q in old_prev]}; resolve them first")
    state = "publishing"
    for _, path in staged:
        if path.exists():
            prev = path.with_name(path.name + ".prev")
            shutil.copy2(path, prev)
            kept[path] = prev
    try:
        for tmp, path in staged:
            # recorded before the move: an interrupt right after it must
            # still see it (round 10 F6); putting back a file that was not
            # moved restores an identical copy
            done.append(path)
            os.replace(tmp, path)          # atomic on one file system
    except BaseException:
        state = "restoring"
        stuck = {}
        for path in done:
            try:
                if path in kept:
                    os.replace(kept[path], path)
                    del kept[path]
                else:
                    path.unlink(missing_ok=True)
            except OSError as exc:
                stuck[path] = exc
        if stuck:
            out_csv.with_name(out_csv.name + ".RECOVER").write_text(json.dumps(
                {"restore_from_prev": {str(p): str(q) for p, q in kept.items() if p in stuck},
                 "errors": [f"{p}: {e}" for p, e in stuck.items()]}, indent=1))
            raise
        state = "restored"
        raise
    state = "done"
finally:
    for tmp, _ in staged:                  # never-published temporaries
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass
    if state in ("staging", "done", "restored"):
        # copies of files now in place, or never replaced
        for prev in kept.values():
            try:
                prev.unlink(missing_ok=True)
            except OSError:
                pass
        lock.unlink()
print(f"[obc] wrote {out_csv}", flush=True)
