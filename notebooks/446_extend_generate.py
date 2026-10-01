# Mesh the sea outside a base mesh, out to a new open boundary.
#
#   python notebooks/446_extend_generate.py RECIPE OUTDIR
#
# Stage 1 of notebooks/445_extend_mesh.py. The base (a finished FVCOM case)
# is treated as land here, so the only sea is the new one; the base's open
# boundary nodes and the new open boundary's nodes are fixed points and edges
# (pfix/egfix), exactly where the inputs put them. Writes OUTDIR/outer_utm.14
# (UTM 54N; open boundary 0 = the new boundary, open boundary 1 = the base
# interface) and OUTDIR/generate.json.
#
# oceanmesh is GPL: it is imported here, in a script run as a subprocess,
# never from the fvcom_mesh_tools package.
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import geopandas as gpd  # noqa: E402
import oceanmesh as om  # noqa: E402
import shapely  # noqa: E402
from oceanmesh import Shoreline  # noqa: E402
from pyproj import Transformer  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from shapely.prepared import prep  # noqa: E402

from fvcom_mesh_tools.algorithms.obc_finish import prune_one_wide_protected  # noqa: E402
from fvcom_mesh_tools.base_recipe import read_open_boundary  # noqa: E402
from fvcom_mesh_tools.dem.sources import sample  # noqa: E402
from fvcom_mesh_tools.extend import band_field, compose_sizing, land_segments  # noqa: E402
from fvcom_mesh_tools.extend_recipe import load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, write_fort14  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402
from fvcom_mesh_tools.obc_band import build_obc_band  # noqa: E402

T0 = time.time()
MESH_EPSG = 32654
DEG = 1.0 / 111e3


def say(msg):
    print(f"[gen] {msg} +{time.time() - T0:.0f}s", flush=True)


recipe = load_extend_recipe(sys.argv[1])
OUT = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
S = recipe["settings"]
to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
to_m = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)

# ------------------------------------------------------------------ base
b = Path(recipe["base"]) / recipe["base_case"]
base = read_fvcom_case(f"{b}_grd.dat", f"{b}_dep.dat", f"{b}_obc.dat")
if len(base.open_boundaries) != 1:
    raise SystemExit(f"the base has {len(base.open_boundaries)} open boundaries; one is needed")
IB = np.asarray(base.open_boundaries[0], np.int64)
iface_ll = np.column_stack(to_ll.transform(base.nodes[IB, 0], base.nodes[IB, 1]))
foot_m = shapely.union_all(shapely.polygons(base.nodes[base.elements][:, :, :2]))
foot_ll = shapely.ops.transform(lambda x, y, z=None: to_ll.transform(x, y), foot_m)
say(f"base {recipe['base_case']}: NP={base.n_nodes:,} NE={base.n_elements:,}, "
    f"interface {len(IB)} nodes")

# ----------------------------------------------------------------- inputs
OBC = np.asarray(read_open_boundary(recipe["open_boundary"]), float)
bb = tuple(recipe["land"]["bbox"])
DATA = Path(__import__("os").environ["DATA_DIR"])
land = gpd.read_file(DATA / "geodata/OSM/land-polygons-split-4326/land_polygons.shp",
                     bbox=bb).clip(bb)
# the base is land for this stage: the only sea to mesh is the new one
land_all = shapely.union_all([*land.geometry.values, foot_ll.buffer(0)])
LAND_SHP = OUT / "land_with_base.shp"
gpd.GeoDataFrame(geometry=list(getattr(land_all, "geoms", [land_all])),
                 crs="EPSG:4326").to_file(LAND_SHP)
# the meshing domain: the new boundary, closed over land to the north
north = bb[3] - 0.05
poly = np.vstack([OBC, [OBC[-1, 0], north], [OBC[0, 0], north]])
if not shapely.Polygon(poly).is_valid:
    raise SystemExit("the domain polygon (boundary closed northward) is not simple")
say(f"inputs: open boundary {len(OBC)} nodes, land {len(land)} polygon(s)")

# ----------------------------------------------------------------- sizing
sh = Shoreline(str(LAND_SHP), poly, S["lattice_m"] * DEG)
sdf = om.signed_distance_function(sh)
fd = om.distance_sizing_function(sh, rate=S["gradation"], max_edge_length=None)
vals = np.ma.filled(np.ma.asarray(fd.values), S["max_edge_m"] * DEG)
d_coast = np.clip((vals / DEG - S["lattice_m"]) / S["gradation"], 0, None)
lon_g, lat_g = fd.create_grid()
lat0, lon0 = float(lat_g.mean()), float(lon_g.mean())
kx = 111e3 * np.cos(np.radians(lat0))
x, y = (lon_g - lon0) * kx, (lat_g - lat0) * 111e3
ambient = np.minimum(S["coast_h_m"] + S["gradation"] * d_coast, S["max_edge_m"])
depth_g, _ = sample(recipe["bathymetry"]["sizing"], lon_g, lat_g)
# sea that no sizing source covers would get depth 0 and no time-step floor;
# refuse it (review F6). Land is depth <= 0 by itself.
wet = sdf.eval(np.column_stack([lon_g.ravel(), lat_g.ravel()])).reshape(lon_g.shape) < 0
if np.isnan(depth_g[wet]).any():
    raise SystemExit(f"{int(np.isnan(depth_g[wet]).sum())} sea lattice point(s) outside every "
                     f"sizing source {recipe['bathymetry']['sizing']}")
depth_g = np.clip(np.nan_to_num(depth_g), 0, None)
floor = S["cfl_dt_s"] * np.sqrt(9.81 * depth_g) / S["cfl_cr"]
say(f"lattice {lon_g.shape}, depth {depth_g.max():.0f} m max, floor up to {floor.max():.0f} m")


def metric(ll):
    return np.column_stack([(ll[:, 0] - lon0) * kx, (ll[:, 1] - lat0) * 111e3])


def spacing(xy):
    e = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    return np.r_[e[0], 0.5 * (e[1:] + e[:-1]), e[-1]]


iface_m = np.column_stack(base.nodes[IB, 0:2].T)          # UTM, the base's own lengths
bands = [band_field(x, y, metric(iface_ll), spacing(iface_m), S["interface_band_m"]),
         band_field(x, y, metric(OBC), spacing(metric(OBC)), S["obc_band_m"])]
for name, line in (("interface", iface_ll), ("open boundary", OBC)):
    # each constrained line gets a ladder, which needs six nodes (review F12)
    if len(line) < 6:
        raise SystemExit(f"the {name} has {len(line)} nodes; its ladder needs 6 or more")
h, srep = compose_sizing(ambient, x, y, grade=S["gradation"], floor=floor, bands=bands)
# the new open boundary was designed at or above the floor (444): a band below
# it there means the design and this sizing disagree (review F4). The base
# interface keeps the base's own spacing, which may be below the floor: the
# base, not the extension, then sets the time step.
if srep.get("band_1_below_floor_cells", 0):
    raise SystemExit(f"{srep['band_1_below_floor_cells']} open-boundary band cell(s) are set "
                     "below the time-step floor; redesign the boundary spacing (444)")
fd.values = h / S["dm_scale"] * DEG
fd.build_interpolant()
say("sizing " + json.dumps(srep))
np.savez_compressed(OUT / "sizing.npz", lon=lon_g, lat=lat_g, h=h, floor=floor, depth=depth_g)

# ------------------------------------------------------------- generation
# Each constrained line gets a LADDER: a second fixed line inside it, one
# local size away (obc_band.build_obc_band, as the base's own open boundary
# has). Without it DistMesh projects free points onto the constrained line,
# and a free node sitting on a fixed edge makes a zero-area sliver there --
# the third build had 13 such on the new boundary's south side.
land_test = prep(land_all)


def ladder(arc_ll, h_m):
    """The inner guide line of ``arc_ll`` (sea on its left), off land."""
    band = build_obc_band(arc_ll, h_m, k_offset=1.25, skip_ends=2, taper="local")
    inner = band["inner_ll"]
    keep = np.array([not land_test.intersects(shapely.Point(q).buffer(0.25 * h / 111e3))
                     for q, h in zip(inner, h_m[2:len(h_m) - 2])])
    return inner, keep


pts = [iface_ll, OBC]
segs = []
off = 0
for arc in pts:
    segs.append(off + np.column_stack([np.arange(len(arc) - 1), np.arange(1, len(arc))]))
    off += len(arc)
ladders = {}
def extension_on_left(arc_ll, h_m):
    """Is the extension (the generation domain) on the left of ``arc_ll``?

    The ladder is built on the left; which way an input line runs is not
    guaranteed (review round 2 F16), so each line is tested: points a
    quarter of the local size to the left of the middle edges must lie in
    the domain (sdf < 0), or to the right if not.
    """
    xy = metric(arc_ll)
    mid = 0.5 * (xy[1:] + xy[:-1])
    t = np.diff(xy, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    left = np.column_stack([-t[:, 1], t[:, 0]])
    d = 0.25 * 0.5 * (h_m[1:] + h_m[:-1])
    probe_m = mid + left * d[:, None]
    probe = np.column_stack([probe_m[:, 0] / kx + lon0, probe_m[:, 1] / 111e3 + lat0])
    inside = sdf.eval(probe) < 0
    if inside.mean() not in (0.0, 1.0) and abs(inside.mean() - 0.5) < 0.25:
        raise SystemExit("cannot tell which side of a constrained line the extension is on")
    return bool(inside.mean() > 0.5)


for name, arc, h_arc in (("interface", iface_ll, spacing(iface_m)),
                         ("open boundary", OBC, spacing(metric(OBC)))):
    flip = not extension_on_left(arc, h_arc)
    a = arc[::-1] if flip else arc
    hh = h_arc[::-1] if flip else h_arc
    inner, keep = ladder(a, hh)
    ladders[name] = {"n_inner": int(len(inner)), "n_kept": int(keep.sum()), "reversed": flip}
    idx = np.flatnonzero(keep)
    pts.append(inner[idx])
    run = [(i, j) for i, j in zip(range(len(idx) - 1), range(1, len(idx)))
           if idx[j] - idx[i] == 1]
    if run:
        segs.append(off + np.array(run))
    off += len(idx)
PFIX = np.vstack(pts)
SEGS = np.vstack(segs)
ni = len(iface_ll)
say("ladders " + json.dumps(ladders))
# cleanup="none": oceanmesh's default clean deletes low-quality boundary
# faces with no protection for fixed points (only faces carrying a fixed edge
# are spared), and the first build lost 14 of the new boundary's nodes that
# way. The finishing chain in 447 repairs quality with the lines held fixed.
p, t = om.generate_mesh(sdf, fd, max_iter=int(S["max_iter"]), seed=int(S["gen_seed"]),
                        pfix=PFIX, egfix=SEGS, cleanup="none")
ne0 = len(t)


def n_lost(p):
    return int((cKDTree(p).query(PFIX)[0] > 1e-8).sum())


def n_degenerate(p, t):
    a = p[t] * [111e3 * np.cos(np.radians(lat0)), 111e3]
    ar = 0.5 * ((a[:, 1, 0] - a[:, 0, 0]) * (a[:, 2, 1] - a[:, 0, 1])
                - (a[:, 2, 0] - a[:, 0, 0]) * (a[:, 1, 1] - a[:, 0, 1]))
    return int((np.abs(ar) < 1e3).sum())


say(f"DistMesh: NP={len(p):,} NE={len(t):,}, fixed nodes lost {n_lost(p)}, "
    f"near-zero-area elements {n_degenerate(p, t)}")

# REPAIR along the constrained lines. DistMesh projects points that step
# outside the domain back onto its boundary -- the open boundary included --
# so a free node can end up on a fixed edge or 0.1 m from a fixed node; each
# makes a zero-area element (the fourth build: 39, at the same places with or
# without the ladders). oceanmesh's default clean removes them but takes fixed
# nodes with them. Here the free nodes within 0.3 of a fixed edge's length of
# that edge are dropped and the rest re-triangulated with the same CGAL
# constrained Delaunay and the same inside test DistMesh uses.
from _constrained_delaunay_class import ConstrainedDelaunayTriangulation as CDT  # noqa: E402
from oceanmesh.fix_mesh import fix_mesh  # noqa: E402

kxy = np.array([111e3 * np.cos(np.radians(lat0)), 111e3])
d_fix, i_fix = cKDTree(p).query(PFIX)
free = np.ones(len(p), bool)
free[i_fix[d_fix < 1e-8]] = False
q = p[free] * kxy
too_close = np.zeros(len(q), bool)
for a_, b_ in SEGS:
    A, B = PFIX[a_] * kxy, PFIX[b_] * kxy
    ab = B - A
    L2 = float(ab @ ab)
    s = np.clip(((q - A) @ ab) / L2, 0, 1)
    dist = np.linalg.norm(q - (A + s[:, None] * ab), axis=1)
    too_close |= dist < 0.3 * np.sqrt(L2)
keep_free = p[free][~too_close]
pts_all = np.vstack([PFIX, keep_free])
dt_ = CDT()
dt_.insert(pts_all.ravel().tolist())
dt_.insert_constraints(np.hstack([PFIX[SEGS[:, 0]], PFIX[SEGS[:, 1]]]).ravel().tolist())
p, t = dt_.get_finite_vertices(), dt_.get_finite_cells()
geps = 1e-12 * float(np.amin(S["lattice_m"] * DEG))
t = t[sdf.eval(p[t].sum(1) / 3) < -geps]
# three consecutive nodes of a straight fixed line form a flat triangle whose
# centroid sits a hair inside the domain; kept, it buries the middle node
# (the fifth build: 33 m2 on 3 km edges, open-boundary node 10 off the loop)
_a = p[t] * kxy
_area = 0.5 * np.abs((_a[:, 1, 0] - _a[:, 0, 0]) * (_a[:, 2, 1] - _a[:, 0, 1])
                     - (_a[:, 2, 0] - _a[:, 0, 0]) * (_a[:, 1, 1] - _a[:, 0, 1]))
_l2 = sum(((_a[:, i] - _a[:, (i + 1) % 3]) ** 2).sum(1) for i in range(3))
flat = 4 * np.sqrt(3) * _area / _l2 < 0.01
t = t[~flat]
p, t, _ = fix_mesh(p, t, dim=2, delete_unused=True)
say(f"repair: dropped {int(too_close.sum())} free node(s) at the fixed lines and "
    f"{int(flat.sum())} flat element(s); "
    f"NP={len(p):,} NE={len(t):,}, fixed nodes lost {n_lost(p)}, "
    f"near-zero-area elements {n_degenerate(p, t)}")
p, t = prune_one_wide_protected(p, t, PFIX)
say(f"after one-wide pruning: NE={len(t):,}, fixed nodes lost {n_lost(p)}")
p, t = om.make_mesh_boundaries_traversable(p, t)
say(f"after boundary cleanup: NE={len(t):,}, fixed nodes lost {n_lost(p)}")
say(f"generated NP={len(p):,} NE={len(t):,} (pruned {ne0 - len(t)})")
d, idx = cKDTree(p).query(PFIX)
if (d > 1e-8).any():
    np.savez(OUT / "generate_failed.npz", p=p, t=t, pfix=PFIX, segs=SEGS)
    lost = np.flatnonzero(d > 1e-8)
    for k in lost:
        which = ("interface" if k < ni else "open boundary" if k < ni + len(OBC)
                 else "ladder")
        print(f"[gen]   lost pfix {k} ({which}) at {PFIX[k].round(5).tolist()}, "
              f"nearest node {d[k] / DEG:.0f} m away", flush=True)
    raise SystemExit(f"{len(lost)} fixed boundary node(s) lost in generation "
                     f"(mesh kept in {OUT / 'generate_failed.npz'})")
chain_i, chain_o = idx[:ni], idx[ni:ni + len(OBC)]
edges = {frozenset(e) for e in np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]).tolist()}
for name, c in (("interface", chain_i), ("open boundary", chain_o)):
    miss = [(int(a), int(b_)) for a, b_ in zip(c[:-1], c[1:]) if frozenset((a, b_)) not in edges]
    if miss:
        raise SystemExit(f"{name}: {len(miss)} fixed edge(s) are not mesh edges")
_e = np.sort(np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1)
_u, _c = np.unique(_e, axis=0, return_counts=True)
on_boundary = set(_u[_c == 1].ravel().tolist())
for name, c in (("interface", chain_i), ("open boundary", chain_o)):
    buried = [k for k, v in enumerate(c) if int(v) not in on_boundary]
    if buried:
        raise SystemExit(f"{name}: node(s) {buried} are not on the mesh boundary")

# ------------------------------------------------------------------ write
xu, yu = to_m.transform(p[:, 0], p[:, 1])
nodes = np.column_stack([xu, yu])
nodes[chain_i] = base.nodes[IB, :2]                         # the base's exact coordinates
dn, _ = sample(recipe["bathymetry"]["sizing"], p[:, 0], p[:, 1])
# every generated node is sea: one no source covers must not get an invented
# depth (review round 2 F19)
if np.isnan(dn).any():
    raise SystemExit(f"{int(np.isnan(dn).sum())} generated node(s) outside every sizing "
                     f"source {recipe['bathymetry']['sizing']}")
dn = np.clip(dn, 2.0, None)
mesh = Fort14Mesh("outer", nodes, dn, t.astype(np.int64),
                  [chain_o.astype(np.int64), chain_i.astype(np.int64)],
                  land_segments(t, [chain_o, chain_i]))
write_fort14(mesh, OUT / "outer_utm.14")
(OUT / "generate.json").write_text(json.dumps({
    "n_nodes": int(len(p)), "n_elements": int(len(t)), "pruned": int(ne0 - len(t)),
    "ladders": ladders,
    "interface_base_nodes": IB.tolist(), "interface_outer_nodes": chain_i.tolist(),
    "open_boundary_outer_nodes": chain_o.tolist(), "sizing": srep,
    "lattice_shape": list(lon_g.shape), "settings": S,
}, indent=1, default=float))
say(f"wrote {OUT / 'outer_utm.14'}")
