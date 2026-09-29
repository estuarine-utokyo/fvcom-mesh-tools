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

from fvcom_mesh_tools.algorithms.obc_finish import prune_one_wide_protected  # noqa: E402
from fvcom_mesh_tools.base_recipe import read_open_boundary  # noqa: E402
from fvcom_mesh_tools.dem.sources import sample  # noqa: E402
from fvcom_mesh_tools.extend import band_field, compose_sizing, land_segments  # noqa: E402
from fvcom_mesh_tools.extend_recipe import load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, write_fort14  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402

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
h, srep = compose_sizing(ambient, x, y, grade=S["gradation"], floor=floor, bands=bands)
fd.values = h / S["dm_scale"] * DEG
fd.build_interpolant()
say("sizing " + json.dumps(srep))

# ------------------------------------------------------------- generation
PFIX = np.vstack([iface_ll, OBC])
ni = len(iface_ll)
SEGS = np.vstack([np.column_stack([np.arange(ni - 1), np.arange(1, ni)]),
                  ni + np.column_stack([np.arange(len(OBC) - 1), np.arange(1, len(OBC))])])
p, t = om.generate_mesh(sdf, fd, max_iter=int(S["max_iter"]), seed=int(S["gen_seed"]),
                        pfix=PFIX, egfix=SEGS)
ne0 = len(t)


def n_lost(p):
    return int((cKDTree(p).query(PFIX)[0] > 1e-8).sum())


say(f"DistMesh: NP={len(p):,} NE={len(t):,}, fixed nodes lost {n_lost(p)}")
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
        which = "interface" if k < ni else "open boundary"
        print(f"[gen]   lost pfix {k} ({which}) at {PFIX[k].round(5).tolist()}, "
              f"nearest node {d[k] / DEG:.0f} m away", flush=True)
    raise SystemExit(f"{len(lost)} fixed boundary node(s) lost in generation "
                     f"(mesh kept in {OUT / 'generate_failed.npz'})")
chain_i, chain_o = idx[:ni], idx[ni:]
edges = {frozenset(e) for e in np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]).tolist()}
for name, c in (("interface", chain_i), ("open boundary", chain_o)):
    miss = [(int(a), int(b_)) for a, b_ in zip(c[:-1], c[1:]) if frozenset((a, b_)) not in edges]
    if miss:
        raise SystemExit(f"{name}: {len(miss)} fixed edge(s) are not mesh edges")

# ------------------------------------------------------------------ write
xu, yu = to_m.transform(p[:, 0], p[:, 1])
nodes = np.column_stack([xu, yu])
nodes[chain_i] = base.nodes[IB, :2]                         # the base's exact coordinates
dn, _ = sample(recipe["bathymetry"]["sizing"], p[:, 0], p[:, 1])
dn = np.clip(np.nan_to_num(dn, nan=2.0), 2.0, None)
mesh = Fort14Mesh("outer", nodes, dn, t.astype(np.int64),
                  [chain_o.astype(np.int64), chain_i.astype(np.int64)],
                  land_segments(t, [chain_o, chain_i]))
write_fort14(mesh, OUT / "outer_utm.14")
(OUT / "generate.json").write_text(json.dumps({
    "n_nodes": int(len(p)), "n_elements": int(len(t)), "pruned": int(ne0 - len(t)),
    "interface_base_nodes": IB.tolist(), "interface_outer_nodes": chain_i.tolist(),
    "open_boundary_outer_nodes": chain_o.tolist(), "sizing": srep,
    "lattice_shape": list(lon_g.shape), "settings": S,
}, indent=1, default=float))
say(f"wrote {OUT / 'outer_utm.14'}")
