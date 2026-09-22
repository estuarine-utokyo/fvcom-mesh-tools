# Local refinement of an existing mesh: cut a hole around a declared region
# and fill it at the target size, leaving everything outside bit-for-bit
# unchanged.  The design and its reasoning are in docs/local_refine.md.
#
# The split of work follows the licence policy (CLAUDE.md): selection, the rim,
# the coastline modes, stitching and verification live in the Apache-2.0
# package (fvcom_mesh_tools.patch); the fill is DistMesh from the GPL
# oceanmesh fork and therefore lives here, in a notebook, exactly as
# notebooks/325_sample_repro.py does.
#
#   python notebooks/420_local_refine.py recipes/refine/futtsu_nori.yaml
#
# Environment:
#   FMESH_LAND   land polygons the base mesh was fitted to (any CRS)
#   LR_OUT       output directory
#   LR_MAX_ITER  DistMesh iterations (default 100)
#   LR_SEED      DistMesh seed (default 0)
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import geopandas as gpd
import oceanmesh as om
import shapely
from shapely.ops import unary_union

from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.patch import (
    ambient_size_field,
    effective_gradation,
    hole_polygon,
    improve_patch,
    patch_sizing,
    rim_constraints,
    select_patch,
    stitch_patch,
    verify_patch,
)
from fvcom_mesh_tools.refine import load_refine, preflight, transition_width_m

MESH_EPSG = 32654
# DistMesh returns bars about this much larger than the sizing field
# (mesh_generator.py L0mult; measured 1.19-1.25 by notebook 394).  target_h_m
# is an ACHIEVED edge length, so the field carries target / this -- the same
# SR_H_TARGET=achieved convention as notebook 325.
DISTMESH_SCALE = 1.2
t0 = time.time()
recipe = Path(sys.argv[1] if len(sys.argv) > 1
              else "recipes/refine/futtsu_nori.yaml").resolve()
OUT = Path(os.environ.get("LR_OUT", f"outputs/refine_{recipe.stem}")).resolve()
OUT.mkdir(parents=True, exist_ok=True)
LAND = Path(os.environ.get("FMESH_LAND",
                           "outputs/sample_repro/land_channel_adj.shp")).resolve()


def say(msg):
    print(f"[lr] {msg} +{time.time() - t0:.0f}s", flush=True)


cfg = load_refine(recipe)
base = read_fort14(cfg["base_mesh"])
say(f"recipe {recipe.name}: base {cfg['base_mesh'].name} "
    f"NP={base.n_nodes:,} NE={base.n_elements:,}")

land_m = unary_union(list(gpd.read_file(LAND).to_crs(MESH_EPSG).geometry))
say(f"land {LAND.name}: {len(getattr(land_m, 'geoms', [land_m]))} polygons")

# ---------------------------------------------------------------- geometry
# The recipe declares the region in lon/lat; everything below is in the mesh
# CRS (metres), which is also what DistMesh is handed.  Working in metres is
# what disables generate_mesh's internal tmerc sandwich -- it only engages for
# a bbox that looks like degrees -- so fd/fh are evaluated in the frame they
# were written in.
from pyproj import Transformer  # noqa: E402

to_m = Transformer.from_crs("EPSG:4326", f"EPSG:{MESH_EPSG}", always_xy=True)


def region_in_metres(region):
    """Project a declared region to the mesh CRS.

    A circle is re-struck about the projected centre rather than projected
    vertex by vertex: the recipe's radius is a distance on the ground, and a
    300 m disc drawn in degrees and then projected is an ellipse.
    """
    g = region.geometry
    c = g.centroid
    cx, cy = to_m.transform(c.x, c.y)
    ring = np.asarray(g.exterior.coords)
    x, y = to_m.transform(ring[:, 0], ring[:, 1])
    r = np.hypot(x - cx, y - cy)
    if float(r.std() / r.mean()) < 0.02:
        return shapely.Point(cx, cy).buffer(float(r.mean()), quad_segs=64)
    return shapely.Polygon(np.column_stack([x, y]))


regions_m = [(region_in_metres(r), r) for r in cfg["refine"]]

# The ambient size is MEASURED on the base mesh around the site, not taken
# from the sizing recipe: what the transition has to reach is the mesh that
# is actually there.  Near Futtsu that is 420 m, not the 350 m nominal.
amb_node = ambient_size_field(base.nodes, base.elements)
ambient = {}
for geom, region in regions_m:
    d = shapely.distance(shapely.points(base.nodes[:, 0], base.nodes[:, 1]), geom)
    near = d < 3000.0
    ambient[region.name] = float(np.median(amb_node[near]))
    say(f"ambient around {region.name}: {ambient[region.name]:.0f} m "
        f"(p90 {np.percentile(amb_node[near], 90):.0f} m)")

# ------------------------------------------------------------- pre-flight
land_ll = unary_union(list(gpd.read_file(LAND).to_crs(4326).geometry))
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402

_mt = Triangulation(base.nodes[:, 0], base.nodes[:, 1], base.elements)
_di = LinearTriInterpolator(_mt, base.depths)


def depth_of(lon, lat):
    x, y = to_m.transform(np.asarray(lon), np.asarray(lat))
    # Outside the base mesh there is no depth to inherit, so NaN is the
    # honest answer; preflight refuses it, which is right for a region that
    # is not over this mesh at all.
    return np.ma.filled(_di(x, y), np.nan)


reports = {"recipe": str(recipe), "base_mesh": str(cfg["base_mesh"]),
           "preflight": [], "land": str(LAND)}
for geom, region in regions_m:
    pf = preflight(region, gradation=cfg["gradation"],
                   dt_expected_s=cfg["dt_expected_s"],
                   ambient_h_m=ambient[region.name], depth_of=depth_of,
                   land=land_ll)
    reports["preflight"].append(pf)
    say(f"preflight {region.name}: transition {pf['transition_m']:.0f} m, "
        f"dt {pf['dt_s']:.2f} s, +{pf['elements_added']:.0f} elements")
    if pf["dt_alert"]:
        say("ALERT " + pf["dt_alert"])
    if pf["land_alert"]:
        say("ALERT " + pf["land_alert"])

# ------------------------------------------------------------------- cut
widths = {region.name: (region.transition_m if region.transition_m is not None
                        else transition_width_m(region.target_h_m,
                                                ambient[region.name],
                                                cfg["gradation"]))
          for _, region in regions_m}
footprint = unary_union([geom.buffer(widths[region.name])
                         for geom, region in regions_m])
grad = effective_gradation(base.nodes, base.elements,
                           [(g, r.target_h_m, widths[r.name]) for g, r in regions_m])
reports["gradation"] = grad
say(f"effective gradation: max {grad['max_effective_gradation']:.3f} against the "
    f"{grad['c4_limit_gradation']:.3f} that C4 allows "
    f"({'within' if grad['within_c4'] else 'OVER'})")
obc_nodes = np.concatenate([np.asarray(s) for s in base.open_boundaries]) \
    if base.open_boundaries else np.empty(0, dtype=np.int64)
sel = select_patch(base.nodes, base.elements, footprint,
                   open_boundary_nodes=obc_nodes,
                   obc_guard_m=float(os.environ.get("LR_OBC_GUARD", 500.0)))
reports["selection"] = sel.report
say("cut: " + json.dumps(sel.report))

# ------------------------------------------------------------------- rim
# The source shoreline for `resample` is picked ONCE per stretch, by the land
# ring the stretch already lies on.  Picking it inside the resampler by
# nearest distance would let a stretch jump to the opposite bank at a strait.
free = np.setdiff1d(np.unique(sel.rim_edges), sel.frozen_nodes)
shore = None
if cfg["coastline"] == "resample" and free.size:
    pts = shapely.MultiPoint(base.nodes[free, :2])
    ringlist = [g.exterior for g in getattr(land_m, "geoms", [land_m])]
    shore = min(ringlist, key=lambda r: shapely.distance(pts, r))
    shore = shapely.LineString(np.asarray(shore.coords))
    say(f"source shoreline: ring of {len(shore.coords):,} points, "
        f"{shapely.distance(pts, shore):.1f} m from the free rim")

target = min(r.target_h_m for _, r in regions_m)
# The coastline inside the hole is cut at the LOCAL size, not at the target.
# h_achieved is the same field DistMesh gets, without the 1.2 field-to-bar
# factor, because these are the bar lengths themselves.
h_achieved = patch_sizing(base.nodes, base.elements,
                          [(g, r.target_h_m, widths[r.name]) for g, r in regions_m],
                          distmesh_scale=1.0)
rc = rim_constraints(base.nodes, sel, size=h_achieved,
                     coastline=cfg["coastline"], shoreline=shore,
                     tolerance_m=cfg["coastline_tolerance_m"])
reports["rim"] = {k: v for k, v in rc.items()
                  if isinstance(v, (int, float, str, bool))}
say("rim: " + json.dumps(reports["rim"]))

hole = hole_polygon(rc["pfix"], rc["egfix"])
say(f"hole {hole.area / 1e6:.3f} km2 ({hole.geom_type}, valid={hole.is_valid})")

# ------------------------------------------------------------------ fill
fh = patch_sizing(base.nodes, base.elements,
                  [(g, r.target_h_m, widths[r.name]) for g, r in regions_m],
                  distmesh_scale=DISTMESH_SCALE)
shapely.prepare(hole)
boundary = shapely.boundary(hole)


def fd(points):
    """Signed distance to the hole: negative inside, in metres."""
    p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
    pt = shapely.points(p[:, 0], p[:, 1])
    d = shapely.distance(pt, boundary)
    return np.where(shapely.contains(hole, pt), -d, d)


xmin, ymin, xmax, ymax = hole.bounds
bbox = (float(xmin), float(xmax), float(ymin), float(ymax))
# The field's floor is the finest target, in FIELD space.  Sampling for it
# would be a lottery: a 300 m core is under 1 % of the hole's bounding box,
# and a seeding lattice anchored at a too-large hmin under-seeds the core.
hmin = target / DISTMESH_SCALE
say(f"fill: bbox {bbox[1] - bbox[0]:.0f} x {bbox[3] - bbox[2]:.0f} m, "
    f"hmin {hmin:.1f} m, pfix {rc['n_pfix']}, egfix {rc['n_egfix']}")

# cleanup='none', and the safe stages run by hand below.  The default clean
# ends in make_mesh_boundaries_traversable, which is the one stage that takes
# no pfix and no egfix: on this patch it deleted 121 of the 235 constrained
# rim points (probe, 2026-09-22).  That is defensible for a mesh whose
# boundary is an output and fatal for one whose boundary is the contract.
p, t = om.generate_mesh(
    fd, fh, bbox=bbox, min_edge_length=hmin,
    max_iter=int(os.environ.get("LR_MAX_ITER", 100)),
    seed=int(os.environ.get("LR_SEED", 0)),
    pfix=rc["pfix"], egfix=rc["egfix"], cleanup="none")
say(f"filled: NP={len(p):,} NE={len(t):,}")

from oceanmesh.mesh_improve import (  # noqa: E402
    collapse_thin_triangles,
    direct_smoother_lur,
)

p, t = collapse_thin_triangles(p, t, min_qual=0.25, pfix=rc["pfix"])
p, t = direct_smoother_lur(p, t, pfix=rc["pfix"])
say(f"cleaned (pfix-protected): NP={len(p):,} NE={len(t):,}")
# bound_connectivity is NOT run: its valence flips are blind to the sizing
# field and on this patch they coarsened the 28.6 m core to 98.7 m.

# Faces whose centroid is outside the hole are the CDT's convex-hull fill;
# they are not part of the patch.  Everything the rim constrains survives
# because stitch_patch checks each fixed point individually afterwards.
cen = p[t].mean(axis=1)
keep = np.asarray(shapely.contains(hole, shapely.points(cen[:, 0], cen[:, 1])))
say(f"outside-hole faces dropped: {int((~keep).sum()):,}")
t = t[keep]
used = np.unique(t)
remap = np.full(len(p), -1, dtype=np.int64)
remap[used] = np.arange(len(used))
p, t = p[used], remap[t]

# ---------------------------------------------------------------- stitch
nodes, elements, depths, node_map, st = stitch_patch(
    base.nodes, base.elements, base.depths, sel, p, t,
    rc["pfix"], rc["pfix_base"])
reports["stitch"] = st
say("stitch: " + json.dumps(st))

# Orientation: fort.14 wants counter-clockwise.  The retained faces already
# are, so only the patch can be wrong, and flipping the whole array would
# break the frozen-connectivity check that runs next.
a = nodes[elements[:, 1]] - nodes[elements[:, 0]]
b = nodes[elements[:, 2]] - nodes[elements[:, 0]]
cw = (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) <= 0
if cw.any():
    elements[cw] = elements[cw][:, [0, 2, 1]]
    say(f"reoriented {int(cw.sum()):,} clockwise elements")

# ------------------------------------------------------------ seam repair
# The base mesh was finished to sit exactly on its gates -- min angle 30.01
# deg, area change 0.500 -- so it has no margin to absorb a new neighbour,
# and the raw stitch came out at 26.1 deg with seven area jumps over the
# limit.  Every offender was a PATCH element, so they are repairable in
# place: improve_patch may flip only patch faces and move only new interior
# nodes, which is what keeps the frozen zone frozen through the repair.
_u, _c = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                      elements[:, [2, 0]]]), axis=1),
                   axis=0, return_counts=True)
on_boundary = np.zeros(len(nodes), dtype=bool)
on_boundary[np.unique(_u[_c == 1])] = True
is_new = np.arange(len(nodes)) >= st["n_nodes_retained"]
movable = is_new & ~on_boundary
slidable = is_new & on_boundary
mutable_faces = np.arange(len(elements)) >= len(sel.retained)
_before_repair = nodes.copy()
# The curves new boundary nodes may slide along are the ones rim_constraints
# cut them from -- one per replaced stretch -- and nothing else.  The base
# mesh's own land_boundaries list will not do: only 820 of its 939
# consecutive pairs are boundary edges and it jumps up to 2,867 m, so a line
# built from it runs through open water and a node projected onto it lands
# in the sea.
slide_on = [shapely.LineString(c) for c in rc["curves"] if len(c) > 1]
nodes, elements, imp = improve_patch(
    nodes, elements, movable, mutable_faces, slidable=slidable,
    slide_on=slide_on, only_below=float(os.environ.get("LR_ONLY_BELOW", 1.15)))
# How faithful the patch's coastline is, measured against the curve it was
# cut from rather than against itself.  Both the nodes and the line between
# them are checked: a node is kept on the curve by construction, but the
# chords between nodes are the coastline the model will actually see, and
# that is what coastline_tolerance_m is a statement about.
_u2, _c2 = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                        elements[:, [2, 0]]]), axis=1),
                     axis=0, return_counts=True)
_bnd = _u2[_c2 == 1]
_new_bnd = _bnd[is_new[_bnd].any(axis=1)]
if len(_new_bnd):
    _a, _b = nodes[_new_bnd[:, 0]], nodes[_new_bnd[:, 1]]
    _f = np.linspace(0.0, 1.0, 9)[:, None, None]
    _samp = (_a[None] + _f * (_b - _a)[None]).reshape(-1, 2)
    _curve = shapely.MultiLineString([np.asarray(ln.coords) for ln in slide_on]) \
        if slide_on else None
    imp["coastline_departure_m"] = float(shapely.distance(
        shapely.points(_samp), _curve).max()) if _curve else 0.0
    # NEW nodes only: the frozen anchors at each end of a stretch sit on the
    # base polyline, which is exactly what `resample` departs from, so
    # including them would report the improvement as an error.
    _nb = np.unique(_new_bnd)
    _nb = _nb[is_new[_nb]]
    imp["coastline_node_departure_m"] = float(shapely.distance(
        shapely.points(nodes[_nb]), _curve).max()) if _curve and len(_nb) else 0.0
else:
    imp["coastline_departure_m"] = 0.0
    imp["coastline_node_departure_m"] = 0.0
reports["improve"] = imp
if imp["coastline_departure_m"] > cfg["coastline_tolerance_m"]:
    raise SystemExit(
        f"the repair moved the coastline {imp['coastline_departure_m']:.0f} m "
        f"from the curve it was cut from, past the "
        f"{cfg['coastline_tolerance_m']:g} m tolerance; a slide may cross an "
        "original vertex and chord off the bend behind it")
say(f"seam repair: {imp['n_flips']} flips, {imp['n_moves']} moves, "
    f"coastline departure {imp['coastline_departure_m']:.1f} m "
    f"(nodes {imp['coastline_node_departure_m']:.2f} m), "
    f"angles {imp['min_angle_deg']:.2f}-{imp['max_angle_deg']:.2f} deg "
    f"({int(movable.sum()):,} movable, {int(slidable.sum()):,} slidable nodes, "
    f"{int(mutable_faces.sum()):,} mutable faces)")

# ---------------------------------------------------------------- verify
ver = verify_patch(base.nodes, base.depths, base.elements, sel,
                   nodes, elements, depths, node_map,
                   open_boundaries=base.open_boundaries)
reports["verify"] = ver
say("verify: " + json.dumps(ver))
if not ver["ok"]:
    raise SystemExit("the frozen zone is not frozen; see verify in the report")

# ------------------------------------------------------------- boundaries
# Rebuilt rather than carried over: the patch changes the boundary node list
# wherever the coastline was re-cut.  The OBC is the exception -- it is an
# input, the cut is forbidden to touch it, and verify_patch has just checked
# that every one of its nodes is where it was -- so it is mapped, not found.
obc = [node_map[np.asarray(s, dtype=np.int64)] for s in base.open_boundaries]
loops = om.boundary_loops(elements)
obc_set = set(np.concatenate(obc).tolist()) if obc else set()
outer = max(loops, key=lambda lp: abs(shapely.Polygon(nodes[lp]).area))
if not obc_set.issubset(set(outer.tolist())):
    raise SystemExit("the open boundary is no longer on the outer loop")
# Split the outer loop at the OBC's two ends: the run between them that stays
# on the OBC is the open string, the complement is the mainland coast.
ring = np.roll(outer, -int(np.where(outer == obc[0][0])[0][0]))
if ring[1] not in obc_set:
    ring = np.roll(ring[::-1], 1)
stop = int(np.where(ring == obc[0][-1])[0][0])
if set(ring[:stop + 1].tolist()) != obc_set:
    raise SystemExit("the outer loop between the OBC ends is not the OBC")
land_bounds = [(20, np.append(ring[stop:], ring[0]))]
land_bounds += [(21, lp) for lp in loops if lp is not outer]
say(f"boundaries: open {len(obc[0])}, mainland {len(land_bounds[0][1])}, "
    f"islands {len(land_bounds) - 1}")
mesh = Fort14Mesh(
    title=f"{base.title} + {recipe.stem}",
    nodes=nodes, depths=depths, elements=elements,
    open_boundaries=[ring[:stop + 1]],
    land_boundaries=land_bounds)
out14 = OUT / f"{Path(cfg['base_mesh']).stem}_{recipe.stem}.14"
write_fort14(mesh, out14)
say(f"wrote {out14.name}: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}")
np.save(OUT / "node_map.npy", node_map)
(OUT / "report.json").write_text(json.dumps(reports, indent=1, default=float))
say("done")
