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
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.patch import (
    ambient_size_field,
    boundary_after_patch,
    effective_gradation,
    field_gradation,
    hole_polygon,
    improve_patch,
    introduced_violations,
    patch_sizing,
    refresh_depths,
    region_conflicts,
    rim_constraints,
    select_patch,
    stitch_patch,
    verify_patch,
)
from fvcom_mesh_tools.qa import run_qa
from fvcom_mesh_tools.refine import (
    limit_rfactor,
    load_refine,
    preflight,
    transition_width_m,
)

# The base's CRS. fort.14 and FVCOM's _grd.dat both carry bare numbers, so
# this is an assumption, not something read from the file; goto2023 is UTM
# 54N and so is everything downstream of it here.
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
# The base is a FINISHED case: an FVCOM `_grd.dat` with the depth file the
# baseline names, or a fort.14 that carries its own. Either way nothing here
# rebuilds it -- the depths are the model's and the refinement inherits them
# (owner, 2026-09-22).
if cfg["base_mesh"].suffix == ".dat":
    base = read_fvcom_case(cfg["base_mesh"], cfg["base_depth"], cfg["base_obc"])
    say(f"recipe {recipe.name}: base {cfg['base_mesh'].name} + "
        f"{cfg['base_depth'].name}")
else:
    base = read_fort14(cfg["base_mesh"])
    say(f"recipe {recipe.name}: base {cfg['base_mesh'].name}")
_be = np.unique(np.sort(np.vstack([base.elements[:, [0, 1]], base.elements[:, [1, 2]],
                                   base.elements[:, [2, 0]]]), axis=1), axis=0)
base_rmax = float((np.abs(base.depths[_be[:, 0]] - base.depths[_be[:, 1]])
                   / (base.depths[_be[:, 0]] + base.depths[_be[:, 1]])).max())
say(f"  NP={base.n_nodes:,} NE={base.n_elements:,}, depth "
    f"{base.depths.min():.3f}-{base.depths.max():.3f} m, r-factor <= "
    f"{base_rmax:.4f}, {len(base.open_boundaries)} open boundary")

# The land polygon is optional and only two things use it: the pre-flight's
# "is the core dry" report, and `coastline: resample`. `preserve` needs
# neither, and the base mesh is a better dryness test anyway -- a sample
# outside it has no depth to inherit, which pre-flight already refuses.
land_m = None
if LAND.exists():
    land_m = unary_union(list(gpd.read_file(LAND).to_crs(MESH_EPSG).geometry))
    say(f"land {LAND.name}: {len(getattr(land_m, 'geoms', [land_m]))} polygons")
elif cfg["coastline"] == "resample":
    raise SystemExit(f"coastline: resample needs a source shoreline; "
                     f"FMESH_LAND={LAND} does not exist")

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
    300 m disc drawn in degrees and then projected is an ellipse.  Which
    regions are circles comes from what the recipe DECLARED, not from the
    spread of vertex radii -- a square's corners are all equidistant from its
    centre, so that guess turned a bbox into a 257-point disc.

    Every ring is projected, interiors included: a fishery boundary read from
    GeoJSON may have holes, and the parser already accepts them.
    """
    if region.kind == "circle":
        lon, lat, radius_m = region.circle
        cx, cy = to_m.transform(lon, lat)
        return shapely.Point(cx, cy).buffer(radius_m, quad_segs=64)

    def ring(coords):
        arr = np.asarray(coords)
        x, y = to_m.transform(arr[:, 0], arr[:, 1])
        return np.column_stack([x, y])

    g = region.geometry
    return shapely.Polygon(ring(g.exterior.coords),
                           [ring(h.coords) for h in g.interiors])


regions_m = [(region_in_metres(r), r) for r in cfg["refine"]]
for _g, _r in regions_m:
    say(f"region {_r.name}: {_r.kind}, {_g.area / 1e6:.4f} km2, "
        f"{len(_g.exterior.coords)} vertices"
        + (f", from {Path(_r.source['file']).name}" if _r.source else ""))

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
land_ll = (unary_union(list(gpd.read_file(LAND).to_crs(4326).geometry))
           if land_m is not None else None)
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402

_mt = Triangulation(base.nodes[:, 0], base.nodes[:, 1], base.elements)
_di = LinearTriInterpolator(_mt, base.depths)


def depth_of(lon, lat):
    x, y = to_m.transform(np.asarray(lon), np.asarray(lat))
    # Outside the base mesh there is no depth to inherit, so NaN is the
    # honest answer; preflight refuses it, which is right for a region that
    # is not over this mesh at all.
    return np.ma.filled(_di(x, y), np.nan)


# The declared regions in the mesh CRS, so a figure can outline what was
# asked for without re-deriving it from a centre and a radius that only a
# circle has.
regions_report = [
    {"name": r.name, "kind": r.kind, "target_h_m": r.target_h_m,
     "source": r.source, "area_m2": float(g.area),
     "xy": np.asarray(g.exterior.coords).round(3).tolist()}
    for g, r in regions_m]
reports = {"recipe": str(recipe), "base_mesh": str(cfg["base_mesh"]),
           "regions": regions_report,
           "base_depth": str(cfg["base_depth"]) if cfg["base_depth"] else None,
           "base_obc": str(cfg["base_obc"]) if cfg["base_obc"] else None,
           "base_rmax": base_rmax, "preflight": [],
           "land": str(LAND) if land_m is not None else None}
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
sized = [(g, r.target_h_m, widths[r.name], r.priority) for g, r in regions_m]
grad = effective_gradation(base.nodes, base.elements, sized)
# Overlapping fisheries are an ordinary input; what is not ordinary is a
# region silently getting a size it did not ask for, so it is measured.
conflicts = region_conflicts(sized, [r.name for _, r in regions_m])
reports["conflicts"] = conflicts
if conflicts["any_overlap"]:
    say("regions overlap; a target is a ceiling, so the shared water takes "
        "the smallest of them")
    for pair in conflicts["overlapping_pairs"]:
        say(f"  {pair['regions'][0]} & {pair['regions'][1]}: "
            f"{pair['overlap_m2'] / 1e6:.4f} km2 shared, "
            f"{pair['effective_target_h_m']:g} m applies")
    for name, f in conflicts["finer_than_declared"].items():
        say(f"  {name}: {100 * f['fraction']:.1f} % of its area "
            f"({f['area_m2'] / 1e6:.4f} km2) comes out at {f['gets_h_m']:g} m "
            f"rather than its own {f['own_target_h_m']:g} m, because of "
            f"{', '.join(f['because_of'])} -- finer than asked, and somebody "
            "pays for it in elements and time step")
    if conflicts["priority_ignored"]:
        say("  NOTE priority differs between overlapping regions and has no "
            "effect here: in a refinement a target is a ceiling, so the "
            "finest always applies. Priority coarsens only in a sizing "
            "recipe, where the whole mesh is rebuilt")
reports["gradation"] = grad
say(f"effective gradation: max {grad['max_effective_gradation']:.3f} against a "
    f"{grad['c4_reference_gradation']:.3f} reference "
    f"({'below' if grad['ramp_below_reference'] else 'ABOVE'}); C4 itself is "
    "gated on the finished mesh")
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
    # Every land ring within reach is offered, and rim_constraints picks the
    # one nearest EACH stretch.  Offering a single nearest ring for the whole
    # free rim is a closest-pair distance, and two stretches on opposite
    # banks of a strait would both be resampled onto whichever bank won.
    pts = shapely.MultiPoint(base.nodes[free, :2])
    reach = pts.buffer(3000.0)
    shore = []
    for g in getattr(land_m, "geoms", [land_m]):
        for r in [g.exterior, *g.interiors]:
            if shapely.intersects(reach, r):
                shore.append(shapely.LineString(np.asarray(r.coords)))
    if not shore:
        raise SystemExit("coastline: resample found no source shoreline within "
                         "3 km of the free rim")
    say(f"source shoreline: {len(shore)} ring(s) within 3 km, nearest "
        f"{min(shapely.distance(pts, ln) for ln in shore):.1f} m from the free rim")

target = min(r.target_h_m for _, r in regions_m)
# The coastline inside the hole is cut at the LOCAL size, not at the target.
# h_achieved is the same field DistMesh gets, without the 1.2 field-to-bar
# factor, because these are the bar lengths themselves.
h_achieved = patch_sizing(base.nodes, base.elements, sized, distmesh_scale=1.0)
rc = rim_constraints(base.nodes, sel, size=h_achieved,
                     coastline=cfg["coastline"], shoreline=shore,
                     tolerance_m=cfg["coastline_tolerance_m"])
reports["rim"] = {k: v for k, v in rc.items()
                  if isinstance(v, (int, float, str, bool))}
say("rim: " + json.dumps(reports["rim"]))

hole = hole_polygon(rc["pfix"], rc["egfix"])
say(f"hole {hole.area / 1e6:.3f} km2 ({hole.geom_type}, valid={hole.is_valid})")

# The slope the field actually has, measured on the hole it will be meshed
# in -- not the per-region formula, which omits the ambient term and says
# nothing about where two regions meet.
fslope = field_gradation(h_achieved, hole)
reports["field_gradation"] = fslope
say(f"field slope: max {fslope['max_slope']:.3f}, p99 {fslope['p99_slope']:.3f} "
    f"against a {fslope['c4_reference_gradation']:.3f} reference; "
    f"{100 * fslope['fraction_above_reference']:.2f} % of samples above it")
if fslope["max_slope"] > fslope["c4_reference_gradation"]:
    say("  ALERT the field is locally steeper than two similar triangles one "
        "size apart can be and still pass C4. Expect the fill to be hard: "
        "widen a transition, move a region, or bring the targets closer")

# ------------------------------------------------------------------ fill
fh = patch_sizing(base.nodes, base.elements, sized,
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

from oceanmesh.mesh_improve import (  # noqa: E402
    collapse_thin_triangles,
    direct_smoother_lur,
)


def attempt(seed):
    """One fill-repair-stitch-verify pass at a given DistMesh seed.

    Returns ``(candidate, report)``; ``candidate`` is None when this seed did
    not produce a mesh that keeps the contract.  The seed is a real knob and
    not a cosmetic one: the repair is greedy, so it stops at a local optimum
    whose quality depends on where the fill started.  Measured on this patch,
    neighbouring configurations finished at 27.5 deg and at 30.01 deg.  The
    driver therefore searches seeds and reports which one it used, rather
    than reporting whichever one it happened to try.
    """
    out = {"seed": seed}
    # cleanup='none', and the safe stages run by hand below.  The default clean
    # ends in make_mesh_boundaries_traversable, which is the one stage that takes
    # no pfix and no egfix: on this patch it deleted 121 of the 235 constrained
    # rim points (probe, 2026-09-22).  That is defensible for a mesh whose
    # boundary is an output and fatal for one whose boundary is the contract.
    p, t = om.generate_mesh(
        fd, fh, bbox=bbox, min_edge_length=hmin,
        max_iter=int(os.environ.get("LR_MAX_ITER", 100)),
        seed=seed,
        pfix=rc["pfix"], egfix=rc["egfix"], cleanup="none")
    say(f"filled: NP={len(p):,} NE={len(t):,}")

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
    out["stitch"] = {k: v for k, v in st.items() if k != "pfix_new"}
    say("stitch: " + json.dumps(out["stitch"]))
    # Before the repair, while a fixed point is still exactly where it was
    # put: the repair slides boundary nodes along their curve, so matching
    # pfix by coordinate afterwards fails.  Node ids do not change in the
    # repair, so the edge set built here stays valid.
    want_boundary = boundary_after_patch(base.elements, sel, rc, node_map,
                                         st["pfix_new"])
    out["want_boundary"] = want_boundary

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

    # Depths follow the nodes.  stitch_patch evaluated the base field at the
    # positions the fill produced; improve_patch then moved some of those nodes,
    # and a depth left behind at the old position is not the base field at the
    # delivered coordinate -- on a 5+x field a node moved from (0.2,0.3) to (1,1)
    # kept 5.2 where 6.0 is right (review finding 9, 2026-09-22).  The frozen
    # depths are untouched because frozen nodes do not move.
    only_below = float(os.environ.get("LR_ONLY_BELOW", 1.15))
    nodes, elements, imp = improve_patch(
        nodes, elements, movable, mutable_faces, slidable=slidable,
        slide_on=slide_on, only_below=only_below)
    if imp["min_angle_deg"] < 30.0 or imp["max_angle_deg"] > 130.0:
        # Strict improvement has stalled.  The soft pass allows a move that
        # holds the worst margin and improves the rest, which is how a run
        # stuck at 22.96 deg got to 27.46; it is tried second because when
        # strict succeeds it succeeds better.
        nodes, elements, imp2 = improve_patch(
            nodes, elements, movable, mutable_faces, slidable=slidable,
            slide_on=slide_on, only_below=only_below, soft=True)
        imp = {**imp2, "n_flips": imp["n_flips"] + imp2["n_flips"],
               "n_moves": imp["n_moves"] + imp2["n_moves"], "soft_pass": True}
    _shifted = is_new & (np.linalg.norm(nodes - _before_repair, axis=1) > 0)
    depths, _n_out = refresh_depths(base.nodes, base.elements, base.depths,
                                    nodes, _shifted, depths)
    imp["n_depths_recomputed"] = int(_shifted.sum())
    imp["n_recomputed_outside_base"] = int(_n_out)

    # Inherit the base's r-factor property, not just its values. The base is
    # m7001tp_rfac0p2_cap300: every one of its edges satisfies r <= 0.2, and
    # interpolation does not carry that across a new edge joining different
    # base elements -- ten new edges came out above it, the worst at 0.3075.
    # Only new nodes' depths move; every base depth is untouched, and the
    # report says how far a new one was pulled.
    if cfg["rfactor_limit"] != "off":
        rmax = base_rmax if cfg["rfactor_limit"] == "base" \
            else float(cfg["rfactor_limit"])
        depths, rinfo = limit_rfactor(
            elements, depths, is_new, rmax,
            depth_min=float(base.depths.min()), depth_max=float(base.depths.max()))
        # A frozen-pair edge over the limit is the base's only if the base
        # has that edge. New connectivity can join two retained nodes that
        # were never neighbours, and that edge is the patch's.
        _b_edges = {tuple(sorted(x)) for x in
                    np.unique(np.sort(np.vstack(
                        [base.elements[:, [0, 1]], base.elements[:, [1, 2]],
                         base.elements[:, [2, 0]]]), axis=1),
                        axis=0).tolist()}
        _inv = np.full(len(nodes), -1, dtype=np.int64)
        _inv[node_map[node_map >= 0]] = np.flatnonzero(node_map >= 0)
        _new_over = [ab for ab in rinfo["frozen_pair_edges_over_rmax"]
                     if tuple(sorted(_inv[list(ab)].tolist())) not in _b_edges]
        rinfo["n_new_frozen_pair_over_rmax"] = len(_new_over)
        rinfo.pop("frozen_pair_edges_over_rmax", None)
        out["rfactor"] = rinfo
        say(f"r-factor <= {rmax:.4f}: {rinfo['n_depths_changed']} new depths "
            f"moved, worst {rinfo['max_depth_change_m']:.2f} m, "
            f"{'converged' if rinfo['converged'] else 'NOT CONVERGED'} in "
            f"{rinfo['rounds']} rounds (base depths moved "
            f"{rinfo['max_frozen_depth_change_m']:.3g} m)")
        if rinfo["n_over_rmax_movable"] or rinfo["n_new_frozen_pair_over_rmax"]:
            # The base guarantees r <= rmax and the patch is supposed to
            # inherit it. Failing to is this seed's problem, not something to
            # print and carry on from.
            say(f"    the patch leaves {rinfo['n_over_rmax_movable']} movable "
                f"and {rinfo['n_new_frozen_pair_over_rmax']} new frozen-pair "
                f"edge(s) above r = {rmax:.4f}")
            return None, out
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
    out["improve"] = imp
    if imp["coastline_departure_m"] > cfg["coastline_tolerance_m"]:
        say(f"  seed {seed}: the repair moved the coastline "
            f"{imp['coastline_departure_m']:.0f} m from the curve it was cut "
            f"from, past the {cfg['coastline_tolerance_m']:g} m tolerance")
        return None, out
    say(f"seam repair: {imp['n_flips']} flips, {imp['n_moves']} moves, "
        f"coastline departure {imp['coastline_departure_m']:.1f} m "
        f"(nodes {imp['coastline_node_departure_m']:.2f} m), "
        f"angles {imp['min_angle_deg']:.2f}-{imp['max_angle_deg']:.2f} deg "
        f"({int(movable.sum()):,} movable, {int(slidable.sum()):,} slidable nodes, "
        f"{int(mutable_faces.sum()):,} mutable faces)")

    # ---------------------------------------------------------------- verify
    ver = verify_patch(base.nodes, base.depths, base.elements, sel,
                       nodes, elements, depths, node_map,
                       open_boundaries=base.open_boundaries,
                       expected_boundary=want_boundary)
    out["verify"] = ver
    say("verify: " + json.dumps(ver))
    if not ver["ok"]:
        return None, out
    return (nodes, elements, depths, node_map), out


def serialise(candidate, out, path):
    """Build the boundary lists, write the fort.14, and read it back.

    Everything before this checked objects in memory.  What the model runs is
    the file, so the file is what is verified: the frozen zone again, through
    the map, and then the whole 21-gate battery by the caller.
    """
    nodes, elements, depths, node_map = candidate
    # ------------------------------------------------------------- boundaries
    # Rebuilt rather than carried over: the patch changes the boundary node list
    # wherever the coastline was re-cut.  The OBC is the exception -- it is an
    # input, the cut is forbidden to touch it, and verify_patch has just checked
    # that every one of its nodes is where it was -- so it is mapped, not found.
    if len(base.open_boundaries) != 1:
        # The split below assumes one open arc on the outer loop: it indexes
        # obc[0] and tests a set.  Two disjoint arcs would silently merge and
        # none would raise an IndexError instead.  Refuse rather than guess.
        raise SystemExit(
            f"this generator handles exactly one open boundary; the base mesh has "
            f"{len(base.open_boundaries)}")
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
    write_fort14(mesh, path)
    say(f"wrote {path.name}: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}")

    # --------------------------------------------------- the delivered artefact
    # Everything above checked objects in memory.  What the model runs is the
    # file, so the file is what is checked: read it back, verify the frozen zone
    # through the map again, and run the whole 21-gate battery.  A generator
    # whose quality claim rests on someone remembering to run QA afterwards does
    # not have a quality claim (review findings 3 and 15, 2026-09-22).
    written = read_fort14(path)
    out["verify_on_disk"] = verify_patch(
        base.nodes, base.depths, base.elements, sel,
        written.nodes, written.elements, written.depths, node_map,
        open_boundaries=base.open_boundaries,
        expected_boundary=out["want_boundary"])
    if not out["verify_on_disk"]["ok"]:
        return None, None
    if not np.array_equal(written.open_boundaries[0], mesh.open_boundaries[0]) or \
            len(written.land_boundaries) != len(mesh.land_boundaries):
        raise SystemExit("the written boundary lists differ from the ones built")
    return written, mesh


# --------------------------------------------------------- the seed search
out14 = OUT / f"{Path(cfg['base_mesh']).stem}_{recipe.stem}.14"
seeds = [int(x) for x in os.environ.get("LR_SEEDS", "0,1,2,3,4").split(",")]
best = None
reports["attempts"] = []
def save_report():
    (OUT / "report.json").write_text(json.dumps(reports, indent=1, default=float))


for seed in seeds:
    say(f"--- seed {seed}")
    try:
        candidate, out = attempt(seed)
    except ValueError as exc:
        # A candidate that cannot be built is one seed's problem, not the
        # search's: stitch_patch raises when the fill loses a constrained
        # point, and a bad first seed used to end the whole run (second
        # review, finding 5).  Anything that is not a candidate failure --
        # a bad recipe, a coding error -- still propagates.
        say(f"    seed {seed}: {exc}")
        reports["attempts"].append({"seed": seed, "error": str(exc)})
        save_report()
        continue
    if candidate is None:
        reports["attempts"].append({k: v for k, v in out.items()
                                    if k != "want_boundary"})
        save_report()
        continue
    written, mesh = serialise(candidate, out, out14)
    if written is None:
        reports["attempts"].append({k: v for k, v in out.items()
                                    if k != "want_boundary"})
        save_report()
        continue
    qa = run_qa(written, name=out14.stem, path=out14, max_offenders=10_000)
    # What the patch is answerable for. A refinement may not be held to a
    # standard its base does not meet: the goto2023 production mesh fails C1
    # at one element 18 km from Futtsu, the contract freezes that element,
    # and an absolute gate blamed every seed for it.
    new_bad = introduced_violations(qa.checks, len(sel.retained), written.elements)
    out["qa"] = {"n_gate_total": qa.n_gate_total,
                 "n_gate_failed": qa.n_gate_failed,
                 "n_introduced": len(new_bad),
                 "introduced": new_bad[:20],
                 "failed": [{"check": c.check_id, "requirement": c.requirement,
                             "observed": c.observed} for c in qa.checks
                            if c.status == "fail"]}
    reports["attempts"].append({k: v for k, v in out.items()
                                if k != "want_boundary"})
    save_report()
    say(f"    seed {seed}: QA {qa.n_gate_total - qa.n_gate_failed}/"
        f"{qa.n_gate_total}, {len(new_bad)} introduced by the patch"
        + ("" if not new_bad else "  " + "; ".join(
            f"{v['check']} at {v['kind']} {v.get('id', v.get('elements'))}"
            for v in new_bad[:4])))
    if best is None or len(new_bad) < best[0]:
        best = (len(new_bad), seed, candidate, out, written, qa)
    if not new_bad:
        break

if best is None:
    save_report()
    raise SystemExit("no seed produced a mesh that keeps the frozen-zone "
                     f"contract; see attempts in {OUT / 'report.json'}")
_, seed, candidate, out, written, qa = best
nodes, elements, depths, node_map = candidate
reports.update({k: v for k, v in out.items()
                if k not in ("seed", "want_boundary")})
reports["seed"] = seed
np.save(OUT / "node_map.npy", node_map)
reports["mesh"] = str(out14)
# Always, not only when the accepted seed is not the first one tried: the
# file on disk is whichever attempt ran last, and when none passed that is
# not the best one.  A failed artefact investigated with another attempt's
# node map and QA is worse than no artefact (second review, finding 4).
written, mesh = serialise(candidate, out, out14)
qa = run_qa(written, name=out14.stem, path=out14, max_offenders=10_000)
say(f"accepted seed {seed}")

_new_bad = introduced_violations(qa.checks, len(sel.retained), written.elements)
reports["qa"] = {"n_gate_total": qa.n_gate_total,
                 "n_gate_failed": qa.n_gate_failed,
                 "n_introduced": len(_new_bad),
                 "introduced": _new_bad[:20],
                 "failed": [{"check": c.check_id, "requirement": c.requirement,
                             "observed": c.observed} for c in qa.checks
                            if c.status == "fail"]}
(OUT / f"{out14.stem}_qa.json").write_text(json.dumps(qa.to_dict(), indent=1,
                                                      default=float))

# Achieved, not predicted: the dt the finished mesh allows, by the minimum
# altitude of a triangle over sqrt(g*H), which is the measure notebook 392
# and coast_fit both use.
_u3 = written.nodes[written.elements[:, 1]] - written.nodes[written.elements[:, 0]]
_v3 = written.nodes[written.elements[:, 2]] - written.nodes[written.elements[:, 0]]
_area = 0.5 * np.abs(_u3[:, 0] * _v3[:, 1] - _u3[:, 1] * _v3[:, 0])
_side = np.stack([
    np.linalg.norm(written.nodes[written.elements[:, (i + 1) % 3]]
                   - written.nodes[written.elements[:, i]], axis=1)
    for i in range(3)], axis=1)
_dt = (2 * _area / _side.max(axis=1)) / np.sqrt(
    9.81 * written.depths[written.elements].max(axis=1))
# Achieved, per region, on the finished mesh. A patch that passes every gate
# and did not deliver the resolution that was asked for is not a success, and
# nothing else here would notice (third review, finding 16).
_e = np.unique(np.sort(np.vstack([written.elements[:, [0, 1]],
                                  written.elements[:, [1, 2]],
                                  written.elements[:, [2, 0]]]), axis=1), axis=0)
_mid = 0.5 * (written.nodes[_e[:, 0], :2] + written.nodes[_e[:, 1], :2])
_len = np.linalg.norm(written.nodes[_e[:, 0], :2] - written.nodes[_e[:, 1], :2],
                      axis=1)
_pts = shapely.points(_mid[:, 0], _mid[:, 1])
per_region = {}
for _g, _r in regions_m:
    _in = np.asarray(shapely.contains(_g, _pts))
    per_region[_r.name] = {
        "target_h_m": _r.target_h_m,
        "n_edges": int(_in.sum()),
        "median_m": float(np.median(_len[_in])) if _in.any() else None,
        "p90_m": float(np.percentile(_len[_in], 90)) if _in.any() else None,
        "max_m": float(_len[_in].max()) if _in.any() else None,
    }
    if _in.any():
        say(f"achieved in {_r.name}: {_in.sum():,} edges, median "
            f"{np.median(_len[_in]):.1f} m against a {_r.target_h_m:g} m target "
            f"(p90 {np.percentile(_len[_in], 90):.1f}, max {_len[_in].max():.1f})")
reports["achieved"] = {
    "dt_min_s": float(_dt.min()),
    "dt_min_element": int(_dt.argmin()),
    "n_nodes": int(written.n_nodes),
    "n_elements": int(written.n_elements),
    "per_region": per_region,
}
# The operational product is an FVCOM case, not a fort.14. Depth control at
# the open boundary is NOT applied: it rewrites OBC depths, and those nodes
# are frozen. The cor column is the node's latitude, which is what the base's
# own TokyoBay_cor.dat holds (checked: worst difference 7e-10 deg).
_lon, _lat = Transformer.from_crs(f"EPSG:{MESH_EPSG}", "EPSG:4326",
                                  always_xy=True).transform(
    written.nodes[:, 0], written.nodes[:, 1])
_case = export_fvcom_case(written, OUT / "fvcom", recipe.stem,
                          cor=_lat, obc_depth_control=False,
                          obc_type=getattr(base, "obc_type", 1))
reports["fvcom_case"] = {k: str(v) for k, v in _case.items()}
_check = read_fvcom_case(_case["grd"], _case["dep"], _case["obc"])
if not np.array_equal(_check.depths[node_map[node_map >= 0]],
                      base.depths[node_map >= 0]):
    raise SystemExit("the written FVCOM case does not carry the base depths")
say(f"wrote the FVCOM case: {', '.join(sorted(_case))} in {OUT / 'fvcom'}")

say(f"QA {qa.n_gate_total - qa.n_gate_failed}/{qa.n_gate_total} "
    f"({reports['qa']['n_introduced']} introduced by the patch), achieved "
    f"dt {_dt.min():.2f} s (predicted {reports['preflight'][0]['dt_s']:.2f} s)")
(OUT / "report.json").write_text(json.dumps(reports, indent=1, default=float))
if reports["qa"]["n_introduced"]:
    raise SystemExit(
        f"the patch introduces {reports['qa']['n_introduced']} QA violation(s) "
        "the base did not have: "
        + "; ".join(f"{v['check']} at {v['kind']} "
                    f"{v.get('id', v.get('elements'))}"
                    for v in reports["qa"]["introduced"][:6]))
say("done")
