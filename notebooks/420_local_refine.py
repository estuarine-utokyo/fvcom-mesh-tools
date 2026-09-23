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

from fvcom_mesh_tools.bathy_patch import (  # noqa: E402
    edge_slopes,
    interface_lines,
    patch_depths,
)
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.patch import (
    _subdivide,  # noqa: E402
    ambient_size_field,
    base_size_field,
    boundary_after_patch,
    effective_gradation,
    field_gradation,
    filter_shoreline,
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
from fvcom_mesh_tools.walls import (  # noqa: E402
    extract_walls,
    node_walls,
    split_along_walls,
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
elif cfg["coastline"] in ("resample", "resolve"):
    raise SystemExit(f"coastline: {cfg['coastline']} needs a source shoreline; "
                     f"FMESH_LAND={LAND} does not exist. On the hires branch "
                     "that is the OSM land polygons (DATA_INVENTORY.md: "
                     "coastline precedence #1).")

# ---------------------------------------------------------------- geometry
# The recipe declares the region in lon/lat; everything below is in the mesh
# CRS (metres), which is also what DistMesh is handed.  Working in metres is
# what disables generate_mesh's internal tmerc sandwich -- it only engages for
# a bbox that looks like degrees -- so fd/fh are evaluated in the frame they
# were written in.
from pyproj import Transformer  # noqa: E402

to_m = Transformer.from_crs("EPSG:4326", f"EPSG:{MESH_EPSG}", always_xy=True)
to_ll = Transformer.from_crs(f"EPSG:{MESH_EPSG}", "EPSG:4326", always_xy=True)


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


def base_depth_of(lon, lat):
    x, y = to_m.transform(np.asarray(lon), np.asarray(lat))
    # Outside the base mesh there is no depth to inherit, so NaN is the
    # honest answer; preflight refuses it, which is right for a region that
    # is not over this mesh at all.
    return np.ma.filled(_di(x, y), np.nan)


# ---------------------------------------------------------------- the branch
# ONE `if`, taken here (owner, 2026-09-23).  `hires` absent and nothing below
# is reached; everything the recipe did before it existed, it still does.
HIRES = cfg["hires"]
if HIRES is not None:
    say(f"hires: coastline {HIRES['coastline']}, bathymetry "
        f"{HIRES['bathymetry']}, scope {HIRES['scope']}, blend {HIRES['blend']}")
    if HIRES["bathymetry"] == "tokyo_bay":
        say("hires: depths come from the ladder, NOT from the base mesh. "
            "Nothing is floored, capped or smoothed here -- the minimum depth "
            "and the r-factor smoothing are the next step.")

_LADDER = HIRES is not None and HIRES["bathymetry"] == "tokyo_bay"


def depth_of(lon, lat):
    """The field pre-flight judges: the one the mesh will actually carry.

    On the hires branch that is the ladder, not the base mesh -- advertising a
    time step computed on a seabed the run will not have is the mistake this
    avoids.  Depths at or above the datum come back as they are: a tidal flat
    is a tidal flat, and `allow_dry` lets pre-flight report it instead of
    refusing it.
    """
    if not _LADDER:
        return base_depth_of(lon, lat)
    from fvcom_mesh_tools.dem import tokyo_bay as _tb

    d, _rung, _dist = _tb.sample(np.asarray(lon), np.asarray(lat))
    return d


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
                   land=land_ll, allow_dry=_LADDER)
    reports["preflight"].append(pf)
    say(f"preflight {region.name}: transition {pf['transition_m']:.0f} m, "
        f"dt {pf['dt_s']:.2f} s, +{pf['elements_added']:.0f} elements")
    if pf["dt_alert"]:
        say("ALERT " + pf["dt_alert"])
    if pf["land_alert"]:
        say("ALERT " + pf["land_alert"])
    if pf.get("core_at_or_above_datum"):
        say(f"    {pf['core_at_or_above_datum']} of {pf['core_samples']} core "
            f"samples ({100 * pf['core_dry_fraction']:.0f} %) are at or above "
            "the datum -- a tidal flat, which needs WET_DRY_ON and a MIN_DEPTH")

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
conflicts = region_conflicts(sized, [r.name for _, r in regions_m],
                             base_size=base_size_field(base.nodes, base.elements))
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

# The interface the blend's weight measures from: the rim of the RETAINED
# region, minus the base mesh's own boundary.  It is computed on the base
# because retained nodes keep their coordinates, so these lines are the same
# lines after stitching -- and it is the cut the selection ACTUALLY made,
# which is not the analytic buffer: select_patch takes whole faces by centroid
# and then grows the selection to repair pinches.
iface_lines = interface_lines(base.nodes, sel.rim_edges, sel.physical_rim) \
    if HIRES is not None else None
if HIRES is not None:
    say(f"interface: {len(iface_lines.geoms)} edge(s), "
        f"{iface_lines.length / 1000:.2f} km")

# ------------------------------------------------------------------- rim
# The source shoreline for `resample` is picked ONCE per stretch, by the land
# ring the stretch already lies on.  Picking it inside the resampler by
# nearest distance would let a stretch jump to the opposite bank at a strait.
free = np.setdiff1d(np.unique(sel.rim_edges), sel.frozen_nodes)
shore = None
if cfg["coastline"] in ("resample", "resolve") and free.size:
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
        raise SystemExit(f"coastline: {cfg['coastline']} found no source "
                         "shoreline within 3 km of the free rim")
    say(f"source shoreline: {len(shore)} ring(s) within 3 km, nearest "
        f"{min(shapely.distance(pts, ln) for ln in shore):.1f} m from the free rim")

_shl = None
if HIRES is not None and cfg["coastline"] == "resolve" and shore:
    # THE JUDGEMENT the declared grid size implies, made once and applied to
    # the whole hole -- the transition is treated like the region (owner,
    # 2026-09-23).  Everything narrower than the target goes, because a mesh
    # of that size cannot carry it; what survives is handed to oceanmesh,
    # which resamples, culls by area and smooths it, and whose signed
    # distance function is then part of the domain the fill sees.
    _h0 = min(r.target_h_m for _, r in regions_m)
    # Only the land the PATCH touches.  The first version filtered every
    # polygon within 3 km of the free rim, which fragmented industrial coast
    # kilometres away -- 27 rings into 189 -- for no benefit to a hole that
    # never reaches it.  Whole polygons, never clipped: clipping a land
    # polygon invents coastline along the cut.
    _foot = shapely.MultiPoint(
        base.nodes[np.unique(base.elements[sel.removed]), :2]).convex_hull
    _foot = _foot.buffer(max(500.0, 3.0 * float(max(ambient.values()))))
    _keep = [g for g in getattr(land_m, "geoms", [land_m])
             if shapely.intersects(_foot, g)]
    # Two elements across for an AREA, and the rest is not deleted but
    # becomes WALLS (docs/linear_structures_design.md, owner 2026-09-23): a
    # 12 m pier cannot be meshed as land at 30 m, but it can as a line.
    _filtered, _frep = filter_shoreline(_keep, _h0, elements_per_feature=2)
    reports["shoreline_filter"] = _frep
    _walls_src, _wrep = extract_walls(_keep, _filtered, _h0)
    _walls_src = node_walls(_walls_src, snap_m=_h0)
    reports["walls_extracted"] = {k: v for k, v in _wrep.items() if k != "dropped"}
    say(f"walls: {_wrep['n_walls']} extracted, {_wrep['wall_length_m'] / 1000:.2f} km, "
        f"{_wrep['n_dropped']} piece(s) dropped as shorter than {_h0:g} m; "
        f"{_wrep['footprint_given_to_water_m2'] / 1e6:.4f} km2 of structure "
        "footprint becomes water")
    say(f"shoreline filter at h0 = {_h0:g} m: "
        f"{_frep['rings_before']} ring(s) -> {_frep['rings_after']}, "
        f"land lost {_frep['land_lost_m2'] / 1e6:.4f} km2, water lost "
        f"{_frep['water_lost_m2'] / 1e6:.4f} km2, perimeter "
        f"{_frep['perimeter_before_m'] / 1000:.2f} -> "
        f"{_frep['perimeter_after_m'] / 1000:.2f} km")
    _shp = OUT / "shoreline_filtered.shp"
    gpd.GeoDataFrame(geometry=[_filtered], crs=MESH_EPSG).explode(
        index_parts=False).to_file(_shp)
    _b = _filtered.bounds
    _pad = 3.0 * _h0
    _shl = om.Shoreline(str(_shp),
                        (_b[0] - _pad, _b[2] + _pad, _b[1] - _pad, _b[3] + _pad),
                        _h0, crs=f"EPSG:{MESH_EPSG}")
    say(f"oceanmesh Shoreline(h0={_h0:g}): mainland {len(_shl.mainland):,} pt, "
        f"inner {len(_shl.inner):,} pt")
    # The rim is cut from what survived, not from the raw OSM.
    shore = []
    for g in getattr(_filtered, "geoms", [_filtered]):
        for r in [g.exterior, *g.interiors]:
            if shapely.intersects(reach, r):
                shore.append(shapely.LineString(np.asarray(r.coords)))
    if not shore:
        raise SystemExit(f"the h0 = {_h0:g} m filter left no shoreline near "
                         "the free rim; the declared size cannot carry this "
                         "coastline at all")

target = min(r.target_h_m for _, r in regions_m)
# The coastline inside the hole is cut at the LOCAL size, not at the target.
# h_achieved is the same field DistMesh gets, without the 1.2 field-to-bar
# factor, because these are the bar lengths themselves.
# `outside` matters only on the hires branch, and there it matters a lot:
# `resolve` moves the boundary onto the source shoreline, so part of the hole
# is outside the base mesh, and the historical "field maximum out there" is a
# cliff -- measured at a slope of 35-42 against a C4 reference of 0.414,
# where the same patch on the default branch measures 0.375.
_OUTSIDE = "nearest" if HIRES is not None else "max"
h_achieved = patch_sizing(base.nodes, base.elements, sized, distmesh_scale=1.0,
                          outside=_OUTSIDE)
# The base mesh's own boundary, so a replaced stretch can be checked against
# the coastline it does NOT own.  `resolve` moves the boundary onto the source
# and a moved boundary can cross the frozen one -- verify_patch does not look
# for that and matplotlib's TriFinder does, which is how a run with a clean
# size field came back "Triangulation is invalid" for a single crossing pair.
_ub, _cb = np.unique(np.sort(np.vstack(
    [base.elements[:, [0, 1]], base.elements[:, [1, 2]],
     base.elements[:, [2, 0]]]), axis=1), axis=0, return_counts=True)
rc = rim_constraints(base.nodes, sel, size=h_achieved,
                     coastline=cfg["coastline"], shoreline=shore,
                     tolerance_m=cfg["coastline_tolerance_m"],
                     boundary_edges=_ub[_cb == 1])
reports["rim"] = {k: v for k, v in rc.items()
                  if isinstance(v, (int, float, str, bool))}
say("rim: " + json.dumps(reports["rim"]))
if rc.get("n_stretches_kept_to_avoid_a_crossing"):
    _why = {"self": "the source polyline touches itself there -- a feature "
                    "narrower than the local element size",
            "other": "it would cross another replaced stretch",
            "frozen": "it would cross the frozen coastline"}
    say(f"    {rc['n_stretches_kept_to_avoid_a_crossing']} stretch(es) were "
        "KEPT on the base polyline: "
        + "; ".join(f"{n} because {_why.get(k, k)}"
                    for k, n in rc.get("kept_because", {}).items()))
if rc.get("n_coastline_nodes_new", 0) < rc.get("n_coastline_nodes_replaced", 0):
    # The coastline is cut at the LOCAL size, and out at the edge of a
    # transition that is the AMBIENT size.  Measured on the first hires run:
    # a 300 m fishery 2 km offshore put its coastline where the field asks
    # for 400-1700 m elements, and `resolve` replaced 17 base nodes with 13 --
    # a coarser coastline than the base's, which is the failure `preserve`
    # was written to avoid.  Resolving a coastline needs the REGION to reach
    # it, not merely the transition.
    say(f"WARNING the resolved coastline is COARSER than the base's: "
        f"{rc['n_coastline_nodes_replaced']} base node(s) replaced by "
        f"{rc['n_coastline_nodes_new']}. The coastline is cut at the local "
        "size, and here that is the transition's, not the target. To refine "
        "a coastline the region has to contain it.")

# The curves each replaced stretch was cut from.  Saved because the fidelity
# question -- how far the delivered coastline is from the source -- can only
# be answered against THESE, and measuring against "any OSM ring nearby"
# answers a different and easier question: on this patch the two differ by
# 821.6 m against 128.2 m.
if rc["curves"]:
    np.savez(OUT / "coastline_curves.npz",
             **{f"c{k}": np.asarray(c, dtype=float)
                for k, c in enumerate(rc["curves"]) if len(c) > 1})

hole = hole_polygon(rc["pfix"], rc["egfix"])
say(f"hole {hole.area / 1e6:.3f} km2 ({hole.geom_type}, valid={hole.is_valid})")

# ------------------------------------------------------------------- walls
# Each wall is clipped to the hole -- never across the frozen interface, and
# kept a local element clear of it -- rooted by INSERTING its root into the
# coastline rim (a root that stops short of the coast leaves a gap the tide
# goes round), resampled at the local size keeping every corner, and handed
# to the fill as interior constrained edges.  After stitching the mesh is
# SPLIT along them, which is what makes them walls (notebook 427).
WALL_PTS = np.zeros((0, 2))
WALL_SEGS = np.zeros((0, 2), dtype=np.int64)
if HIRES is not None and _shl is not None and _walls_src:
    _margin = 0.5 * float(max(ambient.values()))
    _room = shapely.difference(hole, iface_lines.buffer(_margin)) \
        if not iface_lines.is_empty else hole
    _coast = shapely.difference(hole.boundary, iface_lines.buffer(1.0)) \
        if not iface_lines.is_empty else hole.boundary
    # And kept clear of the coast, except where they are rooted in it.  A
    # wall running within a few metres of the shore -- the sliver an area
    # filter leaves along a quay it keeps -- put elements of 3.1 and 4.8 deg
    # between itself and the coast, and closed off pockets the open boundary
    # could not reach (QA: 3 components, 7 elements unreachable).  Removing
    # the part within 0.4 of an element of the shore leaves such a sliver as
    # stubs shorter than L_min; a pier keeps its body, and its end, now
    # 0.4 h from the coast, is rooted by the rule below.
    _room = shapely.difference(_room, _coast.buffer(0.4 * target))
    _pieces = []
    for w in _walls_src:
        g = shapely.intersection(w, _room)
        for q in getattr(g, "geoms", [g]):
            if q.geom_type == "LineString" and q.length >= target:
                _pieces.append(np.asarray(q.coords)[:, :2])
    rim_xy = np.asarray(rc["pfix"], dtype=float)
    rim_eg = np.asarray(rc["egfix"], dtype=np.int64)
    rim_base = np.asarray(rc["pfix_base"], dtype=np.int64)
    _index: dict = {}
    _pts: list = []
    _segs: list = []
    n_rooted = 0

    def _key(xy):
        return (round(float(xy[0]), 6), round(float(xy[1]), 6))

    def _root(xy):
        """Put a wall's root ON the coastline rim; return its pfix row."""
        global rim_xy, rim_eg, rim_base
        h_here = float(h_achieved(np.asarray([xy]))[0])
        a, b = rim_xy[rim_eg[:, 0]], rim_xy[rim_eg[:, 1]]
        ab = b - a
        s = np.clip(np.einsum("ij,ij->i", xy - a, ab)
                    / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-12), 0, 1)
        foot = a + s[:, None] * ab
        k = int(np.argmin(np.linalg.norm(foot - xy, axis=1)))
        # Judged at the FOOT, not at the wall's end.  The end may sit half an
        # element off the coast, so a foot landing ON a rim vertex passed a
        # test made at the end and was inserted a second time -- two fixed
        # points at one position, which the mesher collapses.
        d = np.linalg.norm(rim_xy - foot[k], axis=1)
        if d.min() <= 0.3 * h_here:
            return int(d.argmin())
        i, j = rim_eg[k]
        new = len(rim_xy)
        rim_xy = np.vstack([rim_xy, foot[k]])
        rim_base = np.append(rim_base, -1)
        rim_eg = np.vstack([np.delete(rim_eg, k, axis=0), [[i, new], [new, j]]])
        return new

    for c in _pieces:
        ends = []
        for e in (0, -1):
            h_e = float(h_achieved(np.asarray([c[e]]))[0])
            if float(shapely.distance(_coast, shapely.Point(c[e]))) <= 0.5 * h_e:
                ends.append(_root(c[e]))
                n_rooted += 1
            else:
                ends.append(None)
        walk = _subdivide(c, h_achieved)
        ids = []
        for k_, xy in enumerate(walk):
            if k_ == 0 and ends[0] is not None:
                ids.append(("rim", ends[0]))
                continue
            if k_ == len(walk) - 1 and ends[1] is not None:
                ids.append(("rim", ends[1]))
                continue
            # A point closer than a quarter element to one already placed --
            # on the rim or on another wall -- IS that point.  Two fixed
            # points that close are merged by the mesher anyway, and the
            # first run with walls failed every seed on the resulting edge
            # from a node to itself.
            h_here = float(h_achieved(np.asarray([xy]))[0])
            d_rim = np.linalg.norm(rim_xy - xy, axis=1)
            if d_rim.min() <= 0.25 * h_here:
                ids.append(("rim", int(d_rim.argmin())))
                continue
            if _pts:
                d_w = np.linalg.norm(np.asarray(_pts) - xy, axis=1)
                if d_w.min() <= 0.25 * h_here:
                    ids.append(("wall", int(d_w.argmin())))
                    continue
            key = _key(xy)
            if key not in _index:
                _index[key] = len(_pts)
                _pts.append(xy)
            ids.append(("wall", _index[key]))
        ids = [x for k_, x in enumerate(ids) if k_ == 0 or x != ids[k_ - 1]]
        if len(ids) == 2 and ids[0][0] == "wall" and ids[1][0] == "wall":
            # A detached wall of ONE edge has two free tips and nothing
            # between them: neither end has a second sector to copy into,
            # so the edge stays interior and the split refuses it -- one
            # such edge failed every seed.  A slit needs a node inside it.
            mid = 0.5 * (_pts[ids[0][1]] + _pts[ids[1][1]])
            _pts.append(mid)
            ids.insert(1, ("wall", len(_pts) - 1))
        for (ka, ia), (kb, ib) in zip(ids[:-1], ids[1:]):
            if (ka, ia) != (kb, ib):
                _segs.append(((ka, ia), (kb, ib)))
    rc["pfix"], rc["egfix"], rc["pfix_base"] = rim_xy, rim_eg, rim_base
    n_rim = len(rim_xy)
    WALL_PTS = np.asarray(_pts, dtype=float).reshape(-1, 2)
    WALL_SEGS = np.asarray([[ia if ka == "rim" else n_rim + ia,
                             ib if kb == "rim" else n_rim + ib]
                            for (ka, ia), (kb, ib) in _segs],
                           dtype=np.int64).reshape(-1, 2)
    if len(WALL_SEGS):
        WALL_SEGS = np.unique(np.sort(WALL_SEGS, axis=1), axis=0)
    # A wall edge along the coastline rim is already boundary; it has no
    # second side to split off, and the mesher already honours it.
    _rimset = {tuple(sorted(e)) for e in rim_eg.tolist()}
    WALL_SEGS = np.asarray([e for e in WALL_SEGS.tolist()
                            if tuple(e) not in _rimset],
                           dtype=np.int64).reshape(-1, 2)
    # No two lines may meet at under 30 degrees where a wall is involved.
    # The elements between two lines that meet at an angle cannot be wider
    # than the angle, and C1 wants 30: the crossing stub of an L-shaped
    # breakwater left elements of 3.1 and 4.8 deg even after the pocket it
    # closed was opened.  The WALL edge of an acute pair goes, never the
    # coast's; a detached single edge left behind goes too, since a slit
    # needs a node inside it.
    _all_xy = np.vstack([rim_xy, WALL_PTS])
    _rim_set = {tuple(sorted(e)) for e in rim_eg.tolist()}
    n_acute = 0
    while len(WALL_SEGS):
        inc: dict = {}
        for k, (a, b) in enumerate(WALL_SEGS.tolist()):
            inc.setdefault(a, []).append(("w", k, b))
            inc.setdefault(b, []).append(("w", k, a))
        for a, b in rim_eg.tolist():
            if a in inc:
                inc[a].append(("r", -1, b))
            if b in inc:
                inc[b].append(("r", -1, a))
        drop = None
        for v, lst in inc.items():
            for i in range(len(lst)):
                for j in range(i + 1, len(lst)):
                    if lst[i][0] == "r" and lst[j][0] == "r":
                        continue
                    u1 = _all_xy[lst[i][2]] - _all_xy[v]
                    u2 = _all_xy[lst[j][2]] - _all_xy[v]
                    ang = np.degrees(np.arccos(np.clip(
                        u1 @ u2 / (np.linalg.norm(u1) * np.linalg.norm(u2) + 1e-12), -1, 1)))
                    if ang < 30.0:
                        drop = lst[i][1] if lst[i][0] == "w" else lst[j][1]
                        if lst[i][0] == "w" and lst[j][0] == "w":
                            li = np.linalg.norm(u1)
                            lj = np.linalg.norm(u2)
                            drop = lst[i][1] if li <= lj else lst[j][1]
                        break
                if drop is not None:
                    break
            if drop is not None:
                break
        if drop is None:
            break
        WALL_SEGS = np.delete(WALL_SEGS, drop, axis=0)
        n_acute += 1
        # a lone edge whose ends touch nothing else and are not on the rim
        deg = np.bincount(WALL_SEGS.ravel(), minlength=len(_all_xy)) if len(WALL_SEGS) \
            else np.zeros(len(_all_xy), int)
        lone = [k for k, (a, b) in enumerate(WALL_SEGS.tolist())
                if deg[a] == 1 and deg[b] == 1 and a >= n_rim and b >= n_rim]
        if lone:
            WALL_SEGS = np.delete(WALL_SEGS, lone, axis=0)
    # A wall point no edge uses any more would be a lone fixed point in the
    # water -- not a wall, and a small element waiting to happen.
    _used = np.unique(WALL_SEGS[WALL_SEGS >= n_rim]) if len(WALL_SEGS) else \
        np.zeros(0, dtype=np.int64)
    _remap = np.full(len(_all_xy), -1, dtype=np.int64)
    _remap[:n_rim] = np.arange(n_rim)
    _remap[_used] = n_rim + np.arange(len(_used))
    WALL_PTS = _all_xy[_used] if len(_used) else np.zeros((0, 2))
    WALL_SEGS = _remap[WALL_SEGS] if len(WALL_SEGS) else WALL_SEGS
    reports["walls"] = {"n_pieces_in_hole": len(_pieces), "n_rooted_ends": n_rooted,
                        "n_wall_edges_dropped_for_an_acute_angle": n_acute,
                        "n_wall_points": int(len(WALL_PTS)),
                        "n_wall_edges": int(len(WALL_SEGS))}
    say(f"walls in the hole: {len(_pieces)} piece(s), {n_rooted} end(s) rooted on "
        f"the coast, {len(WALL_PTS)} constrained point(s), {len(WALL_SEGS)} edge(s); "
        f"{n_acute} dropped for meeting another line at under 30 deg")
PFIX_ALL = np.vstack([np.asarray(rc["pfix"], dtype=float), WALL_PTS])
# What the fill was given, kept so a wall's geometry can be inspected
# without re-running the whole cut.
np.savez(OUT / "fill_constraints.npz", pfix=PFIX_ALL,
         egfix=np.vstack([np.asarray(rc["egfix"], dtype=np.int64), WALL_SEGS]),
         n_rim=len(rc["pfix"]))
EGFIX_ALL = np.vstack([np.asarray(rc["egfix"], dtype=np.int64), WALL_SEGS])
PFIX_BASE_ALL = np.concatenate([np.asarray(rc["pfix_base"], dtype=np.int64),
                                np.full(len(WALL_PTS), -1, dtype=np.int64)])

# The slope the field actually has, measured on the hole it will be meshed
# in -- not the per-region formula, which omits the ambient term and says
# nothing about where two regions meet.
fslope = field_gradation(h_achieved, hole)
reports["field_gradation"] = fslope
if not fslope.get("measured", True):
    raise SystemExit("the sizing field could not be measured on this hole "
                     f"({fslope['n_samples']} samples at "
                     f"{fslope['spacing_m']:.1f} m): an unmeasured slope is "
                     "not a gentle one")
say(f"field slope: max {fslope['max_slope']:.3f}, p99 {fslope['p99_slope']:.3f} "
    f"against a {fslope['c4_reference_gradation']:.3f} reference; "
    f"{100 * fslope['fraction_above_reference']:.2f} % of samples above it")
if fslope["max_slope"] > fslope["c4_reference_gradation"]:
    say("  ALERT the field is locally steeper than two similar triangles one "
        "size apart can be and still pass C4. Expect the fill to be hard: "
        "widen a transition, move a region, or bring the targets closer")

# ------------------------------------------------------------------ fill
fh = patch_sizing(base.nodes, base.elements, sized,
                  distmesh_scale=DISTMESH_SCALE, outside=_OUTSIDE)
shapely.prepare(hole)
boundary = shapely.boundary(hole)


def fd(points):
    """Signed distance to the hole: negative inside, in metres.

    The shoreline enters through the RIM, not through here.  Intersecting
    this with oceanmesh's signed distance function for the filtered coastline
    was tried and does not work: the filter moves the coastline -- on the
    Kimitsu patch it closed 4.46 km2 of water -- so base rim points that the
    FROZEN mesh says are water fall on the land side of the new shoreline,
    and DistMesh prunes them.  Measured, 19 of 248 fixed points were dropped,
    one of them 1,434 m from where it was asked for.

    The frozen rim is not negotiable and the shoreline is not either, so the
    reconciliation has to happen where they meet, and it does: every coastal
    rim point is cut from the FILTERED shoreline, and this polygon is built
    from those points.  The domain inside the hole is therefore the source's,
    and the base supplies only the interface -- which is the one thing it
    must supply, because that is where the patch meets what it promised not
    to touch.
    """
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

# What "delivered" means, and it is a statement about WATER, not about edges.
# Edge statistics are counted per edge, so a request whose northern half is
# refined and whose southern half was never touched still shows a median at
# the target: the fine half owns 1,100 edges and the coarse half owns eight.
# On the delivered circle mesh a request with 46 % of its area left at the
# base size passed the edge-median gate (fifth review).
#
# So the region's area is sampled and each sample is asked how big the
# element covering it is, as an equivalent edge. Measured on the five
# delivered regions, against their own targets:
#
#   region            median ratio   area within 1.25x
#   circle                    0.95               0.970
#   polygon                   0.95               0.978
#   two beds, north           0.98               0.998
#   two beds, east            0.74               0.997
#   two beds, channel         0.98               0.980
#
# and on the two requests the fifth review showed were not delivered:
# the two-lobe region 14.87 / 0.449, the request buffered by 300 m
# 1.69 / 0.330. The thresholds below sit in that gap. COVERAGE_TOLERANCE is
# not 1.05: a DistMesh fill leaves individual cells above the target -- the
# channel has 69 % of its area within 1.05 -- and a mesh that delivers the
# resolution must not be rejected for that.
RESOLUTION_TOLERANCE = 1.05        # on the area-weighted median


def achieved_per_region(mesh):
    """What each declared region actually got, and whether that is its target.

    Reports both measures: the area-weighted resolution, which is what the
    gate is, and the edge statistics inside the region, which is what the
    earlier reports quoted and what a reader compares against them.

    A region with no water in it at all -- every sample outside the mesh --
    is a miss, not a pass. A region that reaches over land is not: the
    fraction outside the mesh is reported so it can be seen.

    The two coverage numbers and the import are local so that this function
    can be lifted out and exercised on its own, which is how both reviews
    tested the gate rather than a restatement of it.
    """
    from fvcom_mesh_tools.patch import region_resolution

    coverage_tolerance = 1.25   # a cell this much over target still counts
    coverage_fraction = 0.95    # of the requested water, at least
    e = np.unique(np.sort(np.vstack([mesh.elements[:, [0, 1]],
                                     mesh.elements[:, [1, 2]],
                                     mesh.elements[:, [2, 0]]]), axis=1), axis=0)
    mid = 0.5 * (mesh.nodes[e[:, 0], :2] + mesh.nodes[e[:, 1], :2])
    length = np.linalg.norm(mesh.nodes[e[:, 0], :2] - mesh.nodes[e[:, 1], :2],
                            axis=1)
    pts = shapely.points(mid[:, 0], mid[:, 1])
    per_region, missed = {}, []
    for geom, region in regions_m:
        inside = np.asarray(shapely.contains(geom, pts))
        stat = region_resolution(mesh.nodes, mesh.elements, geom,
                                 region.target_h_m,
                                 coverage_tolerance=coverage_tolerance)
        stat["n_edges"] = int(inside.sum())
        stat["edge_median_m"] = float(np.median(length[inside])) \
            if inside.any() else None
        stat["edge_p90_m"] = float(np.percentile(length[inside], 90)) \
            if inside.any() else None
        stat["edge_max_m"] = float(length[inside].max()) if inside.any() else None
        if stat["median_ratio"] is None:
            stat["miss"] = "none of the region is inside the mesh"
        elif stat["median_ratio"] > RESOLUTION_TOLERANCE:
            stat["miss"] = (f"the median cell is {stat['median_m']:.1f} m "
                            f"against a {region.target_h_m:g} m target")
        elif stat["covered_fraction"] < coverage_fraction:
            stat["miss"] = (
                f"only {100 * stat['covered_fraction']:.1f} % of the water is "
                f"within {coverage_tolerance:g}x the {region.target_h_m:g} m "
                f"target ({100 * coverage_fraction:.0f} % required)")
        if "miss" in stat:
            missed.append(region.name)
        per_region[region.name] = stat
    return per_region, missed


def walls_cut_off(elements, n_nodes, copy_of, wall_edges, obc_nodes, xy=None):
    """Indices of the wall edges whose PIECE borders water cut off from the OBC.

    A piece is a connected run of wall edges, by original node id.  Water is
    cut off when a component of the split mesh holds no open-boundary node.
    """
    tri = np.asarray(elements, dtype=np.int64)
    parent = list(range(n_nodes))

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for a, b, c in tri.tolist():
        for x, y in ((a, b), (b, c)):
            ra, rb = find(x), find(y)
            if ra != rb:
                parent[ra] = rb
    good = {find(n) for n in obc_nodes if n < n_nodes}
    used = np.unique(tri)
    stranded = {int(copy_of[n]) for n in used.tolist() if find(n) not in good}
    if not stranded:
        return []
    we = np.asarray(wall_edges, dtype=np.int64)
    wp = list(range(int(we.max()) + 1)) if len(we) else []

    def wfind(k):
        while wp[k] != k:
            wp[k] = wp[wp[k]]
            k = wp[k]
        return k

    for a, b in we.tolist():
        ra, rb = wfind(a), wfind(b)
        if ra != rb:
            wp[ra] = rb
    del wfind
    # The SHORTEST wall edge that touches the stranded water, one at a time.
    # Withdrawing the whole connected piece removed an entire arm of the
    # L-shaped breakwater to open a pocket its short crossing stub had made.
    touching = [k for k, (a, b) in enumerate(we.tolist())
                if a in stranded or b in stranded]
    if not touching:
        return []
    if xy is None:
        return touching
    lengths = [float(np.linalg.norm(xy[we[k, 0]] - xy[we[k, 1]])) for k in touching]
    return [touching[int(np.argmin(lengths))]]


def wall_pairs(out):
    """Every pair of nodes a wall split made coincident ON PURPOSE."""
    co = out.get("copy_of") if isinstance(out, dict) else None
    if co is None:
        return None
    co = np.asarray(co)
    groups: dict = {}
    for k, s in enumerate(co.tolist()):
        groups.setdefault(s, []).append(k)
    return [(a, b) for g in groups.values() if len(g) > 1
            for i, a in enumerate(g) for b in g[i + 1:]]


def patch_violations(qa_checks, written):
    """The QA failures this patch is answerable for, and the one exception.

    ``introduced_violations`` decides what the patch caused rather than what
    it inherited.  On the hires branch there is a second question: run_qa's
    floor is 2 m and this branch writes the depths as the source gives them,
    so a tidal flat fails it.  That is a tidal flat under wetting and drying,
    not a defect the seed can be blamed for, and gating it would refuse the
    very field the option exists to deliver.

    It is a function, and both the seed loop and the final gate call it,
    because the first version of this branch downgraded the check in the loop
    and not at the end -- so every seed passed and the finished mesh was
    rejected by a rule the search had not been applying.

    Returns ``(blamed, n_reported)``; the second number is the point, because
    the absence of the gate must not become the absence of the report.
    """
    bad = introduced_violations(qa_checks, len(sel.retained), written.elements)
    if not _LADDER:
        return bad, 0
    shallow = [b for b in bad if b.get("check") == "min_depth_clip"]
    return [b for b in bad if b.get("check") != "min_depth_clip"], len(shallow)


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
        pfix=PFIX_ALL, egfix=EGFIX_ALL, cleanup="none")
    say(f"filled: NP={len(p):,} NE={len(t):,}")

    p, t = collapse_thin_triangles(p, t, min_qual=0.25, pfix=PFIX_ALL)
    p, t = direct_smoother_lur(p, t, pfix=PFIX_ALL)
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
        PFIX_ALL, PFIX_BASE_ALL)
    out["stitch"] = {k: v for k, v in st.items() if k != "pfix_new"}
    say("stitch: " + json.dumps(out["stitch"]))
    # Before the repair, while a fixed point is still exactly where it was
    # put: the repair slides boundary nodes along their curve, so matching
    # pfix by coordinate afterwards fails.  Node ids do not change in the
    # repair, so the edge set built here stays valid.
    want_boundary = boundary_after_patch(base.elements, sel, rc, node_map,
                                         st["pfix_new"][:len(rc["pfix"])])
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

    # -------------------------------------------------------------- the split
    copy_of = np.arange(len(nodes))
    wall_node = np.zeros(len(nodes), dtype=bool)
    if len(WALL_SEGS):
        _we = np.asarray(st["pfix_new"], dtype=np.int64)[WALL_SEGS]
        _self = _we[:, 0] == _we[:, 1]
        if _self.any():
            # two fixed points the mesher merged: the edge between them is
            # gone, not a wall, and it is counted rather than raised
            out["wall_edges_merged_by_mesher"] = int(_self.sum())
            _we = _we[~_self]
        # A wall that closes water off -- two walls crossing into a small
        # triangle, or a wall rooted twice on the coast around a pocket --
        # leaves a piece the open boundary cannot reach.  The first harbour
        # with walls had two, one of them two elements of 3.1 and 4.8 deg.
        # Every wall PIECE (connected run of wall edges) touching such a
        # piece is withdrawn and the split redone; the count is reported.
        _obc0 = set(np.concatenate([np.asarray(s) for s in base.open_boundaries]).tolist())
        _obc_new = {int(node_map[n]) for n in _obc0 if node_map[n] >= 0}
        _pre = (nodes, elements)
        withdrawn = 0
        for _round in range(40):
            nodes, elements, copy_of, srep = split_along_walls(_pre[0], _pre[1], _we)
            pieces = walls_cut_off(elements, len(nodes), copy_of, _we, _obc_new,
                                   xy=_pre[0])
            if not pieces:
                break
            keep_e = ~np.isin(np.arange(len(_we)), pieces)
            withdrawn += int((~keep_e).sum())
            _we = _we[keep_e]
            if not len(_we):
                nodes, elements = _pre
                copy_of = np.arange(len(nodes))
                srep = {"n_wall_edges": 0, "n_copies": 0, "n_free_tips": 0,
                        "n_components": 1, "pairs": []}
                break
        out["wall_edges_withdrawn_for_closing_water_off"] = withdrawn
        if withdrawn:
            say(f"walls: {withdrawn} wall edge(s) withdrawn -- they closed water "
                "off from the open boundary")
        depths = depths[copy_of]
        _wn = set(np.unique(_we).tolist())
        wall_node = np.isin(copy_of, list(_wn))
        _ws = {tuple(sorted(e)) for e in _we.tolist()}
        _ub2, _cb2 = np.unique(np.sort(np.vstack(
            [elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]]]), axis=1),
            axis=0, return_counts=True)
        # The expectation was written in PRE-split ids.  A coast node where a
        # wall is rooted is duplicated too, so the coastline edge on one side
        # of the pier now runs to the copy: every seed failed verify with 12
        # boundary edges "unexpected" and the same 12 "missing".  Re-express
        # each expected edge through copy_of; one that has no post-split
        # boundary edge at all stays in its old form, so it is still
        # reported missing rather than quietly dropped.
        _by_orig: dict = {}
        for e in _ub2[_cb2 == 1].tolist():
            _by_orig.setdefault(tuple(sorted(copy_of[e].tolist())), []).append(
                tuple(sorted(e)))
        _want = set()
        for e in set(want_boundary) | _ws:
            found = _by_orig.get(tuple(sorted(e)))
            _want.update(found if found else [tuple(sorted(e))])
        want_boundary = _want
        out["want_boundary"] = want_boundary
        out["walls"] = {k: v for k, v in srep.items() if k != "pairs"}
        say(f"walls split: {srep['n_wall_edges']} edge(s), {srep['n_copies']} "
            f"node(s) duplicated, {srep['n_free_tips']} free tip(s), "
            f"{srep['n_components']} component(s)")
    out["copy_of"] = copy_of

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
    # A wall node does not move: it is on a structure, not on a curve it
    # could slide along, and the nearest coastline curve is not its own.
    movable = is_new & ~on_boundary & ~wall_node
    slidable = is_new & on_boundary & ~wall_node
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

    if _LADDER:
        # The ladder, sampled at the FINAL coordinates -- after improve_patch,
        # because the repair slides boundary nodes along their coastline curve
        # and a depth sampled before the move belongs to a coordinate the mesh
        # no longer has.  Only new nodes: every retained node keeps its base
        # depth bit for bit, which is what makes the frozen contract checkable.
        _idx = np.flatnonzero(is_new)
        _lo, _la = to_ll.transform(nodes[_idx, 0], nodes[_idx, 1])
        _d, binfo = patch_depths(
            np.column_stack([_lo, _la]), nodes[_idx], depths[_idx],
            regions=[g for g, _ in regions_m], interface=iface_lines,
            scope=HIRES["scope"], blend=HIRES["blend"])
        depths = depths.copy()
        depths[_idx] = _d
        binfo["slopes"] = edge_slopes(nodes, elements, depths, is_new)
        out["bathymetry"] = binfo
        say(f"hires depths: {binfo['n_from_m7001']} m7001 / "
            f"{binfo['n_from_grid30']} grid30 / {binfo['n_from_kanto']} kanto / "
            f"{binfo['n_extrapolated']} extrapolated, "
            f"{binfo['depth_min_m']:.2f}..{binfo['depth_max_m']:.2f} m, "
            f"{binfo['n_at_or_below_zero']} at or above the datum, "
            f"worst |source-base| {binfo['max_source_minus_base_m']:.2f} m")
        if binfo["n_extrapolated"]:
            say(f"    {binfo['n_extrapolated']} node(s) were covered by no "
                f"product and were extrapolated up to "
                f"{binfo['extrapolated_distance_max_m']:.0f} m")
        _sl = binfo["slopes"]["seam"]
        say(f"    seam: {_sl['n']} retained-to-new edges, slope max "
            f"{_sl.get('slope_max', 0):.4f} m/m"
            + (f", r max {_sl['r_max']:.4f}" if _sl.get("r_max") is not None
               else " (r undefined on an intertidal pair)"))

    # Inherit the base's r-factor property, not just its values. The base is
    # m7001tp_rfac0p2_cap300: every one of its edges satisfies r <= 0.2, and
    # interpolation does not carry that across a new edge joining different
    # base elements -- ten new edges came out above it, the worst at 0.3075.
    # Only new nodes' depths move; every base depth is untouched, and the
    # report says how far a new one was pulled.
    if cfg["rfactor_limit"] != "off" and not _LADDER:
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
    # BOTH ends new.  A chord with one retained end runs out past the end of
    # the substring the stretch was cut from, and the distance from a curve
    # that has stopped is not a fidelity measure: on the first hires run that
    # put one junction chord 821.6 m "from the source" while every chord
    # inside the resolved stretch was within 78.5 m.  The junction chords are
    # measured separately, against the whole source, below.
    # Walls are boundary now, and they are not coastline: measuring them
    # against the coastline curves reported a 775.6 m "departure" that was a
    # breakwater's distance from the shore.
    _bnd = _bnd[~wall_node[_bnd].all(axis=1)]
    _new_bnd = _bnd[is_new[_bnd].all(axis=1)]
    _junction = _bnd[is_new[_bnd].any(axis=1) & ~is_new[_bnd].all(axis=1)]
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
    # Against the SOURCE, which on the resolve branch the curve is not: there
    # the curve is the delivered polyline, so a departure from it is zero by
    # construction and says nothing.  This is the fidelity number.
    if len(_new_bnd) and shore:
        _srcs = shapely.MultiLineString([np.asarray(ln.coords) for ln in shore])
        _nn = np.unique(_new_bnd)
        _nn = _nn[is_new[_nn]]
        imp["coastline_source_departure_m"] = float(shapely.distance(
            shapely.points(nodes[_nn]), _srcs).max()) if len(_nn) else 0.0
        imp["coastline_source_departure_median_m"] = float(np.median(
            shapely.distance(shapely.points(nodes[_nn]), _srcs))) \
            if len(_nn) else 0.0
    else:
        imp["coastline_source_departure_m"] = 0.0
        imp["coastline_source_departure_median_m"] = 0.0
    imp["n_coastline_chords"] = int(len(_new_bnd))
    imp["n_junction_chords"] = int(len(_junction))
    # The junctions, against the WHOLE source rather than one stretch's
    # substring, because that is the only reference that reaches them.
    if len(_junction) and shore:
        _ja, _jb = nodes[_junction[:, 0]], nodes[_junction[:, 1]]
        _jf = np.linspace(0.0, 1.0, 9)[:, None, None]
        _js = (_ja[None] + _jf * (_jb - _ja)[None]).reshape(-1, 2)
        imp["junction_departure_m"] = float(shapely.distance(
            shapely.points(_js),
            shapely.MultiLineString([np.asarray(ln.coords) for ln in shore])).max())
    else:
        imp["junction_departure_m"] = 0.0
    out["improve"] = imp
    if imp["coastline_departure_m"] > cfg["coastline_tolerance_m"]:
        say(f"  seed {seed}: the repair moved the coastline "
            f"{imp['coastline_departure_m']:.0f} m from the curve it was cut "
            f"from, past the {cfg['coastline_tolerance_m']:g} m tolerance")
        return None, out
    say(f"seam repair: {imp['n_flips']} flips, {imp['n_moves']} moves, "
        f"coastline departure {imp['coastline_departure_m']:.1f} m over "
        f"{imp['n_coastline_chords']} chord(s) (from OSM: median "
        f"{imp['coastline_source_departure_median_m']:.1f} m, max "
        f"{imp['coastline_source_departure_m']:.1f} m), junctions "
        f"{imp['junction_departure_m']:.1f} m over {imp['n_junction_chords']} "
        f"(nodes {imp['coastline_node_departure_m']:.2f} m), "
        f"angles {imp['min_angle_deg']:.2f}-{imp['max_angle_deg']:.2f} deg "
        f"({int(movable.sum()):,} movable, {int(slidable.sum()):,} slidable nodes, "
        f"{int(mutable_faces.sum()):,} mutable faces)")

    # ---------------------------------------------------------------- verify
    ver = verify_patch(base.nodes, base.depths, base.elements, sel,
                       nodes, elements, depths, node_map,
                       open_boundaries=base.open_boundaries,
                       expected_boundary=want_boundary, copy_of=copy_of)
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
        expected_boundary=out["want_boundary"], copy_of=out.get("copy_of"))
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
                                    if k not in ("want_boundary", "copy_of")})
        save_report()
        continue
    written, mesh = serialise(candidate, out, out14)
    if written is None:
        reports["attempts"].append({k: v for k, v in out.items()
                                    if k not in ("want_boundary", "copy_of")})
        save_report()
        continue
    qa = run_qa(written, name=out14.stem, path=out14, max_offenders=10_000,
                allowed_duplicate_pairs=wall_pairs(out))
    # What the patch is answerable for. A refinement may not be held to a
    # standard its base does not meet: the goto2023 production mesh fails C1
    # at one element 18 km from Futtsu, the contract freezes that element,
    # and an absolute gate blamed every seed for it.
    new_bad, _n_shallow = patch_violations(qa.checks, written)
    out["min_depth_reported_not_gated"] = _n_shallow
    if _n_shallow:
        say(f"    seed {seed}: min-depth is REPORTED not gated on this branch "
            f"-- {_n_shallow} offender(s), minimum {written.depths.min():.2f} m")
    per_region, missed = achieved_per_region(written)
    out["achieved_per_region"] = per_region
    for _name in missed:
        say(f"    seed {seed}: {_name} did not get what it asked for -- "
            f"{per_region[_name]['miss']}")
    out["qa"] = {"n_gate_total": qa.n_gate_total,
                 "n_gate_failed": qa.n_gate_failed,
                 "n_introduced": len(new_bad),
                 "introduced": new_bad[:20],
                 "failed": [{"check": c.check_id, "requirement": c.requirement,
                             "observed": c.observed} for c in qa.checks
                            if c.status == "fail"]}
    reports["attempts"].append({k: v for k, v in out.items()
                                if k not in ("want_boundary", "copy_of")})
    save_report()
    say(f"    seed {seed}: QA {qa.n_gate_total - qa.n_gate_failed}/"
        f"{qa.n_gate_total}, {len(new_bad)} introduced by the patch"
        + ("" if not new_bad else "  " + "; ".join(
            f"{v['check']} at {v['kind']} {v.get('id', v.get('elements'))}"
            for v in new_bad[:4])))
    if best is None or (len(missed), len(new_bad)) < best[:2]:
        best = (len(missed), len(new_bad), seed, candidate, out, written, qa)
    if not new_bad and not missed:
        break

if best is None:
    save_report()
    raise SystemExit("no seed produced a mesh that keeps the frozen-zone "
                     f"contract; see attempts in {OUT / 'report.json'}")
_n_missed, _, seed, candidate, out, written, qa = best
if _n_missed:
    save_report()
    raise SystemExit(
        "no seed delivered the resolution that was asked for: "
        + "; ".join(f"{k}: {v['miss']}" for k, v
                    in best[4]["achieved_per_region"].items() if "miss" in v)
        + f"; see attempts in {OUT / 'report.json'}")
nodes, elements, depths, node_map = candidate
reports.update({k: v for k, v in out.items()
                if k not in ("seed", "want_boundary", "copy_of")})
reports["seed"] = seed
np.save(OUT / "node_map.npy", node_map)
reports["mesh"] = str(out14)
# Always, not only when the accepted seed is not the first one tried: the
# file on disk is whichever attempt ran last, and when none passed that is
# not the best one.  A failed artefact investigated with another attempt's
# node map and QA is worse than no artefact (second review, finding 4).
written, mesh = serialise(candidate, out, out14)
qa = run_qa(written, name=out14.stem, path=out14, max_offenders=10_000,
            allowed_duplicate_pairs=wall_pairs(out))
say(f"accepted seed {seed}")

_new_bad, _n_shallow = patch_violations(qa.checks, written)
if _n_shallow:
    say(f"min-depth is REPORTED not gated on this branch -- {_n_shallow} "
        f"offender(s), minimum {written.depths.min():.2f} m. The depths are "
        "the source's; the floor and the smoothing are the next step.")
reports["qa"] = {"n_gate_total": qa.n_gate_total,
                 "n_gate_failed": qa.n_gate_failed,
                 "min_depth_reported_not_gated": _n_shallow,
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
# A dry element has no wave speed, and sqrt(g * H) for H <= 0 is NaN -- which
# is what the first coastal hires run reported as its achieved time step.  The
# floor is a DIAGNOSTIC one, not a depth written anywhere: shallower water
# gives a LARGER dt, so an element clamped to it cannot become the binding one,
# and the number this reports is still the real constraint (owner, 2026-09-23:
# "the wave speed can be computed if you give it a sensible minimum depth").
DT_MIN_DEPTH_M = 0.05
_hmax = np.maximum(written.depths[written.elements].max(axis=1), DT_MIN_DEPTH_M)
_n_dry_elements = int((written.depths[written.elements].max(axis=1)
                       <= DT_MIN_DEPTH_M).sum())
_dt = (2 * _area / _side.max(axis=1)) / np.sqrt(9.81 * _hmax)
# Achieved, per region, on the finished mesh -- the same measure the seed
# loop gated on, recomputed on the file that was actually written.
per_region, missed = achieved_per_region(written)
if missed:
    raise SystemExit("the written mesh does not deliver the resolution that "
                     "was asked for: "
                     + "; ".join(f"{k}: {per_region[k]['miss']}" for k in missed))
for _name, _st in per_region.items():
    say(f"achieved in {_name}: median cell {_st['median_m']:.1f} m against a "
        f"{_st['target_h_m']:g} m target, {100 * _st['covered_fraction']:.1f} % "
        f"of the water within {_st['coverage_tolerance']:g}x "
        f"({_st['n_samples']:,} area samples at {_st['spacing_m']:.1f} m; "
        f"edges: {_st['n_edges']:,}, median {_st['edge_median_m']:.1f} m, "
        f"p90 {_st['edge_p90_m']:.1f}, max {_st['edge_max_m']:.1f})")
reports["achieved"] = {
    "dt_min_s": float(_dt.min()),
    "dt_min_element": int(_dt.argmin()),
    "dt_wave_speed_floor_m": DT_MIN_DEPTH_M,
    "n_elements_at_or_above_datum": _n_dry_elements,
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
    f"dt {_dt.min():.2f} s (predicted {reports['preflight'][0]['dt_s']:.2f} s"
    + (f", {_n_dry_elements} element(s) at or above the datum, wave speed "
       f"floored at {DT_MIN_DEPTH_M:g} m for this diagnostic)"
       if _n_dry_elements else ")"))
(OUT / "report.json").write_text(json.dumps(reports, indent=1, default=float))
if reports["qa"]["n_introduced"]:
    raise SystemExit(
        f"the patch introduces {reports['qa']['n_introduced']} QA violation(s) "
        "the base did not have: "
        + "; ".join(f"{v['check']} at {v['kind']} "
                    f"{v.get('id', v.get('elements'))}"
                    for v in reports["qa"]["introduced"][:6]))
say("done")
