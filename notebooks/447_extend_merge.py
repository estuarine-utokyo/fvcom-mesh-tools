# Finish the outer mesh, merge it onto the base, give the new nodes depths,
# write the FVCOM case and run QA.
#
#   python notebooks/447_extend_merge.py RECIPE GENDIR OUTDIR
#
# Stage 2 of notebooks/445_extend_mesh.py. The base comes through unchanged
# (nodes, elements and depths, checked bit for bit); only the new nodes get
# depths, from the recipe's bathymetry sources, floored, optionally capped,
# and r-factor limited with the base depths held fixed.
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import geopandas as gpd  # noqa: E402
import shapely  # noqa: E402
from pyproj import Transformer  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from fvcom_mesh_tools.algorithms.obc_finish import finish_obc_mesh  # noqa: E402
from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local  # noqa: E402
from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast  # noqa: E402
from fvcom_mesh_tools.dem.m7001 import node_edges  # noqa: E402
from fvcom_mesh_tools.dem.sources import non_tp_count, sample  # noqa: E402
from fvcom_mesh_tools.extend import (  # noqa: E402
    check_island_holes,
    check_no_overlap,
    land_segments,
    merge_outer,
    rfactor_smooth_free,
    round_depths_inside,
    trim_lone_corners,
    verify_frozen_base,
)
from fvcom_mesh_tools.extend_recipe import check_expected, load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.io.fort14 import read_fort14, write_fort14  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case  # noqa: E402
from fvcom_mesh_tools.patch import improve_patch  # noqa: E402
from fvcom_mesh_tools.qa import run_qa  # noqa: E402

T0 = time.time()
MESH_EPSG = 32654


def say(msg):
    print(f"[merge] {msg} +{time.time() - T0:.0f}s", flush=True)


recipe = load_extend_recipe(sys.argv[1])
check_expected(recipe)          # the recipe the driver recorded (review round 5 F3)
GEN = Path(sys.argv[2]).resolve()
OUT = Path(sys.argv[3]).resolve()
S, D = recipe["settings"], recipe["depths"]
CASE = recipe["case"]
b = Path(recipe["base"]) / recipe["base_case"]
base = read_fvcom_case(f"{b}_grd.dat", f"{b}_dep.dat", f"{b}_obc.dat")
IB = np.asarray(base.open_boundaries[0], np.int64)

# --------------------------------------------------------------- finishing
outer = read_fort14(GEN / "outer_utm.14")
land_utm = shapely.union_all(list(gpd.read_file(GEN / "land_with_base.shp")
                                  .to_crs(MESH_EPSG).geometry))
outer, info = finish_obc_mesh(outer, seed=int(S["fin_seed"]), land_union=land_utm,
                              one_wide="forbid")
say("finish: " + json.dumps({k: v for k, v in info.items() if not isinstance(v, (list, dict))},
                            default=str)[:600])
fixed = np.concatenate([np.asarray(c, int) for c in outer.open_boundaries])
cf = fit_boundary_to_coast(outer.nodes, outer.elements, land_utm, fixed=fixed,
                           depths=outer.depths)
outer.nodes[:, :2] = cf.nodes[:, :2]
say(f"coast fit: {cf.summary()}")

# the interface must still be the base's nodes and edges, exactly
d, io_ = cKDTree(outer.nodes[:, :2]).query(base.nodes[IB, :2])
if d.max() > 1e-6:
    raise SystemExit(f"interface node(s) moved in finishing (up to {d.max():.3g} m)")
edges = {frozenset(e) for e in np.vstack([outer.elements[:, [0, 1]], outer.elements[:, [1, 2]],
                                           outer.elements[:, [2, 0]]]).tolist()}
lost = [i for i in range(len(io_) - 1) if frozenset((int(io_[i]), int(io_[i + 1]))) not in edges]
if lost:
    raise SystemExit(f"interface edge(s) {lost} were split or flipped in finishing")
open_new = [np.asarray(c, int) for c in outer.open_boundaries
            if not set(np.asarray(c, int).tolist()) <= set(io_.tolist())]
if len(open_new) != 1:
    raise SystemExit(f"expected one new open boundary after finishing, found {len(open_new)}")

# ------------------------------------------------------------------ merge
merged = merge_outer(base, outer.nodes[:, :2], outer.elements, io_, IB, open_new[0])
contract = verify_frozen_base(merged, base, IB)
say("frozen base: " + json.dumps(contract))

# ----------------------------------------------------------------- repair
# The sixth build left 9 QA failures, all on new elements: a lone corner at a
# cape, three coastal angles under 30 deg, three area jumps, one valence of 9
# and one open-boundary node at 27 deg. They are repaired here, restricted so
# that no base node moves and no base element changes.
NB, EB = base.n_nodes, base.n_elements
obc_nodes = np.asarray(merged.open_boundaries[0], np.int64)
# a cape tip one element wide: not resolved (resolution principle); bisecting
# it instead (walls.open_lone_corners) left two 27.5 deg slivers (7th build)
elems, mutable, lone = trim_lone_corners(merged.elements, np.arange(merged.n_elements) >= EB,
                                         keep_nodes=obc_nodes)
used = np.zeros(merged.n_nodes, bool)
used[elems.ravel()] = True
if not used[:NB].all():
    raise SystemExit("trimming a lone corner orphaned a base node")
renum = np.cumsum(used) - 1
merged = replace(merged, nodes=merged.nodes[used], elements=renum[elems],
                 depths=merged.depths[used], open_boundaries=[renum[obc_nodes]])
obc_nodes = renum[obc_nodes]
nodes, elems = merged.nodes, merged.elements
_u, _c = np.unique(np.sort(np.vstack([elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]]),
                           axis=1), axis=0, return_counts=True)
on_boundary = np.zeros(len(nodes), bool)
on_boundary[np.unique(_u[_c == 1])] = True
# the elements on the open boundary and the nodes just inside it are left as
# finishing made them: flipping or moving them there broke the perpendicular
# partner of 8 open-boundary nodes (39 deg, 7th build)
at_obc = np.isin(elems, obc_nodes).any(axis=1)
next_to_obc = np.zeros(len(nodes), bool)
next_to_obc[elems[at_obc].ravel()] = True
movable = (np.arange(len(nodes)) >= NB) & ~on_boundary & ~next_to_obc
mutable = mutable & ~at_obc
nodes, elems, imp = improve_patch(nodes, elems, movable, mutable, only_below=1.15)
if imp["min_angle_deg"] < 30.0 or imp["max_angle_deg"] > 130.0:
    nodes, elems, imp2 = improve_patch(nodes, elems, movable, mutable, only_below=1.15,
                                       soft=True)
    imp = {**imp2, "n_flips": imp["n_flips"] + imp2["n_flips"],
           "n_moves": imp["n_moves"] + imp2["n_moves"], "soft_pass": True}
merged = replace(merged, nodes=nodes, elements=elems)
merged, perp = align_open_boundary_local(merged)
verify_frozen_base(merged, base, IB)
repair = {"lone_corners": lone, "improve": imp,
          "perpendicularity": {k: v for k, v in perp.items() if not isinstance(v, list)},
          "perp_remaining": len(perp.get("remaining", []))}
say("repair: " + json.dumps(repair, default=str)[:700])

# ----------------------------------------------------------------- depths
to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
lon, lat = (np.asarray(v) for v in to_ll.transform(merged.nodes[:, 0], merged.nodes[:, 1]))
new = np.arange(merged.n_nodes) >= base.n_nodes
raw, which = sample(recipe["bathymetry"]["depths"], lon[new], lat[new])
if np.isnan(raw).any():
    raise SystemExit(f"{int(np.isnan(raw).sum())} new node(s) outside every bathymetry source")
h = merged.depths.copy()
h[new] = np.maximum(raw, D["min_m"])
if D["max_m"] is not None:
    h[new] = np.minimum(h[new], D["max_m"])
ei, ej = node_edges(merged.elements)
# the cap is applied inside the limiter, so the result meets both (review F3)
h, iters, _ = rfactor_smooth_free(h, ei, ej, new, rmax=D["rfactor"], hmin=D["min_m"],
                                  hmax=D["max_m"])
# six decimals, inside the recipe's bounds (review round 10 F12); the r
# check below is on these written depths
h[new] = round_depths_inside(h[new], D["min_m"], D["max_m"])
# r over every edge with a new end, on the depths that are written
touch = new[ei] | new[ej]
r_after = float((np.abs(h[ei] - h[ej]) / (h[ei] + h[ej]))[touch].max()) if touch.any() else 0.0
if r_after > D["rfactor"] + 1e-6:
    raise SystemExit(f"final r-factor {r_after:.4f} exceeds {D['rfactor']}")
merged.depths = h
verify_frozen_base(merged, base, IB)
names = recipe["bathymetry"]["depths"]
depth_report = {
    "sources": names,
    "nodes_per_source": {n: int((which == k).sum()) for k, n in enumerate(names)},
    "raw_min_m": float(raw.min()), "raw_max_m": float(raw.max()),
    "new_min_m": float(h[new].min()), "new_max_m": float(h[new].max()),
    "rfactor_iterations": iters, "r_max_on_free_edges": r_after,
    "nodes_not_on_tp": non_tp_count(names, which)[0],
}
say("depths: " + json.dumps(depth_report))

# ----------------------------------------------------------------- export
# no hole in the new sea without land in it (review round 11 F7), checked
# before anything is written (round 12 F2)
islands = check_island_holes(merged, land_utm, base.n_nodes)
say("islands: " + json.dumps(islands))
merged.land_boundaries = land_segments(merged.elements, merged.open_boundaries)
written = export_fvcom_case(merged, OUT, CASE, cor=lat, obc_depth_control=False)
write_fort14(merged, OUT / f"{CASE}.14")
# the fort.14 copy carries the same mesh, bit for bit (review round 3 F8)
_f14 = read_fort14(OUT / f"{CASE}.14")
if not (np.array_equal(_f14.nodes[:, :2], merged.nodes[:, :2])
        and np.array_equal(_f14.elements, merged.elements)
        and np.array_equal(_f14.depths, merged.depths)):
    raise SystemExit(f"{CASE}.14 does not carry the mesh that was built")
back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
if not np.array_equal(back.depths, merged.depths):
    raise SystemExit("the written case does not carry the depths that were built")
# the frozen contract on what was written, not only on what was built: the
# writer rounds coordinates (review F9)
contract = verify_frozen_base(back, base, IB)
overlap = check_no_overlap(back, base.n_elements)      # anywhere, not only the seam (r2 F9)

# whether the extension limits the time step (review F4; a warning, not a
# failure -- owner, 2026-10-01): the smallest edge / sqrt(g H) over the new
# elements against the base's
def _dt_allow(mesh, elems):
    xy = mesh.nodes[elems, :2]
    edge = np.linalg.norm(xy - np.roll(xy, 1, axis=1), axis=2).min(axis=1)
    return edge / np.sqrt(9.81 * np.maximum(mesh.depths[elems].max(axis=1), 1e-9))


dt_base = float(_dt_allow(back, back.elements[:base.n_elements]).min())
dt_new = float(_dt_allow(back, back.elements[base.n_elements:]).min())
say(f"time-step allowance: base {dt_base:.2f} s, new elements {dt_new:.2f} s")
qa = run_qa(back, name=CASE, path=written["grd"], max_offenders=10_000)
(OUT / f"{CASE}_qa.json").write_text(json.dumps(qa.to_dict(), indent=1, default=float))
say(f"QA {qa.n_gate_total - qa.n_gate_failed}/{qa.n_gate_total}")
for c in qa.checks:
    if c.status == "fail":
        say(f"  FAIL {c.check_id} {c.requirement} | {c.observed}")
problems = []
if qa.n_gate_failed:
    problems.append(f"QA {qa.n_gate_failed} gate(s) failed")
warnings_ = []
if dt_new < dt_base:
    # reported, not refused: the time step is settled in the depth stage
    # (maximum depth, smoothing), the mesh is made from the real depths
    # (owner, 2026-10-01)
    warnings_.append(f"new elements limit the time step ({dt_new:.2f} s < base "
                     f"{dt_base:.2f} s, raw depths)")
(OUT / "merge.json").write_text(json.dumps({
    "finish": {k: v for k, v in info.items() if not isinstance(v, (list, dict))},
    "coast_fit": cf.to_dict(), "frozen_base": contract, "repair": repair,
    "islands": islands,
    "depths": depth_report,
    "qa": {"n_gate_total": qa.n_gate_total, "n_gate_failed": qa.n_gate_failed},
    "dt_allowance_s": {"base": dt_base, "new": dt_new}, "problems": problems,
    "warnings": warnings_,
    "overlap": overlap,
    "n_nodes": merged.n_nodes, "n_elements": merged.n_elements,
    "n_open_boundary_nodes": int(len(merged.open_boundaries[0])),
}, indent=1, default=str))
say(f"wrote {', '.join(sorted(written))} + {CASE}.14 in {OUT}")
# the reports stay for diagnosis, but a failed build must not exit as a
# success (review F2)
for w_ in warnings_:
    say("WARNING: " + w_)
if problems:
    raise SystemExit("the build failed: " + "; ".join(problems))
