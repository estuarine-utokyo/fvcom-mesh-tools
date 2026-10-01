"""Give the new nodes of an extended case depths from another bathymetry stack.

    python notebooks/453_redepth_extended.py RECIPE BUILT_DIR OUTDIR \\
        --sources cao_shutochokka_2025,m7001,srtm15plus [--case-name TokyoBayEnshuM7001]

A sensitivity tool: the mesh (nodes, elements, open boundary) of the case
445/447 built from RECIPE into BUILT_DIR is kept, the base nodes keep their depths (checked
bit for bit), and every new node is given its depth again from ``--sources``
with the recipe's floor, cap and r-factor limit -- 447's depth stage, nothing
else. Writes the case, ``redepth.json`` and a map of the change.
"""

from __future__ import annotations

import argparse
import atexit
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyproj import Transformer  # noqa: E402

from fvcom_mesh_tools.dem.m7001 import node_edges  # noqa: E402
from fvcom_mesh_tools.dem.sources import SOURCES, non_tp_count, sample, source_files  # noqa: E402
from fvcom_mesh_tools.extend import (  # noqa: E402
    check_no_overlap,
    land_segments,
    rfactor_smooth_free,
    verify_frozen_base,
)
from fvcom_mesh_tools.extend_recipe import check_case_name, load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case  # noqa: E402
from fvcom_mesh_tools.outdir import reserve  # noqa: E402
from fvcom_mesh_tools.provenance import (  # noqa: E402
    changed_files,
    changed_inventory,
    collect,
)
from fvcom_mesh_tools.qa import run_qa  # noqa: E402

MESH_EPSG = 32654

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("recipe")
p.add_argument("built_dir", type=Path)
p.add_argument("outdir", type=Path)
p.add_argument("--sources", required=True, help="comma-separated, first covering source wins")
p.add_argument("--case-name", default=None)
p.add_argument("--allow-failing-gates", action="store_true",
               help="write a variant that fails QA or the time-step gate (a deliberate "
                    "sensitivity case; the report says so)")
a = p.parse_args()

recipe = load_extend_recipe(a.recipe)
D = recipe["depths"]
names = a.sources.split(",")
unknown = [n for n in names if n not in SOURCES]
if unknown:
    raise SystemExit(f"unknown source(s) {unknown}; known {sorted(SOURCES)}")
src_dir = a.built_dir.resolve()
a.outdir = a.outdir.resolve()
# never write over the source case or another experiment (review F21)
if a.outdir == src_dir or src_dir in a.outdir.parents or a.outdir in src_dir.parents:
    raise SystemExit(f"OUTDIR {a.outdir} overlaps the built case {src_dir}")
if a.case_name is not None:
    check_case_name(a.case_name)                      # review F28
a.outdir = reserve(a.outdir)          # atomically, before any work (review r2 F1)
case = recipe["case"]
# A report is written whatever happens, from the moment the output is
# reserved (review round 3 F6, round 4 F7); a failure before the provenance
# is taken says so ("provenance": null).
STATE = {"done": False, "provenance": None}


def _on_exit():
    if not STATE["done"]:
        (a.outdir / "redepth.json").write_text(json.dumps(
            {"from": str(src_dir / case), "sources": names, "status": "failed",
             "provenance": STATE["provenance"],
             "provenance_complete": STATE["provenance"] is not None},
            indent=1, default=str) + "\n")


atexit.register(_on_exit)
# provenance captured before any work; the inputs read below (recipe, built
# case, base) are hashed again at the end and a change fails the run (round 4 F8)
b = Path(recipe["base"]) / recipe["base_case"]
INPUTS = {"recipe": recipe["recipe_path"],
          "built_case": [str(src_dir / f"{case}_{k}.dat") for k in ("grd", "dep", "obc")],
          "base": [f"{b}_{k}.dat" for k in ("grd", "dep", "obc")]}
CODE = {"fvcom_mesh_tools": str(ROOT / "src" / "fvcom_mesh_tools" / "__init__.py"),
        "driver": __file__}
STATE["provenance"] = PROV = collect(
    code=CODE,
    files={**INPUTS,
           **{f"bathymetry_{k}": [str(q) for q in v] for k, v in source_files(names).items()}})
# the recipe hashed is the one parsed above (review round 5 F3)
if PROV["files"]["recipe"]["sha256"] != recipe["recipe_sha256"]:
    raise SystemExit("the recipe changed between reading and recording it")
base = read_fvcom_case(f"{b}_grd.dat", f"{b}_dep.dat", f"{b}_obc.dat")
mesh = read_fvcom_case(src_dir / f"{case}_grd.dat", src_dir / f"{case}_dep.dat",
                       src_dir / f"{case}_obc.dat")
NB = base.n_nodes
IB = np.asarray(base.open_boundaries[0], np.int64)
# the full frozen contract on the input, connectivity included (review r2 F3)
verify_frozen_base(mesh, base, IB)

to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
lon, lat = (np.asarray(v) for v in to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1]))
new = np.arange(mesh.n_nodes) >= NB
raw, which = sample(names, lon[new], lat[new])
if np.isnan(raw).any():
    raise SystemExit(f"{int(np.isnan(raw).sum())} new node(s) outside every source")
h = mesh.depths.copy()
h[new] = np.maximum(raw, D["min_m"])
if D["max_m"] is not None:
    h[new] = np.minimum(h[new], D["max_m"])
ei, ej = node_edges(mesh.elements)
# the cap is applied inside the limiter, so the result meets both (review F3)
h, iters, _ = rfactor_smooth_free(h, ei, ej, new, rmax=D["rfactor"], hmin=D["min_m"],
                                  hmax=D["max_m"])
h[new] = np.round(h[new], 6)
# r over every edge with a new end, on the depths that are written
touch = new[ei] | new[ej]
r_after = float((np.abs(h[ei] - h[ej]) / (h[ei] + h[ej]))[touch].max()) if touch.any() else 0.0
if r_after > D["rfactor"] + 1e-6:
    raise SystemExit(f"final r-factor {r_after:.4f} exceeds {D['rfactor']}")
if not np.array_equal(h[:NB], base.depths):
    raise SystemExit("a base depth changed")
old = mesh.depths.copy()
mesh.depths = h
mesh.land_boundaries = land_segments(mesh.elements, mesh.open_boundaries)
name = a.case_name or f"{case}_redepth"
a.outdir.mkdir(parents=True, exist_ok=True)
written = export_fvcom_case(mesh, a.outdir, name, cor=lat, obc_depth_control=False)
back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
if not np.array_equal(back.depths, h):
    raise SystemExit("the written case does not carry the depths that were built")
# the acceptance of a build, repeated on the variant (review round 2 F3):
# the frozen base on what was written, no overlap, QA, the time-step gate
verify_frozen_base(back, base, IB)
check_no_overlap(back, base.n_elements)
qa = run_qa(back, name=name, path=written["grd"], max_offenders=10_000)


def _dt_allow(m_, elems):
    xy = m_.nodes[elems, :2]
    edge = np.linalg.norm(xy - np.roll(xy, 1, axis=1), axis=2).min(axis=1)
    return edge / np.sqrt(9.81 * np.maximum(m_.depths[elems].max(axis=1), 1e-9))


dt_base = float(_dt_allow(back, back.elements[:base.n_elements]).min())
dt_new = float(_dt_allow(back, back.elements[base.n_elements:]).min())
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

# the inputs this run read must still be the ones its provenance names
# (round 4 F8); no override accepts a variant built from moving inputs
# every input the provenance names, the bathymetry included (review round 6
# F8), and the datasets listed again (round 7 F9)
changed = sorted(set(changed_files(PROV, list(PROV["files"]))) | set(changed_inventory(
    PROV, {f"bathymetry_{k}": v for k, v in source_files(names).items()})))
# the code that ran must still be the code recorded (round 7 F10)
if collect(code=CODE, libraries=())["code"] != PROV["code"]:
    changed.append("code")

d = h - old
rel = d[new] / np.maximum(old[new], 1.0)
report = {
    "from": str(src_dir / case), "sources": names,
    "nodes_per_source": {n: int((which == k).sum()) for k, n in enumerate(names)},
    "rfactor_iterations": iters, "r_max_on_free_edges": r_after,
    "nodes_not_on_tp": non_tp_count(names, which)[0],
    "change_m": {"mean": float(d[new].mean()), "abs_mean": float(np.abs(d[new]).mean()),
                 "min": float(d[new].min()), "max": float(d[new].max())},
    "relative_change": {"mean": float(rel.mean()), "abs_mean": float(np.abs(rel).mean()),
                        "p05": float(np.percentile(rel, 5)), "p95": float(np.percentile(rel, 95))},
    "case": name, "n_nodes": mesh.n_nodes,
    "depth_controls": D, "allow_failing_gates": bool(a.allow_failing_gates),
    "qa": {"n_gate_total": qa.n_gate_total, "n_gate_failed": qa.n_gate_failed},
    "dt_allowance_s": {"base": dt_base, "new": dt_new}, "problems": problems,
    "warnings": warnings_, "inputs_changed_during_run": changed,
    "status": "failed" if changed else "ok" if not problems else (
        "accepted sensitivity variant" if a.allow_failing_gates else "failed"),
    # what made it, captured before the work (review rounds 2 F18, 3 F6)
    "provenance": PROV,
}


def _write_report():
    (a.outdir / "redepth.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
    STATE["done"] = True
    print("[453] " + json.dumps({k: v for k, v in report.items() if k != "provenance"},
                                default=str), flush=True)


if changed or (problems and not a.allow_failing_gates):
    _write_report()
    if changed:
        raise SystemExit(f"input(s) changed during the run: {changed}")
    raise SystemExit("the variant fails: " + "; ".join(problems)
                     + " (--allow-failing-gates writes it as a sensitivity case)")
# The report says the run succeeded only once the figure is written too; a
# failure in between leaves the exit handler's failure report (round 7 F12).

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

tri = mtri.Triangulation(mesh.nodes[:, 0] / 1e3, mesh.nodes[:, 1] / 1e3, mesh.elements)
edges = np.sort(np.vstack([mesh.elements[:, [0, 1]], mesh.elements[:, [1, 2]],
                           mesh.elements[:, [2, 0]]]), axis=1)
u, cnt = np.unique(edges, axis=0, return_counts=True)
bnd = u[cnt == 1]
obset = set(np.concatenate([np.asarray(o) for o in mesh.open_boundaries]).tolist())
is_open = np.array([i in obset and j in obset for i, j in bnd])
xk, yk = mesh.nodes[:, 0] / 1e3, mesh.nodes[:, 1] / 1e3
segs = np.stack([np.c_[xk[bnd[:, 0]], yk[bnd[:, 0]]], np.c_[xk[bnd[:, 1]], yk[bnd[:, 1]]]], axis=1)
fig, ax = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
for axi, val, cmap, lim, title in (
        (ax[0], h, "viridis_r", (0, 2000), f"depth, {'+'.join(names)} (m)"),
        (ax[1], 100 * np.r_[np.zeros(NB), rel], "RdBu_r", (-30, 30),
         "depth change on new nodes (% of the old depth)")):
    pc = axi.tripcolor(tri, val, cmap=cmap, vmin=lim[0], vmax=lim[1], shading="gouraud")
    axi.add_collection(LineCollection(segs[~is_open], colors="k", linewidths=0.5))
    axi.add_collection(LineCollection(segs[is_open], colors="r", linewidths=1.0))
    axi.set_aspect("equal")
    axi.set_title(title, fontsize=10)
    fig.colorbar(pc, ax=axi, shrink=0.8)
fig.savefig(a.outdir / "redepth.png", dpi=120)
_write_report()
print(f"[453] wrote {a.outdir}", flush=True)
