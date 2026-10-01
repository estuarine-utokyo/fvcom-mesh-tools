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
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyproj import Transformer  # noqa: E402

from fvcom_mesh_tools.dem.m7001 import node_edges  # noqa: E402
from fvcom_mesh_tools.dem.sources import SOURCES, non_tp_count, sample  # noqa: E402
from fvcom_mesh_tools.extend import land_segments, rfactor_smooth_free  # noqa: E402
from fvcom_mesh_tools.extend_recipe import load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case  # noqa: E402

MESH_EPSG = 32654

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("recipe")
p.add_argument("built_dir", type=Path)
p.add_argument("outdir", type=Path)
p.add_argument("--sources", required=True, help="comma-separated, first covering source wins")
p.add_argument("--case-name", default=None)
a = p.parse_args()

recipe = load_extend_recipe(a.recipe)
D = recipe["depths"]
names = a.sources.split(",")
unknown = [n for n in names if n not in SOURCES]
if unknown:
    raise SystemExit(f"unknown source(s) {unknown}; known {sorted(SOURCES)}")
src_dir = a.built_dir.resolve()
case = recipe["case"]
b = Path(recipe["base"]) / recipe["base_case"]
base = read_fvcom_case(f"{b}_grd.dat", f"{b}_dep.dat", f"{b}_obc.dat")
mesh = read_fvcom_case(src_dir / f"{case}_grd.dat", src_dir / f"{case}_dep.dat",
                       src_dir / f"{case}_obc.dat")
NB = base.n_nodes
if not (np.array_equal(mesh.nodes[:NB, :2], base.nodes[:, :2])
        and np.array_equal(mesh.depths[:NB], base.depths)):
    raise SystemExit("the case does not start with the base nodes and depths")

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
h, iters, r_after = rfactor_smooth_free(h, ei, ej, new, rmax=D["rfactor"], hmin=D["min_m"])
if D["max_m"] is not None:
    h[new] = np.minimum(h[new], D["max_m"])
h[new] = np.round(h[new], 6)
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
}
(a.outdir / "redepth.json").write_text(json.dumps(report, indent=1) + "\n")
print("[453] " + json.dumps(report), flush=True)

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
print(f"[453] wrote {a.outdir}", flush=True)
