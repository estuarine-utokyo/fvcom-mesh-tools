# Figure for a local refinement: the finished mesh, with the part that changed.
#
#   python notebooks/421_local_refine_map.py outputs/refine_futtsu_nori
#
# One mesh, not two. The claim this operation makes is about which elements
# were replaced, so the figure draws the DELIVERED mesh and colours those
# elements red; everything else is drawn exactly as it is, because it is
# exactly as it was. No base underlay and no node markers: both hid the mesh
# they were meant to explain (owner, 2026-09-22).
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapely
from matplotlib.tri import Triangulation
from pyproj import Transformer

from fvcom_mesh_tools.io.fort14 import read_fort14
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case
from fvcom_mesh_tools.qa import run_qa

KEPT = "0.45"
CHANGED = "tab:red"

OUT = Path(sys.argv[1] if len(sys.argv) > 1
           else "outputs/refine_futtsu_nori").resolve()
rep = json.loads((OUT / "report.json").read_text())
if Path(rep["base_mesh"]).suffix == ".dat":
    base = read_fvcom_case(rep["base_mesh"], rep["base_depth"], rep.get("base_obc"))
else:
    base = read_fort14(rep["base_mesh"])
# By name from the report, not by globbing: an output directory reused for a
# different base keeps the old .14 beside the new one, and a figure drawn
# from the wrong mesh looks entirely plausible.
patched = read_fort14(rep["mesh"])
ner = rep["selection"]["n_elements_retained"]
pf = rep["preflight"][0]

to_m = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
cx, cy = to_m.transform(*json.loads(os.environ.get("LR_CENTRE", "[139.7881, 35.3228]")))
radius = float(os.environ.get("LR_RADIUS", 300.0))
core_xy = np.asarray(shapely.Point(cx, cy).buffer(radius, quad_segs=128).exterior.coords)

# Retained elements come first and in base order, so the split is exact.
kept = Triangulation(patched.nodes[:, 0], patched.nodes[:, 1],
                     patched.elements[:ner])
changed = Triangulation(patched.nodes[:, 0], patched.nodes[:, 1],
                        patched.elements[ner:])

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3)

# --- the whole bay --------------------------------------------------------
ax = fig.add_subplot(gs[:, 0])
ax.triplot(kept, color=KEPT, lw=0.25)
ax.triplot(changed, color=CHANGED, lw=0.3)
ax.set_aspect(1)
ax.set_title(f"{Path(rep['base_mesh']).name}\n"
             f"NE {base.n_elements:,} -> {patched.n_elements:,}", fontsize=9)
ax.set_xticks([])
ax.set_yticks([])

# --- three zooms ----------------------------------------------------------
W = float(pf["transition_m"]) + radius + 900.0
for cell, title, zoom, lw in [
        (gs[0, 1], "the patch", W, 0.45),
        (gs[0, 2], "the seam", 0.45 * W, 0.6),
        (gs[1, 1], "the core at the target size", 3.0 * radius, 0.7)]:
    ax = fig.add_subplot(cell)
    ax.triplot(kept, color=KEPT, lw=lw)
    ax.triplot(changed, color=CHANGED, lw=lw)
    ax.plot(core_xy[:, 0], core_xy[:, 1], color="k", lw=1.0, ls="--")
    ax.set_xlim(cx - zoom, cx + zoom)
    ax.set_ylim(cy - zoom, cy + zoom)
    ax.set_aspect(1)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])

# --- achieved edge length against distance from the core ------------------
# The delivered mesh only, split the same way: an edge is "replaced" when a
# replaced element uses it.
ax = fig.add_subplot(gs[1, 2])


def edges(tri):
    return np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                        tri[:, [2, 0]]]), axis=1), axis=0)


all_e = edges(patched.elements)
new_e = {tuple(x) for x in edges(patched.elements[ner:]).tolist()}
is_new = np.array([tuple(x) in new_e for x in all_e.tolist()], dtype=bool)
mid = 0.5 * (patched.nodes[all_e[:, 0], :2] + patched.nodes[all_e[:, 1], :2])
r = np.linalg.norm(mid - [cx, cy], axis=1)
length = np.linalg.norm(patched.nodes[all_e[:, 0], :2]
                        - patched.nodes[all_e[:, 1], :2], axis=1)
keep = r < W
ax.plot(r[keep & ~is_new], length[keep & ~is_new], ".", ms=2, color=KEPT,
        label="unchanged")
ax.plot(r[keep & is_new], length[keep & is_new], ".", ms=2, color=CHANGED,
        label="replaced")
ax.axvline(radius, color="k", lw=1.0, ls="--")
ax.axvline(radius + pf["transition_m"], color="k", lw=0.8, ls=":")
ax.axhline(pf["target_h_m"], color="k", lw=0.8, ls=":")
ax.set_xlabel("distance from the core centre (m)", fontsize=8)
ax.set_ylabel("edge length (m)", fontsize=8)
ax.set_title("achieved size by distance", fontsize=9)
ax.legend(fontsize=7, markerscale=3)
ax.tick_params(labelsize=7)

qa = run_qa(patched, name="patched")
ver = rep["verify"]
fig.suptitle(
    f"{OUT.name}: {pf['name']} target {pf['target_h_m']:g} m, transition "
    f"{pf['transition_m']:.0f} m  |  red = replaced "
    f"({patched.n_elements - ner:,} of {patched.n_elements:,} elements)  |  "
    f"frozen nodes {ver['n_frozen_nodes']:,}, moved {ver['n_frozen_moved']}, "
    f"depth change {ver['max_frozen_depth_change_m']:g} m  |  dt "
    f"{rep['achieved']['dt_min_s']:.2f} s  |  QA "
    f"{qa.n_gate_total - qa.n_gate_failed}/{qa.n_gate_total}, "
    f"{rep['qa']['n_introduced']} introduced by the patch", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
png = OUT / f"421_{OUT.name}.png"
fig.savefig(png, dpi=130)
print(f"[map] wrote {png}")
