# Figures for a local refinement: what the patch replaced, what it left alone.
#
#   python notebooks/421_local_refine_map.py outputs/refine_futtsu_nori
#
# Every panel draws the mesh itself.  The point of this operation is a claim
# about which elements changed, and a filled contour cannot show that -- only
# the edges can.
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
from fvcom_mesh_tools.qa import run_qa

OUT = Path(sys.argv[1] if len(sys.argv) > 1
           else "outputs/refine_futtsu_nori").resolve()
rep = json.loads((OUT / "report.json").read_text())
base = read_fort14(rep["base_mesh"])
patched = read_fort14(next(OUT.glob("*.14")))
nr = rep["stitch"]["n_nodes_retained"]
ner = rep["selection"]["n_elements_retained"]
pf = rep["preflight"][0]

to_m = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
cx, cy = to_m.transform(*json.loads(os.environ.get("LR_CENTRE", "[139.7881, 35.3228]")))
radius = float(os.environ.get("LR_RADIUS", 300.0))
core = shapely.Point(cx, cy).buffer(radius, quad_segs=128)
core_xy = np.asarray(core.exterior.coords)

tb = Triangulation(base.nodes[:, 0], base.nodes[:, 1], base.elements)
tp = Triangulation(patched.nodes[:, 0], patched.nodes[:, 1], patched.elements)
is_patch = np.arange(patched.n_elements) >= ner

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3)

# --- the whole bay, with the site marked --------------------------------
ax = fig.add_subplot(gs[:, 0])
ax.triplot(tp, color="0.7", lw=0.25)
ax.triplot(Triangulation(patched.nodes[:, 0], patched.nodes[:, 1],
                         patched.elements[is_patch]), color="tab:blue", lw=0.3)
ax.plot(cx, cy, "o", ms=9, mfc="none", mec="tab:red", mew=2)
ax.set_aspect(1)
ax.set_title(f"{Path(rep['base_mesh']).name}\n"
             f"NE {base.n_elements:,} -> {patched.n_elements:,}", fontsize=9)
ax.set_xticks([])
ax.set_yticks([])

# --- the hole: base under, patched over ---------------------------------
W = float(pf["transition_m"]) + radius + 900.0
for cell, title, zoom in [
        (gs[0, 1], "cut and fill (base in grey)", W),
        (gs[0, 2], "the seam", 0.45 * W),
        (gs[1, 1], "the core at the target size", 3.0 * radius)]:
    ax = fig.add_subplot(cell)
    ax.triplot(tb, color="0.75", lw=0.9)
    ax.triplot(tp, color="tab:blue", lw=0.45)
    ax.plot(patched.nodes[nr:, 0], patched.nodes[nr:, 1], ".",
            ms=1.6, color="tab:orange")
    ax.plot(core_xy[:, 0], core_xy[:, 1], color="tab:red", lw=1.4)
    ax.set_xlim(cx - zoom, cx + zoom)
    ax.set_ylim(cy - zoom, cy + zoom)
    ax.set_aspect(1)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])

# --- achieved edge length against distance from the core ----------------
ax = fig.add_subplot(gs[1, 2])
e = np.unique(np.sort(np.vstack([patched.elements[:, [0, 1]],
                                 patched.elements[:, [1, 2]],
                                 patched.elements[:, [2, 0]]]), axis=1), axis=0)
mid = 0.5 * (patched.nodes[e[:, 0], :2] + patched.nodes[e[:, 1], :2])
r = np.linalg.norm(mid - [cx, cy], axis=1)
L = np.linalg.norm(patched.nodes[e[:, 0], :2] - patched.nodes[e[:, 1], :2], axis=1)
eb = np.unique(np.sort(np.vstack([base.elements[:, [0, 1]], base.elements[:, [1, 2]],
                                  base.elements[:, [2, 0]]]), axis=1), axis=0)
midb = 0.5 * (base.nodes[eb[:, 0], :2] + base.nodes[eb[:, 1], :2])
rb = np.linalg.norm(midb - [cx, cy], axis=1)
Lb = np.linalg.norm(base.nodes[eb[:, 0], :2] - base.nodes[eb[:, 1], :2], axis=1)
keep = rb < W
ax.plot(rb[keep], Lb[keep], ".", ms=2, color="0.7", label="base")
ax.plot(r[r < W], L[r < W], ".", ms=2, color="tab:blue", label="patched")
ax.axvline(radius, color="tab:red", lw=1.2)
ax.axvline(radius + pf["transition_m"], color="tab:red", lw=1.0, ls="--")
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
    f"{pf['transition_m']:.0f} m  |  QA "
    f"{qa.n_gate_total - qa.n_gate_failed}/{qa.n_gate_total}  |  frozen nodes "
    f"{ver['n_frozen_nodes']:,}, moved {ver['n_frozen_moved']}, retained faces "
    f"missing {ver['n_retained_faces_missing']}, interface splits "
    f"{ver['n_interface_segments_split']}  |  dt {pf['dt_s']:.2f} s "
    f"(expected {pf['dt_expected_s']:g} s)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
png = OUT / f"421_{OUT.name}.png"
fig.savefig(png, dpi=130)
print(f"[map] wrote {png}")
