# Where are the QA failures of a mesh with walls, relative to the walls?
#
#   python notebooks/429_wall_qa_map.py <refinement output dir>
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import boundary_legend, draw_mesh  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
m = read_fort14(next(OUT.glob("*.14")))
qa = json.loads(next(OUT.glob("*_qa.json")).read_text())
xy, tri = m.nodes[:, :2], m.elements
# wall nodes: coincident positions
_, inv, cnt = np.unique(np.round(xy, 6), axis=0, return_inverse=True, return_counts=True)
wall = cnt[inv.ravel()] > 1
e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
u, c = np.unique(e, axis=0, return_counts=True)
bnd_node = np.zeros(len(xy), bool)
bnd_node[np.unique(u[c == 1])] = True
el_on_wall = wall[tri].any(axis=1)
# A wall EDGE is a boundary segment met from both sides, i.e. present twice
# by position.  Judging it by its two nodes being duplicated drew the edge
# to a free tip -- whose tip is NOT duplicated -- as coastline.
_b = u[c == 1]
_key = np.round(np.sort(np.stack([xy[_b[:, 0]], xy[_b[:, 1]]], axis=1)
                        .view([("x", float), ("y", float)]), axis=1)
                .view(float).reshape(len(_b), 4), 3)
_, _kinv, _kc = np.unique(_key, axis=0, return_inverse=True, return_counts=True)
wall_edge = _kc[_kinv.ravel()] == 2

def angles(k):
    p = xy[tri[k]]
    out = []
    for i in range(3):
        a, b = p[(i + 1) % 3] - p[i], p[(i + 2) % 3] - p[i]
        c = a @ b / np.linalg.norm(a) / np.linalg.norm(b)
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return out

fails = [ch for ch in qa["checks"] if ch.get("status") == "fail"]
print("failing checks:", [(ch["check_id"], ch.get("observed")) for ch in fails])
bad_el = []
for ch in fails:
    kinds = Counter()
    for off in ch.get("offender_ids") or ch.get("offenders") or []:
        els = off.get("elements") or ([off["id"]] if off.get("kind") == "element" else [])
        for k in els:
            k = int(k)
            if k < len(tri):
                bad_el.append(k)
                kinds["touches a wall node" if el_on_wall[k] else "no wall node"] += 1
        if off.get("kind") == "node":
            kinds["node on a wall" if wall[int(off["id"])] else "node off walls"] += 1
    print(f"  {ch['check_id']}: {dict(kinds)}")
    if ch["check_id"] != "min_depth_clip":
        offs = ch.get("offender_ids") or ch.get("offenders") or []
        ks = sorted({int(k) for off in offs
                     for k in (off.get("elements")
                               or ([off["id"]] if off.get("kind") == "element" else []))
                     if int(k) < len(tri)})
        for k in ks[:8]:
            cc = xy[tri[k]].mean(axis=0)
            print(f"      element {k} at ({cc[0]:.0f}, {cc[1]:.0f}) min angle "
                  f"{min(angles(k)):.1f}, wall={bool(el_on_wall[k])}")
bad_el = np.unique(bad_el)
wa = [min(angles(k)) for k in bad_el if el_on_wall[k]]
print(f"min-angle of offending wall elements: {np.round(sorted(wa)[:15], 1)}")
cen = xy[tri].mean(axis=1)
c0 = cen[bad_el].mean(axis=0) if len(bad_el) else xy.mean(axis=0)
fig, axes = plt.subplots(1, 2, figsize=(17, 8.5))
for ax, half in zip(axes, (1100.0, 350.0)):
    cx, cy = (393010.0, 3909480.0) if half > 500 else (392500.0, 3909000.0)
    draw_mesh(ax, xy, tri, m.open_boundaries, mesh_lw=0.25)
    ax.plot(cen[bad_el, 0], cen[bad_el, 1], "x", color="tab:orange", ms=8, mew=2)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
axes[0].plot([], [], "x", color="tab:orange", label=f"QA offender ({len(bad_el)})")
boundary_legend(axes[0], loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "walls_qa.png", dpi=150)
print(f"wrote {OUT / 'walls_qa.png'}")
