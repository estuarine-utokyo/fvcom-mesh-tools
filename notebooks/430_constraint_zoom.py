# The constraints and the delivered mesh around one point.
#
#   python notebooks/430_constraint_zoom.py <refine out dir> <x> <y> <half_m>
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
x, y, half = (float(a) for a in sys.argv[2:5])
z = np.load(OUT / "fill_constraints.npz")
pf, eg, nr = z["pfix"], z["egfix"], int(z["n_rim"])
m = read_fort14(next(OUT.glob("*.14")))
near = np.flatnonzero(np.hypot(pf[:, 0] - x, pf[:, 1] - y) < half)
print(f"fixed points within {half:g} m: {len(near)}")
for k in near:
    kind = "rim" if k < nr else "wall"
    print(f"  pfix {k:5d} {kind:4s} ({pf[k, 0]:.2f}, {pf[k, 1]:.2f})")
for a, b in eg:
    if a in near or b in near:
        kind = "rim" if (a < nr and b < nr) else "wall"
        print(f"  egfix {a}-{b} {kind} length {np.linalg.norm(pf[a] - pf[b]):.1f} m")
# the elements under 30 deg in the window, and their smallest angle
xy, tri = m.nodes[:, :2], m.elements
cen = xy[tri].mean(axis=1)
win = np.flatnonzero((np.abs(cen[:, 0] - x) < half) & (np.abs(cen[:, 1] - y) < half))
bad = []
for k in win:
    p3 = xy[tri[k]]
    ang = []
    for i in range(3):
        a, b = p3[(i + 1) % 3] - p3[i], p3[(i + 2) % 3] - p3[i]
        ang.append(np.degrees(np.arccos(np.clip(
            a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))
    if min(ang) < 30:
        bad.append(k)
        print(f"  element {k}: angles {np.round(ang, 1)}, vertices "
              + ", ".join(f"({xy[v, 0]:.1f},{xy[v, 1]:.1f})" for v in tri[k]))
fig, ax = plt.subplots(figsize=(8, 8))
ax.triplot(m.nodes[:, 0], m.nodes[:, 1], m.elements, lw=0.4, color="0.6")
for a, b in eg:
    col = "tab:blue" if (a < nr and b < nr) else "tab:red"
    ax.plot(pf[[a, b], 0], pf[[a, b], 1], color=col, lw=1.8)
ax.plot(pf[near, 0], pf[near, 1], "o", color="k", ms=4)
for k in bad:
    q = xy[np.append(tri[k], tri[k][0])]
    ax.fill(q[:, 0], q[:, 1], color="tab:orange", alpha=0.6)
ax.set_xlim(x - half, x + half)
ax.set_ylim(y - half, y + half)
ax.set_aspect("equal")
fig.savefig(OUT / f"zoom_{x:.0f}_{y:.0f}.png", dpi=150)
print("wrote", OUT / f"zoom_{x:.0f}_{y:.0f}.png")
