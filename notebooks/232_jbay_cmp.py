# JBAY final visual comparison: our seed3 (all-pass, V=2.075) vs
# MATLAB golden, atlas grid A-K x 1-9 (0.02 deg).
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

JB_GRID = GridRef(-73.97, 40.50, -73.75, 40.68, 0.02, 0.02)
OUT = Path("outputs/om2d_examples/jbay")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
f14 = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/5_JBAY_w_weirs_mesh.14"))
with open(f14) as f:
    f.readline()
    ne, npn = map(int, f.readline().split()[:2])
    gp = np.empty((npn, 2))
    for i in range(npn):
        q = f.readline().split()
        gp[i] = float(q[1]), float(q[2])
    gt = np.empty((ne, 3), dtype=int)
    for i in range(ne):
        q = f.readline().split()
        gt[i] = int(q[2]), int(q[3]), int(q[4])
gt -= 1
print(f"ours NP={len(p):,} golden NP={npn:,}", flush=True)

views = {
    "full": (-73.97, -73.75, 40.50, 40.68, 0.15),
    "weir_west": (-73.95, -73.92, 40.545, 40.582, 0.5),
    "weir_east": (-73.895, -73.868, 40.558, 40.588, 0.5),
    "marsh": (-73.87, -73.82, 40.60, 40.64, 0.35),
}
for name, (x0, x1, y0, y1, lw) in views.items():
    asp = 1 / np.cos(np.deg2rad(0.5 * (y0 + y1)))
    h = 8.5 * (y1 - y0) * asp / (x1 - x0)
    fig, axes = plt.subplots(1, 2, figsize=(17, max(h, 5)))
    for ax, (pp, tt, title, n_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python, seed3)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw, color="steelblue")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        add_atlas_grid(ax, crs="EPSG:4326", grid=JB_GRID)
        ax.set_title(f"{title} NP={n_:,}")
    fig.suptitle(
        f"Example_5b JBAY w/ weirs (h0=15 m) — {name}")
    fig.tight_layout()
    fig.savefig(OUT / f"cmp_{name}.png", dpi=170,
                bbox_inches="tight")
    print("saved", name, flush=True)
