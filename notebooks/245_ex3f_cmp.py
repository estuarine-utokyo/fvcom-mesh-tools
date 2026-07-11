import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
import sys
sys.path.insert(0, os.path.expanduser(
    "~/Github/fvcom-mesh-tools/src"))
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

OUT = Path("outputs/om2d_examples/ex3full")
p = np.load(OUT/"p.npy"); t = np.load(OUT/"t.npy")
m = loadmat(str(OUT/"matlab_ex3.mat"))
gp, gt = m["pc"], m["tc"].astype(int) - 1
print(f"ours NP={len(p):,} ml NP={len(gp):,}")
quad = np.array([[-71.6, 42.7], [-64, 30], [-80, 24],
                 [-85, 38], [-71.6, 42.7]])
nest2 = np.array([[-74.25, 40.5], [-73.75, 40.55],
                  [-73.75, 41], [-74, 41], [-74.25, 40.5]])
G1 = GridRef(-85, 24, -64, 43, 1.0, 1.0)
G2 = GridRef(-74.25, 40.5, -73.75, 41.0, 0.05, 0.05)
views = {
    "full": (-85.5, -63.5, 23.5, 43.2, 0.1, G1),
    "ny_nest": (-74.3, -73.7, 40.45, 41.05, 0.25, G2),
    "nest_edge_w": (-74.28, -74.05, 40.45, 40.75, 0.4, G2),
}
for name, (x0, x1, y0, y1, lw, G) in views.items():
    asp = 1/np.cos(np.deg2rad(0.5*(y0+y1)))
    h = 8.5*(y1-y0)*asp/(x1-x0)
    fig, axes = plt.subplots(1, 2, figsize=(17, max(h, 5)))
    for ax, (pp, tt, ttl, n_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw,
                   color="steelblue")
        ax.plot(quad[:, 0], quad[:, 1], "g-", lw=1.5,
                label="nest1 boubox")
        ax.plot(nest2[:, 0], nest2[:, 1], "m-", lw=2.0,
                label="nest2 boubox")
        ax.legend(loc="lower left", fontsize=8)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        add_atlas_grid(ax, crs="EPSG:4326", grid=G)
        ax.set_title(f"{ttl} NP={n_:,}")
    fig.suptitle(f"Example_3 FULL (2 nests) — {name}")
    fig.tight_layout()
    fig.savefig(OUT/f"cmp_{name}.png", dpi=170,
                bbox_inches="tight")
    print("saved", name, flush=True)
