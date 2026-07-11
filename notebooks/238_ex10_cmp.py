import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

G = GridRef(-1.0, -1.0, 2.0, 2.0, 0.25, 0.25)
OUT = Path("outputs/om2d_examples/ex10")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
m = loadmat(str(OUT / "matlab_ex10.mat"))
gp, gt = m["pc"], m["tc"].astype(int) - 1
inner_box = np.array(
    [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], float)
outer_box = np.array(
    [[-1, -1], [2, -1], [2, 2], [-1, 2], [-1, -1]], float)
views = {
    "full": (-1.05, 2.05, -1.05, 2.05, 0.25),
    "nest_edge": (-0.35, 0.55, 0.35, 1.25, 0.4),
    "inner_corner": (-0.2, 0.35, -0.25, 0.3, 0.5),
}
for name, (x0, x1, y0, y1, lw) in views.items():
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.5))
    for ax, (pp, tt, ttl, n_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw, color="steelblue")
        ax.plot(inner_box[:, 0], inner_box[:, 1], "m-", lw=2.0,
                label="inner nest bbox [0,1]$^2$")
        ax.plot(outer_box[:, 0], outer_box[:, 1], "g-", lw=2.0,
                label="outer bbox [-1,2]$^2$")
        ax.legend(loc="lower left", fontsize=8)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(1)
        add_atlas_grid(ax, crs="EPSG:4326", grid=G)
        ax.set_title(f"{ttl} NP={n_:,}")
    fig.suptitle(f"Example_10 Multiscale Smoother — {name}")
    fig.tight_layout()
    fig.savefig(OUT / f"cmp_{name}.png", dpi=170,
                bbox_inches="tight")
    print("saved", name, flush=True)
