import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
sys.path.insert(0, os.path.expanduser(
    "~/Github/fvcom-mesh-tools/src"))
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

OUT = Path("outputs/om2d_examples/ex9tbay")
m = loadmat(str(OUT/"matlab_ex9.mat"))
box1 = np.array([[-83, 27], [-82, 27], [-82, 28.5], [-83, 28.5],
                 [-83, 27]])
box2 = np.array([[-82.8, 27.25], [-82.4, 27.25], [-82.4, 28.25],
                 [-82.8, 28.25], [-82.8, 27.25]])
cases = {
    "min": (np.load(OUT/"p_min.npy"), np.load(OUT/"t_min.npy"),
            m["pc"], m["tc"].astype(int)-1),
    "nomin": (np.load(OUT/"p_nomin.npy"), np.load(OUT/"t_nomin.npy"),
              m["pn"], m["tn"].astype(int)-1),
}
views = {
    "full": (-83.05, -81.95, 26.95, 28.55, 0.2,
             GridRef(-83, 27, -82, 28.5, 0.1, 0.1)),
    "bay": (-82.75, -82.35, 27.5, 28.0, 0.35,
            GridRef(-82.75, 27.5, -82.35, 28.0, 0.05, 0.05)),
}
for tag, (p, t, gp, gt) in cases.items():
    for name, (x0, x1, y0, y1, lw, G) in views.items():
        if tag == "nomin" and name == "bay":
            continue
        asp = 1/np.cos(np.deg2rad(0.5*(y0+y1)))
        h = 8.5*(y1-y0)*asp/(x1-x0)
        fig, axes = plt.subplots(1, 2, figsize=(17, min(max(h,5),12)))
        for ax, (pp, tt, ttl, n_) in zip(
            axes,
            [(p, t, "improved oceanmesh (Python)", len(p)),
             (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
        ):
            ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw,
                       color="steelblue")
            ax.plot(box1[:, 0], box1[:, 1], "g-", lw=1.5,
                    label="nest1 boubox")
            ax.plot(box2[:, 0], box2[:, 1], "m-", lw=1.8,
                    label="nest2 boubox")
            ax.legend(loc="lower left", fontsize=8)
            ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
            ax.set_aspect(asp)
            add_atlas_grid(ax, crs="EPSG:4326", grid=G)
            ax.set_title(f"{ttl} NP={n_:,}")
        fig.suptitle(f"Example_9 TBAY enforceMin={tag} — {name}")
        fig.tight_layout()
        fig.savefig(OUT/f"cmp_{tag}_{name}.png", dpi=170,
                    bbox_inches="tight")
        plt.close(fig)
        print("saved", tag, name, flush=True)
