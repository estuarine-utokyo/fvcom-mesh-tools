# Side-by-side: improved oceanmesh (Python) vs OM2D golden (MATLAB)
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from fvcom_mesh_tools.io import read_fort14

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
OUT = Path("outputs/om2d_examples/ex1_nz")
g = read_fort14(str(OM2D / "Examples/1_NZ_mesh.14"))
gp, gt = g.nodes[:, :2], g.elements
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")

views = {
    "full": (166.0, 176.0, -48.0, -40.0, 0.15),
    "marlborough": (173.5, 174.6, -41.6, -40.7, 0.5),
    "fiordland": (166.3, 167.3, -46.3, -45.0, 0.5),
    "banks": (172.4, 173.3, -44.2, -43.5, 0.5),
}
for name, (x0, x1, y0, y1, lw) in views.items():
    asp = 1 / np.cos(np.deg2rad(0.5 * (y0 + y1)))
    h = 8 * (y1 - y0) * asp / (x1 - x0)
    fig, axes = plt.subplots(1, 2, figsize=(16, max(h, 4)))
    for ax, (pp, tt, title, np_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw, color="steelblue")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        ax.set_title(f"{title} NP={np_:,}")
    fig.suptitle(f"Example_1_NZ — {name}")
    fig.tight_layout()
    fig.savefig(OUT / f"cmp_{name}.png", dpi=190, bbox_inches="tight")
    print("saved", name, flush=True)
