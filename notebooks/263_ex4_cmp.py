import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.rcParams["agg.path.chunksize"] = 100000
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
sys.path.insert(0, os.path.expanduser(
    "~/Github/fvcom-mesh-tools/src"))
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid
from oceanmesh import DEM

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex4prvi")
p = np.load(OUT/"p.npy"); t = np.load(OUT/"t.npy")
m = loadmat(str(OUT/"matlab_ex4.mat"))
gp, gt = m["pc"], m["tc"].astype(int) - 1
print(f"ours NP={len(p):,} ml NP={len(gp):,}", flush=True)

nest1 = np.array([[-100, 5], [-53, 5], [-53, 52.5],
                  [-100, 52.5], [-100, 5]])
boxes = [("nest1", nest1, "green")]
for nm, dem_p, col in (
    ("PR", DS/"PR_1arcsec/pr_1s.nc", "magenta"),
    ("USVI", DS/"USVI_1arcsec/usvi_1_mhw_2014.nc", "orange"),
    ("SJ", DS/"SanJuan_1-9arcsec/san_juan_19_prvd02_2015.nc",
     "red"),
):
    bb = DEM(str(dem_p)).bbox
    ring = np.array([[bb[0], bb[2]], [bb[1], bb[2]],
                     [bb[1], bb[3]], [bb[0], bb[3]],
                     [bb[0], bb[2]]])
    boxes.append((nm, ring, col))
    print(nm, bb, flush=True)

views = {
    "full": (-100.5, -52.5, 4.5, 53, 0.03, GridRef(-100, 5, -53, 52.5, 2, 2), 150),
    "prvi": (-68.3, -64.7, 16.8, 19.2, 0.08, GridRef(-68.25, 16.8, -64.7, 19.2, 0.25, 0.25), 170),
    "san_juan": (-66.24, -65.90, 18.36, 18.56, 0.25, GridRef(-66.24, 18.36, -65.90, 18.56, 0.02, 0.02), 170),
}
for name, (x0, x1, y0, y1, lw, G, dpi) in views.items():
    asp = 1/np.cos(np.deg2rad(0.5*(y0+y1)))
    h = 8.5*(y1-y0)*asp/(x1-x0)
    fig, axes = plt.subplots(1, 2, figsize=(17, min(max(h, 5), 12)))
    for ax, (pp, tt, ttl, n_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        sel = tt
        if name != "full":
            c = pp[tt].mean(axis=1)
            keep = ((c[:, 0] > x0-0.1) & (c[:, 0] < x1+0.1)
                    & (c[:, 1] > y0-0.1) & (c[:, 1] < y1+0.1))
            sel = tt[keep]
        ax.triplot(pp[:, 0], pp[:, 1], sel, lw=lw,
                   color="steelblue", rasterized=True)
        for nm, ring, col in boxes:
            ax.plot(ring[:, 0], ring[:, 1], "-", color=col,
                    lw=1.6, label=f"{nm} boubox")
        ax.legend(loc="lower left", fontsize=7)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        add_atlas_grid(ax, crs="EPSG:4326", grid=G)
        ax.set_title(f"{ttl} NP={n_:,}")
    fig.suptitle(f"Example_4 PRVI (4 nests) — {name}")
    fig.tight_layout()
    fig.savefig(OUT/f"cmp_{name}.png", dpi=dpi,
                bbox_inches="tight")
    plt.close(fig)
    print("saved", name, flush=True)
