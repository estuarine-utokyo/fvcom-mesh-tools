# East-weir window: our fh vs OM2D fh side by side + ratio, with
# the weir rings and shoreline overlaid.
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.weirs import build_weir_geometry

OUT = Path("outputs/om2d_examples/jbay")
DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
WST = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/weirs_struct.mat"))
DEG = 1.0/111e3
w = loadmat(str(WST), squeeze_me=True)["weirs"]
rings, crests, weir_cfg = [], [], []
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
    crests.append(cl)
    weir_cfg.append({"crestline": cl,
                     "min_ele_m": float(wi["min_ele"]),
                     "width_m": float(wi["width"])})
reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  15.0 * DEG)
shore.inner = np.vstack(
    [np.asarray(shore.inner).reshape(-1, 2),
     np.array([[np.nan, np.nan]])] + rings)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3 * DEG)
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=15.0,
    max_edge_length=1e3, gradation=0.15,
    courant={"timestep": 2.0}, weirs=weir_cfg)

with h5py.File(str(OUT / "matlab_jbay_fh.mat"), "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    vals = np.asarray(f["vals"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
x0, x1, y0, y1 = -73.905, -73.862, 40.552, 40.602
xi = np.linspace(x0, x1, 500)
yi = np.linspace(y0, y1, 560)
X, Y = np.meshgrid(xi, yi)
q = np.column_stack([X.ravel(), Y.ravel()])
ours = grid.eval(q).reshape(X.shape) / DEG
from scipy.interpolate import RegularGridInterpolator
Fml = RegularGridInterpolator((xg, yg), vals, bounds_error=False,
                              fill_value=np.nan)
ml = Fml(q).reshape(X.shape)
ml = ml / DEG if np.nanmedian(vals) < 1 else ml
fig, axes = plt.subplots(1, 3, figsize=(21, 8))
for ax, (Z, ttl) in zip(
    axes,
    [(ours, "our fh [m]"), (ml, "OM2D fh [m]"),
     (np.clip(ours/ml, 0, 2), "ratio ours/ml")],
):
    pc = ax.pcolormesh(X, Y, Z, shading="auto",
                       cmap="viridis" if "fh" in ttl else "RdBu_r",
                       vmin=(15 if "fh" in ttl else 0),
                       vmax=(400 if "fh" in ttl else 2))
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ml_ = np.asarray(shore.mainland)
    inn_ = np.asarray(shore.inner)
    for arr, c in ((ml_, "w"), (inn_, "r")):
        if arr.size:
            ax.plot(arr[:, 0], arr[:, 1], c, lw=0.6)
    for cl in crests:
        ax.plot(cl[:, 0], cl[:, 1], "m-", lw=1.5)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1/np.cos(np.deg2rad(40.58)))
    ax.set_title(ttl)
fig.suptitle("East weir window: sizing fields (white=mainland, "
             "red=inner incl. weir, magenta=crest)")
fig.tight_layout()
fig.savefig(OUT / "eastweir_fh.png", dpi=150, bbox_inches="tight")
print("[ew] saved", flush=True)
