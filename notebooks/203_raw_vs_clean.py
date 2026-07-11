# Does the om2d cleanup delete whole weakly-connected wetland
# lobes? Compare RAW (cleanup='none') vs cleaned mesh in C5/D5.
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

NY_GRID = GridRef(-74.25, 40.50, -73.75, 41.00, 0.05, 0.05)
DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = Path("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3
min_el, max_el, max_el_ns = 30.0, 1e3, 240.0

t0 = time.time()
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  min_el * DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=max_el * DEG)
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=min_el,
    max_edge_length=max_el, max_edge_length_nearshore=max_el_ns,
    gradation=0.20, courant={"timestep": 2.0},
)
p_raw, t_raw = om.generate_mesh(sdf, grid, max_iter=100, seed=0,
                                cleanup="none")
print(f"[raw] NP={len(p_raw):,} +{time.time()-t0:.0f}s", flush=True)
np.save(OUT / "p_raw.npy", p_raw)
np.save(OUT / "t_raw.npy", t_raw)

p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
x0, x1, y0, y1 = -74.16, -74.06, 40.73, 40.81
asp = 1 / np.cos(np.deg2rad(40.77))
fig, axes = plt.subplots(1, 2, figsize=(17, 9))
for ax, (pp, tt, title) in zip(
    axes,
    [(p_raw, t_raw,
      f"BEFORE cleanup (raw) NP={len(p_raw):,}"),
     (p, t, f"AFTER om2d cleanup NP={len(p):,}")],
):
    ax.triplot(pp[:, 0], pp[:, 1], tt, lw=0.3, color="steelblue")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(asp)
    add_atlas_grid(ax, crs="EPSG:4326", grid=NY_GRID)
    ax.set_title(title)
fig.suptitle("Example_2_NY C4/C5-D5: raw vs cleaned (ours)")
fig.tight_layout()
fig.savefig(OUT / "cmp_c5_raw_vs_clean.png", dpi=170,
            bbox_inches="tight")
print("[raw] saved", flush=True)
