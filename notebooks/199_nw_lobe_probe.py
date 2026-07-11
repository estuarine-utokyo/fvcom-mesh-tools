# Why is the C5 NW wetland lobe unmeshed? SDF sign + sizing there.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = Path("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 * DEG)
sdf = om.signed_distance_function(shore)
print(f"[nw] inner pts total={len(np.asarray(shore.inner)):,}",
      flush=True)

x0, x1, y0, y1 = -74.13, -74.07, 40.755, 40.795
xg = np.linspace(x0, x1, 800)
yg = np.linspace(y0, y1, 600)
X, Y = np.meshgrid(xg, yg)
d = sdf.eval(np.column_stack([X.ravel(), Y.ravel()]))
D = d.reshape(X.shape)
frac = float((D < 0).mean())
print(f"[nw] water fraction in NW-lobe window: {frac:.3f}", flush=True)
# feature sizing value at lobe centre
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3*DEG)
pts = np.array([[-74.105, 40.780], [-74.10, 40.777],
                [-74.095, 40.770], [-74.11, 40.765]])
fv = feat.eval(pts) / DEG
sv = sdf.eval(pts)
for q, f_, s_ in zip(pts, fv, sv):
    print(f"[nw] point ({q[0]:.3f},{q[1]:.3f}): fh={f_:.0f} m "
          f"sdf={'water' if s_ < 0 else 'LAND'} ({s_:.2e})",
          flush=True)
fig, ax = plt.subplots(figsize=(12, 8))
ax.pcolormesh(X, Y, (D < 0), cmap="coolwarm_r", shading="auto")
inner = np.asarray(shore.inner)
main = np.asarray(shore.mainland)
if inner.size:
    ax.plot(inner[:, 0], inner[:, 1], "m-", lw=0.6)
if main.size:
    ax.plot(main[:, 0], main[:, 1], "g-", lw=0.6)
ax.set_aspect(1 / np.cos(np.deg2rad(40.775)))
ax.set_title("NW lobe SDF sign (blue=water); inner m / mainland g")
fig.savefig(OUT / "nw_lobe_sign.png", dpi=160, bbox_inches="tight")
print("[nw] saved", flush=True)
