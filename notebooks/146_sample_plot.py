import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID as G

gd = open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')).read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
cells = np.array([[int(w) for w in gd[2+i].split()[1:4]]
                  for i in range(ne)]) - 1
P = np.array([[float(w) for w in gd[2+ne+i].split()[1:3]]
              for i in range(nn)])
obc = [int(l.split()[1]) - 1 for l in open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_obc.dat')
    ).read().strip().split('\n')[1:]]
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(P[:, 0], P[:, 1])
print(f"sample: NP={nn} NE={ne} OBC={len(obc)} nodes", flush=True)
print("extent: lon", round(lon.min(),3), "-", round(lon.max(),3),
      "lat", round(lat.min(),3), "-", round(lat.max(),3), flush=True)
COAST = (139.0, 34.5, 141.3, 36.2)
views = {
    "full": (lon.min()-0.03, lat.min()-0.03, lon.max()+0.03,
             lat.max()+0.03, 0.45),
    "mouth": (139.60, 34.95, 140.05, 35.30, 0.6),
}
for name, (x0, y0, x1, y1, lw) in views.items():
    fig, ax = plt.subplots(figsize=(11, 11*(y1-y0)*1.22/(x1-x0)))
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon, lat, cells, lw=lw, color="steelblue")
    ax.plot(lon[obc], lat[obc], color="red", lw=2.5, zorder=6,
            marker="o", ms=4, label=f"open boundary ({len(obc)} nodes)")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_aspect(1/np.cos(np.deg2rad(0.5*(y0+y1))))
    ax.legend(loc="lower right")
    ax.set_title(f"goto2023 SAMPLE (TokyoBay_grd.dat) NP={nn:,} - {name}")
    fig.savefig(f'outputs/figures/sample_{name}.png', dpi=210,
                bbox_inches="tight")
    print("saved", name, flush=True)
