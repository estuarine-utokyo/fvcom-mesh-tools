# H8 lower-left (Keihin canal, Kawasaki): sample vs ours zoom +
# measured through-connectivity along the canal chain.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from pyproj import Transformer
from shapely.ops import unary_union
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
COSW = float(np.cos(np.deg2rad(35.35)))
CX, CY, HALF = 139.751, 35.523, 0.026

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])

m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
T = m.elements
land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")

fig, axes = plt.subplots(1, 2, figsize=(15, 8))
for ax, (ttl, lo, la, tt, col) in zip(axes, [
        ("goto2023 sample", lon_s, lat_s, Ts, "0.35"),
        ("ours (anchor-rule fix)", lon_o, lat_o, T,
         "steelblue")]):
    gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5,
              zorder=1)
    ax.triplot(lo, la, tt, lw=0.6, color=col, zorder=3)
    ax.set_xlim(CX - HALF / COSW, CX + HALF / COSW)
    ax.set_ylim(CY - HALF, CY + HALF)
    ax.set_aspect(1 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(ttl)
fig.suptitle("Shinagawa-Kawasaki canal band (G8): sample omits, width rule keeps")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("outputs/figures/g8_band_zoom.png", dpi=180,
            bbox_inches="tight")
print("[g8] saved", flush=True)

raise SystemExit(0)
