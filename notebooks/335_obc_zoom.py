# OBC ladder showcase + the single residual C4 site
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
gd = os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')
g2 = open(gd).read().split('\n')
nn = int(g2[0].split('=')[1]); ne = int(g2[1].split('=')[1])
Ts = np.array([[int(w) for w in g2[2+i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in g2[2+ne+i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
obc_s = [int(l.split()[1]) - 1 for l in open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_obc.dat')
    ).read().strip().split('\n')[1:]]

m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
ob = np.asarray(m.open_boundaries[0], int)

COAST = (139.0, 34.5, 141.3, 36.2)
x0, y0, x1, y1 = 139.63, 34.955, 139.82, 35.16
fig, axes = plt.subplots(1, 2, figsize=(19, 11))
for ax, (tag, lo, la, tt, obn) in zip(axes, [
        ("goto2023 SAMPLE", lon_s, lat_s, Ts, obc_s),
        ("repro FINAL (ladder OBC)", lon, lat, m.elements, ob)]):
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lo, la, tt, lw=0.8, color="steelblue")
    ax.plot(np.asarray(lo)[obn], np.asarray(la)[obn], color="red",
            lw=2.5, zorder=6, marker="o", ms=4)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_aspect(1 / np.cos(np.deg2rad(35.05)))
    ax.set_title(tag)
fig.savefig("outputs/figures/obc_ladder_zoom.png", dpi=200,
            bbox_inches="tight")
print("saved obc zoom", flush=True)
