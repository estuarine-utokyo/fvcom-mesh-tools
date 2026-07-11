import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
gd = os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')
g2 = open(gd).read().split('\n')
nn = int(g2[0].split('=')[1]); ne = int(g2[1].split('=')[1])
Ts = np.array([[int(w) for w in g2[2+i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in g2[2+ne+i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
COAST = (139.0, 34.5, 141.3, 36.2)
fig, axes = plt.subplots(1, 2, figsize=(18, 11))
for ax, (tag, lo, la, tt) in zip(axes, [
        ("goto2023 SAMPLE", lon_s, lat_s, Ts),
        ("repro FINAL", lon, lat, m.elements)]):
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lo, la, tt, lw=0.8, color="steelblue")
    ax.set_xlim(140.02, 140.14)
    ax.set_ylim(35.53, 35.68)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_aspect(1 / np.cos(np.deg2rad(35.6)))
    ax.set_title(f"{tag} - Hanami/Kemigawa area")
fig.savefig("outputs/figures/hanami_zoom.png", dpi=190,
            bbox_inches="tight")
print("saved", flush=True)
