# The 5 REAL w/h<1 throat sites, large zooms, sample side-by-side.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.diagnostics import channel_width_metric
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
P = np.column_stack([lon, lat])
T = m.elements
gd = os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')
g2 = open(gd).read().split('\n')
nn = int(g2[0].split('=')[1]); ne = int(g2[1].split('=')[1])
Ts = np.array([[int(w) for w in g2[2+i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in g2[2+ne+i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])

ch = channel_width_metric(m, coords="metric")
r = ch["w_h_ratio"]
bad = set(np.where(np.isfinite(r) & (r < 1.0))[0].tolist())

SITES = [
    ("Keihin canal (Tsurumi) - through, w/h 0.54 x7",
     139.681, 35.470, 0.018),
    ("Chiba inlet (Kemigawa) - dead-end, x3",
     140.113, 35.592, 0.016),
    ("Ichihara breakwater channel - through, x2",
     139.811, 35.537, 0.018),
    ("Chiba port inner channel - dead-end, x2",
     139.982, 35.670, 0.016),
    ("Kawasaki canal mouth - x1",
     139.775, 35.507, 0.016),
]
COAST = (139.0, 34.5, 141.3, 36.2)
fig, axes = plt.subplots(len(SITES), 2, figsize=(17, 7.2 * len(SITES)))
for row, (name, cx, cy, R) in enumerate(SITES):
    for col, (tag, lo, la, tt) in enumerate([
            ("goto2023 SAMPLE", lon_s, lat_s, Ts),
            ("repro FINAL", lon, lat, T)]):
        ax = axes[row, col]
        _add_coast(ax, COAST, "EPSG:4326")
        ax.triplot(lo, la, tt, lw=0.9, color="steelblue")
        if col == 1:
            cent = P[T].mean(axis=1)
            近 = [ei for ei in bad
                  if abs(cent[ei, 0] - cx) < R
                  and abs(cent[ei, 1] - cy) < R]
            for ei in 近:
                ax.fill(P[T[ei], 0], P[T[ei], 1], color="crimson",
                        alpha=0.55, zorder=5)
        ax.set_xlim(cx - R, cx + R)
        ax.set_ylim(cy - R * 0.8, cy + R * 0.8)
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_aspect(1 / np.cos(np.deg2rad(cy)))
        ax.set_title(f"{tag}\n{name}", fontsize=10)
fig.savefig("outputs/figures/throat_zooms.png", dpi=180,
            bbox_inches="tight")
print("saved throat zooms", flush=True)
