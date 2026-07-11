import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

m = read_fort14('outputs/tb_varres_3r/tb_varres_3regions.14')
lon, lat = m.nodes[:, 0], m.nodes[:, 1]
COAST = (139.0, 34.5, 141.3, 36.2)
views = {"wide": (136.5, 32.5, 142.5, 37.0, 0.15),
         "bay": (139.55, 35.15, 140.20, 35.75, 0.35),
         "mouth": (139.60, 34.95, 140.05, 35.30, 0.45)}
for name, (x0, y0, x1, y1, lw) in views.items():
    fig, ax = plt.subplots(figsize=(11, 11*(y1-y0)*1.2/(x1-x0)))
    if name != "wide":
        _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon, lat, m.elements, lw=lw, color="steelblue")
    for seg in m.open_boundaries:
        seg = [int(v) for v in seg]
        ax.plot(lon[seg], lat[seg], color="red", lw=1.6, zorder=6)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1/np.cos(np.deg2rad(0.5*(y0+y1))))
    if name != "wide":
        add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(f"tb_varres_3regions (1:1 translation) "
                 f"NP={len(m.nodes):,} - {name}")
    fig.savefig(f'outputs/figures/trans_{name}.png', dpi=200,
                bbox_inches="tight")
    print("saved", name, flush=True)
