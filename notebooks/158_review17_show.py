import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

m = read_fort14('outputs/pipeline_v5u/tokyo_bay_v5u_final.14')
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
COAST = (139.0, 34.5, 141.3, 36.2)
views = {"full": (139.40, 34.90, 140.30, 35.90, 0.3),
         "bay": (139.55, 35.15, 140.20, 35.75, 0.4),
         "mouth": (139.60, 34.95, 140.05, 35.30, 0.5),
         "urayasu": (139.82, 35.58, 140.00, 35.70, 0.6)}
for name, (x0, y0, x1, y1, lw) in views.items():
    fig, ax = plt.subplots(figsize=(11, 11*(y1-y0)*1.2/(x1-x0)))
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon, lat, m.elements, lw=lw, color="steelblue")
    for seg in m.open_boundaries:
        seg = [int(v) for v in seg]
        ax.plot(lon[seg], lat[seg], color="red", lw=2.0, zorder=6)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1/np.cos(np.deg2rad(0.5*(y0+y1))))
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(f"review17-class mesh (v5u final, CDT boundary) "
                 f"NP={len(m.nodes):,} - {name}")
    fig.savefig(f'outputs/figures/r17_{name}.png', dpi=210,
                bbox_inches="tight")
    print("saved", name, flush=True)
