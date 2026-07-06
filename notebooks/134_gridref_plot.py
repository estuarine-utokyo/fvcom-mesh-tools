import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID as G
from fvcom_mesh_tools.plotting import _add_coast

m = read_fort14('outputs/pipeline_v6r/tokyo_bay_v6_final.14')
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
COAST = (139.0, 34.5, 141.3, 36.2)  # true OSM land via xcoast —
# land_opened.shp contains the artificial OBC wall band and must
# not be used for display

LANDMARKS = {
    "Tokyo": (139.79, 35.635), "Kawasaki": (139.76, 35.51),
    "Yokohama": (139.66, 35.44), "Yokosuka": (139.67, 35.28),
    "Kannonzaki": (139.74, 35.255), "Futtsu": (139.79, 35.31),
    "Banzu": (139.92, 35.40), "Chiba": (140.03, 35.58),
    "Uraga OBC": (139.71, 35.06),
}

def draw(fname, x0, y0, x1, y1, lw, dpi, title):
    fig, ax = plt.subplots(figsize=(12, 12 * (y1 - y0) * 1.22 /
                                    (x1 - x0)))
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon, lat, m.elements, lw=lw, color="steelblue")
    for seg in m.open_boundaries:
        seg = [int(v) for v in seg]
        ax.plot(lon[seg], lat[seg], color="red", lw=2.4, zorder=6,
                label="open boundary")
    if m.open_boundaries:
        ax.legend(loc="lower right", fontsize=9)
    # grid lines + labels
    for i in range(G.ncol + 1):
        gx = G.lon0 + i * G.dlon
        if x0 - G.dlon < gx < x1 + G.dlon:
            ax.axvline(gx, color="crimson", lw=0.5, alpha=0.55)
    for j in range(G.nrow + 1):
        gy = G.lat1 - j * G.dlat
        if y0 - G.dlat < gy < y1 + G.dlat:
            ax.axhline(gy, color="crimson", lw=0.5, alpha=0.55)
    for i in range(G.ncol):
        gx = G.lon0 + (i + 0.5) * G.dlon
        if x0 < gx < x1:
            for yy in (y0 + 0.006, y1 - 0.012):
                ax.text(gx, yy, G.col_letter(i), color="crimson",
                        ha="center", fontsize=11, fontweight="bold")
    for j in range(1, G.nrow + 1):
        gy = G.lat1 - (j - 0.5) * G.dlat
        if y0 < gy < y1:
            for xx in (x0 + 0.004, x1 - 0.012):
                ax.text(xx, gy, str(j), color="crimson",
                        va="center", fontsize=11, fontweight="bold")
    for name, (px, py) in LANDMARKS.items():
        if x0 < px < x1 and y0 < py < y1:
            ax.annotate(name, (px, py), fontsize=9, color="black",
                        fontweight="bold",
                        bbox=dict(fc="lemonchiffon", alpha=0.75,
                                  ec="none", pad=1.2))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (y0 + y1))))
    ax.set_title(title)
    fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("saved", fname, flush=True)

draw('outputs/figures/v6_gridref_full.png',
     139.40, 34.90, 140.30, 35.90, 0.25, 200,
     "tokyo_bay_v6 - grid reference sheet (cols A-R west->east, "
     "rows 1-20 north->south, cell 0.05deg ~ 4.5x5.5 km)")
draw('outputs/figures/v6_gridref_bay.png',
     139.55, 35.15, 140.20, 35.75, 0.35, 230,
     "tokyo_bay_v6 - bay interior grid reference "
     "(quadrants: a=NW b=NE c=SW d=SE, e.g. F12c)")
