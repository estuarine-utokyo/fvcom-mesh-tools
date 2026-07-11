# Side-by-side: oceanmesh output (BEFORE) vs FVCOM-finished (AFTER)
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

OUT = "outputs/figures"
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
panels = []
for tag, path in [("oceanmesh output (BEFORE)",
                   "outputs/sample_repro/sample_repro_utm.14"),
                  ("FVCOM finished (AFTER)",
                   "outputs/sample_repro/sample_repro_final.14")]:
    m = read_fort14(path)
    lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
    panels.append((tag, lon, lat, m))

COAST = (139.0, 34.5, 141.3, 36.2)
views = {
    "full": (139.57, 34.93, 140.15, 35.78, 0.35),
    "mouth": (139.60, 34.95, 140.05, 35.30, 0.55),
    "tama": (139.73, 35.48, 139.87, 35.58, 0.7),
}
for name, (x0, y0, x1, y1, lw) in views.items():
    fig, axes = plt.subplots(
        1, 2, figsize=(19, 9.5 * (y1 - y0) * 1.22 / (x1 - x0)))
    for ax, (tag, lon, lat, m) in zip(axes, panels):
        _add_coast(ax, COAST, "EPSG:4326")
        ax.triplot(lon, lat, m.elements, lw=lw, color="steelblue")
        ob = np.asarray(m.open_boundaries[0], int)
        ax.plot(lon[ob], lat[ob], color="red", lw=2.5, zorder=6,
                marker="o", ms=3, label=f"OBC ({len(ob)} nodes)")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (y0 + y1))))
        ax.legend(loc="lower right")
        ax.set_title(f"{tag}  NP={m.n_nodes:,} NE={len(m.elements):,}")
    fig.savefig(f"{OUT}/ba_cmp_{name}.png", dpi=200,
                bbox_inches="tight")
    print("saved", name, flush=True)
