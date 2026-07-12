# Zoom the C1/C4 quality-tail site (element 3894 area, Chiba kept
# channel) with mesh + pipeline land + original land.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
from pyproj import Transformer
from shapely.ops import unary_union
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
CX, CY = tr.transform(373450.0, 3940300.0)
m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
T = m.elements
land0 = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
landp = unary_union(list(gpd.read_file(
    "outputs/sample_repro/land_channel_adj.shp").geometry))
COSW = float(np.cos(np.deg2rad(35.35)))
fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
for ax, (ttl, ld) in zip(axes, [("original OSM land", land0),
                                ("pipeline land (carved)", landp)]):
    gpd.GeoSeries([ld], crs="EPSG:4326").plot(
        ax=ax, color="0.9", edgecolor="0.5", lw=0.6, zorder=1)
    ax.triplot(lon, lat, T, lw=0.6, color="steelblue", zorder=3)
    for e in (5126 - 1, 4704 - 1, 2666 - 1, 4845 - 1):
        ring = np.column_stack([lon, lat])[np.append(T[e], T[e][0])]
        ax.plot(ring[:, 0], ring[:, 1], color="crimson", lw=1.6,
                zorder=4)
        c = ring[:3].mean(axis=0)
        ax.annotate(str(e + 1), c, ha="center", fontsize=11,
                    color="darkred", zorder=5)
    half = 0.014
    ax.set_xlim(CX - half / COSW, CX + half / COSW)
    ax.set_ylim(CY - half, CY + half)
    ax.set_aspect(1 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(ttl)
fig.suptitle(f"C1/C4 tail site at ({CX:.4f}, {CY:.4f})")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("outputs/figures/c4_tail_site.png", dpi=180,
            bbox_inches="tight")
print("[probe] saved outputs/figures/c4_tail_site.png", flush=True)
