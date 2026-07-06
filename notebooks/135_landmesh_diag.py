import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast

m = read_fort14('outputs/pipeline_v6r/tokyo_bay_v6_final.14')
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
land_open = gpd.read_file('outputs/pipeline_v6r/prep/land_opened.shp')
land_solid = gpd.read_file('outputs/pipeline_v6r/prep/land_solid.shp')

x0, y0, x1, y1 = 139.78, 35.32, 139.95, 35.44   # I11/J11 + margin
fig, ax = plt.subplots(figsize=(12, 11))
_add_coast(ax, (139.0, 34.5, 141.3, 36.2), "EPSG:4326")
ax.triplot(lon, lat, m.elements, lw=0.4, color="steelblue")
land_open.boundary.plot(ax=ax, color="darkgreen", lw=1.2,
                        label="land_opened (engine input)")
land_solid.boundary.plot(ax=ax, color="purple", lw=0.8,
                         linestyle="--", label="land_solid (pre-prep)")
for gx in np.arange(139.80, 139.95, 0.05):
    ax.axvline(gx, color="crimson", lw=0.6, alpha=0.6)
for gy in np.arange(35.35, 35.45, 0.05):
    ax.axhline(gy, color="crimson", lw=0.6, alpha=0.6)
ax.text(139.825, 35.428, "I", color="crimson", fontsize=14, fontweight="bold")
ax.text(139.875, 35.428, "J", color="crimson", fontsize=14, fontweight="bold")
ax.text(139.788, 35.372, "11", color="crimson", fontsize=14, fontweight="bold")
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
ax.set_aspect(1 / np.cos(np.deg2rad(35.38)))
ax.legend(loc="lower right")
ax.set_title("I11/J11 diagnosis: mesh vs true OSM land (gray) vs "
             "land_opened input (green)")
fig.savefig('outputs/figures/diag_I11J11.png', dpi=230,
            bbox_inches="tight")
print("saved", flush=True)

# quantify: element centroids inside TRUE OSM land here
import xcoast
mask = xcoast.load((139.78, 35.32, 139.95, 35.44))
gdf = None
for attr in ("gdf", "land", "polygons"):
    if hasattr(mask, attr):
        gdf = getattr(mask, attr)
        break
print("xcoast object:", type(mask), "->", type(gdf), flush=True)
