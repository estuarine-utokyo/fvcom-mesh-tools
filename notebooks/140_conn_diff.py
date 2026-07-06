import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
from fvcom_mesh_tools.prep import fetch_true_land

land = fetch_true_land((139.40, 34.90, 140.30, 35.90),
                       min_water_area_deg2=1e-6)
views = {"tamagawa": (139.70, 35.48, 139.80, 35.58),
         "arakawa": (139.82, 35.60, 139.90, 35.68),
         "obitsu": (139.78, 35.32, 139.95, 35.44)}
import geopandas as gpd
old = gpd.read_file('outputs/pipeline_v6r/prep/land_solid.shp')
for name, (x0, y0, x1, y1) in views.items():
    fig, ax = plt.subplots(figsize=(9, 8))
    land.plot(ax=ax, color="0.8", edgecolor="0.4", lw=0.3)
    old.boundary.plot(ax=ax, color="orange", lw=0.8, linestyle=":")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (y0 + y1))))
    ax.set_title(f"sea-connected subtraction: {name} "
                 "(gray=new land, orange=previous fclass land)")
    fig.savefig(f'outputs/figures/conn_{name}.png', dpi=200,
                bbox_inches="tight")
    print("saved", name, flush=True)
