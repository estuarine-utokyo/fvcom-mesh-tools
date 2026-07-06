# OSM true-land for the varres nest-2 window (the .m's bbox_02
# pentagon extent) — the project's engineered OSM coastline
# replaces coastline_2 in the translation (user policy).
from pathlib import Path
from fvcom_mesh_tools.prep import fetch_true_land

OUT = Path("outputs/tb_varres_3r")
OUT.mkdir(parents=True, exist_ok=True)
land = fetch_true_land((138.19, 33.52, 141.30, 35.94),
                       min_water_area_deg2=5e-5)
land.to_file(OUT / "land_osm_wide.shp")
print(f"[osm] wide true-land: {len(land)} polygons", flush=True)
