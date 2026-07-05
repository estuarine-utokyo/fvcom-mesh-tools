# PoC #115: does lowering min_water_area_deg2 open the Arakawa and
# Tamagawa river channels in the xcoast true-land product?
from shapely import unary_union
from shapely.geometry import Point

from fvcom_mesh_tools.prep import fetch_true_land

BBOX = (139.50, 35.40, 140.05, 35.75)  # river-mouth region only
PROBES = {
    "arakawa_mouth": (139.8565, 35.6303),
    "arakawa_up2km": (139.8600, 35.6500),
    "tamagawa_mouth": (139.7605, 35.5303),
    "tamagawa_up2km": (139.7420, 35.5480),
    "sumidagawa": (139.7967, 35.6236),
    "edogawa_mouth": (139.9050, 35.6480),
}
for thr in (1e-5, 1e-6, 1e-7):
    gdf = fetch_true_land(BBOX, min_water_area_deg2=thr, force=True)
    u = unary_union(list(gdf.geometry))
    states = {k: ("LAND" if u.contains(Point(*p)) else "water")
              for k, p in PROBES.items()}
    print(f"thr={thr:g}: " + "  ".join(f"{k}={v}" for k, v in states.items()),
          flush=True)
