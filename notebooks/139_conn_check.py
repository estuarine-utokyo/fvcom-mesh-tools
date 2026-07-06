# quick semantic check of sea-connectivity subtraction:
# Tamagawa mouth must be WATER, Obitsu hinterland must be LAND,
# an inland lake must be LAND.
import shapely
from shapely.ops import unary_union
from fvcom_mesh_tools.prep import fetch_true_land

land = fetch_true_land((139.40, 34.90, 140.30, 35.90),
                       min_water_area_deg2=1e-6)
u = unary_union(list(land.geometry))
probes = {
    "Tamagawa mouth (want WATER)": (139.7355, 35.5320, False),
    "Arakawa mouth (want WATER)": (139.8570, 35.6430, False),
    "Obitsu hinterland (want LAND)": (139.8600, 35.3700, True),
    "Lake Inba-ish inland (want LAND)": (140.2000, 35.7500, True),
}
ok = True
for name, (lon, lat, want_land) in probes.items():
    is_land = bool(shapely.contains_xy(u, lon, lat))
    mark = "OK" if is_land == want_land else "NG"
    ok &= (is_land == want_land)
    print(f"[check] {mark} {name}: land={is_land}", flush=True)
print("[check] ALL OK" if ok else "[check] FAILURES", flush=True)
