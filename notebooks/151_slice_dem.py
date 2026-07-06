# Recreate the .m's "SRTM15_V2.3_sliced" intent: a coarse basin DEM
# for nest 1 (10-50 km sizing needs no 15-arcsec data) and a
# full-res Kanto slice for nests 2-3.
import xarray as xr
from pathlib import Path
SRC = Path.home() / "Github/OceanMesh2D/datasets/SRTM15+.nc"
OUT = Path.home() / "Github/OceanMesh2D/datasets/TokyoBay/dem"
OUT.mkdir(parents=True, exist_ok=True)
ds = xr.open_dataset(SRC)
lon = [c for c in ds.coords if c.lower() in ("lon", "longitude", "x")][0]
lat = [c for c in ds.coords if c.lower() in ("lat", "latitude", "y")][0]
big = ds.sel({lon: slice(114, 168), lat: slice(17, 76)})
big16 = big.coarsen({lon: 16, lat: 16}, boundary="trim").mean()
big16.to_netcdf(OUT / "SRTM15_pacific_4min.nc")
print("pacific 4min:", {d: big16.sizes[d] for d in big16.sizes}, flush=True)
kanto = ds.sel({lon: slice(137.8, 141.6), lat: slice(33.2, 36.3)})
kanto.to_netcdf(OUT / "SRTM15_kanto_15s.nc")
print("kanto 15s:", {d: kanto.sizes[d] for d in kanto.sizes}, flush=True)
