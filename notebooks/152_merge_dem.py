# Kanto DEM for the varres translation: M7001-derived depth_0090
# burned onto the SRTM15 kanto lattice (project bathymetry policy:
# M7001 wherever it covers; SRTM elsewhere).
import os
import numpy as np
import xarray as xr
from pathlib import Path

OUT = Path.home() / "Github/OceanMesh2D/datasets/TokyoBay/dem"
srtm = xr.open_dataset(OUT / "SRTM15_kanto_15s.nc")
m = xr.open_dataset(os.path.expandvars(
    "$DATA_DIR/geodata/bathymetry/tokyo_bay/depth_0090-05+06+07.nc"))
var_s = [v for v in srtm.data_vars if srtm[v].ndim == 2][0]
print("srtm var:", var_s, "| m7001 Band1 range:",
      float(m.Band1.min()), float(m.Band1.max()), flush=True)
mi = m.Band1.interp(lon=srtm.lon, lat=srtm.lat, method="linear")
z = xr.where(np.isfinite(mi), mi, srtm[var_s])
ds = z.to_dataset(name=var_s)
ds.to_netcdf(OUT / "kanto_M7001_srtm_15s.nc")
n_m = int(np.isfinite(mi).sum())
print(f"merged: {n_m} cells from M7001, "
      f"{int(z.size) - n_m} from SRTM -> kanto_M7001_srtm_15s.nc",
      flush=True)
