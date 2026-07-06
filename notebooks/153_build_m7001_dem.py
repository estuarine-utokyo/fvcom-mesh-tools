# WIDE M7001 DEM: grid the full J-BIRD N-mark point set (3.95M pts,
# Southern Kanto) onto the kanto 15-arcsec lattice; SRTM fills
# outside the M7001 hull and on land (SRTM z > +5 m guard, since
# contour-vertex interpolation spans land between shores).
import os
import numpy as np
import xarray as xr
from pathlib import Path
from scipy.interpolate import LinearNDInterpolator

ASCII = os.path.expandvars(
    "$DATA_DIR/bathymetry/M7001/ascii/M7001_関東南部_Ver.2.4")
OUT = Path.home() / "Github/OceanMesh2D/datasets/TokyoBay/dem"
srtm = xr.open_dataset(OUT / "SRTM15_kanto_15s.nc")
var = [v for v in srtm.data_vars if srtm[v].ndim == 2][0]

lo, la, de = [], [], []
with open(ASCII, "r", encoding="cp932", errors="replace") as fh:
    for line in fh:
        if len(line) < 36 or line[0] != "N":
            continue
        try:
            d = float(line[10:17])
            lat = float(line[17:26])
            lon = float(line[26:36])
        except ValueError:
            continue
        lo.append(lon); la.append(lat); de.append(d)
lo = np.array(lo); la = np.array(la); de = np.array(de)
print(f"[m7001] N-mark points: {len(de):,} | lon "
      f"{lo.min():.3f}-{lo.max():.3f} lat {la.min():.3f}-{la.max():.3f} "
      f"| depth {de.min():.0f}..{de.max():.0f} m", flush=True)

itp = LinearNDInterpolator(np.column_stack([lo, la]), -de)
LON, LAT = np.meshgrid(srtm.lon.values, srtm.lat.values)
zi = itp(np.column_stack([LON.ravel(), LAT.ravel()])).reshape(LON.shape)
zs = srtm[var].values
use = np.isfinite(zi) & (zs < 5.0)
z = np.where(use, zi, zs)
print(f"[m7001] lattice cells from M7001: {int(use.sum()):,} / {z.size:,}",
      flush=True)
ds = xr.Dataset({var: (("lat", "lon"), z)},
                coords={"lat": srtm.lat, "lon": srtm.lon})
ds.to_netcdf(OUT / "kanto_M7001full_srtm_15s.nc")
print("[m7001] wrote kanto_M7001full_srtm_15s.nc", flush=True)
