import time, numpy as np
P = "/home/pj24001722/ku40000343/Github/OceanMesh2D/datasets/PostSandyNCEI/PostSandyNCEI.nc"
bbox = (-73.97, -73.75, 40.5, 40.68)
t0=time.time()
import rasterio
from rasterio.windows import from_bounds
with rasterio.open(P) as src:
    win = from_bounds(bbox[0], bbox[2], bbox[1], bbox[3],
                      transform=src.transform)
    a = src.read(1, window=win, masked=True)
print(f"[bench] rasterio windowed: {time.time()-t0:.1f}s shape={a.shape}", flush=True)
t0=time.time()
import netCDF4
ds = netCDF4.Dataset(P)
names = list(ds.variables)
lat = None; lon = None; zv = None
for n in names:
    v = ds.variables[n]
    if n.lower() in ("lat","latitude","y"): lat=v
    elif n.lower() in ("lon","longitude","x"): lon=v
    elif v.ndim==2: zv=v
la = lat[:]; lo = lon[:]
i0,i1 = np.searchsorted(lo, [bbox[0], bbox[1]])
la_asc = la[1] > la[0]
if la_asc: j0,j1 = np.searchsorted(la, [bbox[2], bbox[3]])
else:
    j0 = len(la)-np.searchsorted(la[::-1], bbox[3]); j1 = len(la)-np.searchsorted(la[::-1], bbox[2])
z = zv[j0:j1, i0:i1]
print(f"[bench] netCDF4 sliced: {time.time()-t0:.1f}s shape={z.shape}", flush=True)
t0=time.time()
np.save("/tmp/demcache_test.npy", np.asarray(z))
z2 = np.load("/tmp/demcache_test.npy")
print(f"[bench] npy save+load: {time.time()-t0:.1f}s", flush=True)
