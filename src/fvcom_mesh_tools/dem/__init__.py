"""DEM I/O helpers (rasterio-bound).

This subpackage isolates every module that needs ``rasterio`` /
``netCDF4`` so that pure-fort.14 callers (``algorithms``, ``io``,
``mesh_compose``, ``cli/perpfix``, ``cli/meshcombine``) do not pull
the compiled raster stack in as a transitive dependency.

Modules:

* :mod:`fvcom_mesh_tools.dem.subset` -- clip a global DEM to a bbox
  and emit a CF-tagged GeoTIFF (``fmesh-subset-dem`` backend).
* :mod:`fvcom_mesh_tools.dem.interp` -- bilinear / nearest sampling
  of a DEM at mesh-node coordinates (post-mesh depth interpolation).
* :mod:`fvcom_mesh_tools.dem.bbox` -- read the lon/lat extent of a
  raster, used by ``fmesh-buildmesh`` for boundary classification.

The runtime dependency is declared via the ``[dem]`` extra in
``pyproject.toml``. Importing any submodule without rasterio /
netCDF4 installed raises ``ImportError``.
"""
