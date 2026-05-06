"""I/O routines for unstructured mesh formats used by FVCOM/ADCIRC.

This subpackage stays free of compiled raster dependencies (rasterio /
netCDF4 / GDAL); modules under :mod:`fvcom_mesh_tools.dem` cover that
ground. shapely / geopandas / pyproj are imported lazily by
``coastline``, ``geom_filter``, and ``rivers`` so a minimal
``import fvcom_mesh_tools.io`` brings in only stdlib + numpy.
"""

from fvcom_mesh_tools.io.coastline import load_coastline_as_lines
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.geom_filter import filter_multipolygon_by_area
from fvcom_mesh_tools.io.rivers import load_river_points

__all__ = [
    "Fort14Mesh",
    "filter_multipolygon_by_area",
    "load_coastline_as_lines",
    "load_river_points",
    "read_fort14",
    "write_fort14",
]
