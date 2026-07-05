"""Shoreline preprocessing for mesh generation (the v5 front end).

Promotes the PoC #90/#91 one-off scripts into parameterized,
machine-independent modules:

* :mod:`.shoreline` — OSM true-land fetch (via xcoast) and the
  LAND-opening simplification (erase piers/breakwaters/islets
  thinner than the mesh scale without ever cutting water
  connectivity);
* :mod:`.skeleton` — water medial-axis seed lines for
  narrow-but-essential water (harbor basins, river mouths) that the
  DistMesh seed lattice would otherwise drop.
"""

from fvcom_mesh_tools.prep.shoreline import (
    auto_utm_epsg,
    fetch_true_land,
    open_land,
)
from fvcom_mesh_tools.prep.skeleton import water_skeleton_lines

__all__ = [
    "auto_utm_epsg",
    "fetch_true_land",
    "open_land",
    "water_skeleton_lines",
]
