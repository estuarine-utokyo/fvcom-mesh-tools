"""Mesh-modifying algorithms (perpendicularity, smoothing, ...)."""

from fvcom_mesh_tools.algorithms.perpendicularity import (
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
    unique_edges,
)

__all__ = [
    "align_open_boundary_first_ring",
    "open_bdy_perpendicularity",
    "signed_areas",
    "unique_edges",
]
