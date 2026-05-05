"""Mesh-modifying algorithms (perpendicularity, smoothing, ...)."""

from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    classify_boundaries_by_bbox,
    classify_outer_loop_by_bbox,
    outer_loop,
)
from fvcom_mesh_tools.algorithms.perpendicularity import (
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
    unique_edges,
)

__all__ = [
    "align_open_boundary_first_ring",
    "boundary_edges_from_tris",
    "chain_edges_to_loops",
    "classify_boundaries_by_bbox",
    "classify_outer_loop_by_bbox",
    "open_bdy_perpendicularity",
    "outer_loop",
    "signed_areas",
    "unique_edges",
]
