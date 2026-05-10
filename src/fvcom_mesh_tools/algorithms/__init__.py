"""Mesh-modifying algorithms (perpendicularity, smoothing, ...)."""

from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    classify_boundaries_by_bbox,
    classify_outer_loop_by_bbox,
    outer_loop,
)
from fvcom_mesh_tools.algorithms.edge_swap import (
    swap_edges_for_quality,
    swap_edges_for_valence,
)
from fvcom_mesh_tools.algorithms.perpendicularity import (
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
    unique_edges,
)
from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    edge_lengths_planar,
    min_interior_angle,
)
from fvcom_mesh_tools.algorithms.refine import refine_bad_triangles
from fvcom_mesh_tools.algorithms.rivers import add_river_inflow_segments
from fvcom_mesh_tools.algorithms.smoothing import laplacian_smooth

__all__ = [
    "add_river_inflow_segments",
    "align_open_boundary_first_ring",
    "alpha_quality",
    "boundary_edges_from_tris",
    "chain_edges_to_loops",
    "classify_boundaries_by_bbox",
    "classify_outer_loop_by_bbox",
    "edge_lengths_planar",
    "laplacian_smooth",
    "min_interior_angle",
    "open_bdy_perpendicularity",
    "outer_loop",
    "refine_bad_triangles",
    "signed_areas",
    "swap_edges_for_quality",
    "swap_edges_for_valence",
    "unique_edges",
]
