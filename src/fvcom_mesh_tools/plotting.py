"""Shared matplotlib defaults for fvcom-mesh-tools."""

from __future__ import annotations

MESH_PNG_DPI: int = 600
"""Default raster DPI for ``*.png`` outputs that visualise a mesh.

Mesh triangulations stay readable when zoomed in only at high DPI; 600 dpi
is the project-wide default for any figure showing the triangulation
itself (``triplot`` / ``tripcolor`` plots, boundary maps, side-by-side
mesh comparisons). Histograms and other non-spatial plots can keep the
matplotlib default. Pass this value explicitly to
``fig.savefig(..., dpi=MESH_PNG_DPI)`` so that the choice is visible at
the call site and consistent across notebooks and scripts.
"""

__all__ = ["MESH_PNG_DPI"]
