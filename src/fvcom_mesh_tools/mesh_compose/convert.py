"""Conversion helpers between :class:`Fort14Mesh` and OCSMesh ``MeshData``.

OCSMesh's combine helpers (``ocsmesh.ops.combine_mesh``) operate on
``MeshData`` objects defined in ``ocsmesh.internal``. We bridge to and
from our :class:`Fort14Mesh` here so the rest of ``mesh_compose`` does
not need to care about either side's internals. Boundaries do *not*
round-trip - ``MeshData`` has no boundary structure, so combine
strategies that go through OCSMesh return a mesh with empty boundary
lists (the caller is expected to re-classify or re-attach boundaries
afterwards).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def fort14_to_meshdata(mesh: Fort14Mesh, *, crs: str = "EPSG:4326") -> Any:
    """Wrap a :class:`Fort14Mesh` as an OCSMesh ``MeshData``.

    Depths are passed through as nodal ``values``. Boundaries are
    dropped (``MeshData`` does not represent them).
    """
    from ocsmesh.internal import MeshData  # lazy import

    return MeshData(
        coords=np.asarray(mesh.nodes, dtype=float),
        tria=np.asarray(mesh.elements, dtype=int),
        values=np.asarray(mesh.depths, dtype=float),
        crs=crs,
    )


def meshdata_to_fort14(
    md: Any,
    *,
    title: str = "fmesh-mesh-combine output",
) -> Fort14Mesh:
    """Wrap an OCSMesh ``MeshData`` back into a :class:`Fort14Mesh`.

    Boundary lists are empty - call ``classify_boundaries_by_bbox`` (or
    pass through ``omesh14-edit-bdy``) downstream to re-establish them.
    """
    nodes = np.asarray(md.coords, dtype=float)
    elements = np.asarray(md.tria, dtype=np.int64)
    raw_values = np.asarray(md.values).reshape(-1)
    if raw_values.shape[0] != nodes.shape[0]:
        depths = np.zeros(nodes.shape[0], dtype=float)
    else:
        depths = raw_values.astype(float)
    return Fort14Mesh(
        title=title,
        nodes=nodes,
        depths=depths,
        elements=elements,
        open_boundaries=[],
        land_boundaries=[],
    )
