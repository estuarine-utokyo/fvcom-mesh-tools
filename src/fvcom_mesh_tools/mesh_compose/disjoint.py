"""``combine_disjoint``: concatenate non-overlapping meshes into one.

Produces a :class:`Fort14Mesh` whose nodes / elements / depths are the
straight concatenation of the inputs, with element and boundary
indices renumbered so they continue to point at the right vertices.
Boundary segments are carried forward verbatim - we never invent or
delete a boundary, which is the right thing for the disjoint case
(two basins on opposite ends of the country don't share a coastline).

This strategy is pure numpy; OCSMesh is not required.
"""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def combine_disjoint(
    meshes: list[Fort14Mesh],
    *,
    title: str | None = None,
) -> Fort14Mesh:
    """Concatenate ``meshes`` with index offsets; preserve all boundaries."""
    if len(meshes) < 2:
        raise ValueError("combine_disjoint requires at least 2 meshes")

    # Stack nodes / depths / elements with vertex-index offsets.
    nodes = np.vstack([m.nodes for m in meshes]).astype(float)
    depths = np.concatenate([m.depths for m in meshes]).astype(float)
    elements_chunks: list[np.ndarray] = []
    open_bnds: list[np.ndarray] = []
    land_bnds: list[tuple[int, np.ndarray]] = []

    offset = 0
    for m in meshes:
        elements_chunks.append(np.asarray(m.elements, dtype=np.int64) + offset)
        for seg in m.open_boundaries:
            open_bnds.append(np.asarray(seg, dtype=np.int64) + offset)
        for ibtype, seg in m.land_boundaries:
            land_bnds.append((int(ibtype), np.asarray(seg, dtype=np.int64) + offset))
        offset += m.n_nodes

    elements = np.vstack(elements_chunks).astype(np.int64)

    if title is None:
        title = " + ".join(m.title or "(untitled)" for m in meshes)

    return Fort14Mesh(
        title=title,
        nodes=nodes,
        depths=depths,
        elements=elements,
        open_boundaries=open_bnds,
        land_boundaries=land_bnds,
    )
