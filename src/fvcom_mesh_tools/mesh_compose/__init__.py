"""Multi-mesh composition: stitch independent ``fort.14`` files together.

The ``oceanmesh`` engine excels at single-mesh generation but does not
provide post-hoc mesh-combination utilities. ``ocsmesh.ops.combine_mesh``
does, with two flavours:

* ``merge_overlapping_meshes`` carves a high-resolution foreground
  mesh into a coarser background mesh, blending the seam.
* ``merge_neighboring_meshes`` snaps two meshes whose boundaries
  already coincide (within tolerance), with no carving.

This module wraps both and adds a third trivial-but-useful strategy:

* ``disjoint`` simply renumbers and concatenates two non-overlapping
  meshes (e.g. Tokyo Bay + Osaka Bay) into one ``fort.14`` with
  preserved boundary segments. Pure numpy; no OCSMesh required.

Use the ``combine`` dispatcher for a uniform CLI/API:

    from fvcom_mesh_tools.mesh_compose import combine
    out_mesh = combine("disjoint", [tokyo_mesh, osaka_mesh])
"""

from __future__ import annotations

from collections.abc import Iterable

from fvcom_mesh_tools.io import Fort14Mesh

STRATEGIES = ("disjoint", "overlap", "neighbor")


def combine(
    strategy: str,
    meshes: Iterable[Fort14Mesh],
    **kwargs,
) -> Fort14Mesh:
    """Combine ``meshes`` using the named ``strategy``.

    Parameters
    ----------
    strategy
        One of ``"disjoint"``, ``"overlap"``, ``"neighbor"``.
    meshes
        Two or more :class:`Fort14Mesh` instances. ``disjoint`` and
        ``neighbor`` are commutative; ``overlap`` is order-sensitive
        (the first mesh is the background, subsequent ones are
        foregrounds carved into it).
    **kwargs
        Passed through to the underlying strategy. See:
        :func:`fvcom_mesh_tools.mesh_compose.disjoint.combine_disjoint`,
        :func:`fvcom_mesh_tools.mesh_compose.overlap.combine_overlap`,
        :func:`fvcom_mesh_tools.mesh_compose.overlap.combine_neighbor`.
    """
    meshes = list(meshes)
    if len(meshes) < 2:
        raise ValueError("combine requires at least 2 meshes")
    if strategy == "disjoint":
        from fvcom_mesh_tools.mesh_compose.disjoint import combine_disjoint
        return combine_disjoint(meshes, **kwargs)
    if strategy == "overlap":
        from fvcom_mesh_tools.mesh_compose.overlap import combine_overlap
        return combine_overlap(meshes, **kwargs)
    if strategy == "neighbor":
        from fvcom_mesh_tools.mesh_compose.overlap import combine_neighbor
        return combine_neighbor(meshes, **kwargs)
    raise ValueError(f"unknown strategy: {strategy!r}; choose from {STRATEGIES}")


__all__ = ["STRATEGIES", "combine"]
