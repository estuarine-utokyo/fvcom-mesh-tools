"""OCSMesh-backed combine strategies: overlap and neighbor.

* :func:`combine_overlap` carves the second-and-later meshes into the
  first one (the "background") and stitches the seam with a buffer
  remesh. Use this for nested-resolution scenarios: a coarse outer
  mesh with one or more high-resolution inner regions.
* :func:`combine_neighbor` snaps two meshes whose edges already
  coincide (within ``snap_tol``). Use this when you generated two
  half-domains separately and need a single ``fort.14``.

Both wrap :mod:`ocsmesh.ops.combine_mesh`. Boundaries are *not*
preserved - the resulting :class:`Fort14Mesh` has empty boundary lists
and is intended to feed back into the standard fmesh post-processing
chain (boundary classification + river inflow + perpfix).
"""

from __future__ import annotations

from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_compose.convert import (
    fort14_to_meshdata,
    meshdata_to_fort14,
)


def combine_overlap(
    meshes: list[Fort14Mesh],
    *,
    crs: str = "EPSG:4326",
    adjacent_layers: int = 0,
    buffer_size: float = 0.0075,
    buffer_domain: float = 0.002,
    min_int_ang: int = 30,
    clip_final: bool = True,
    title: str | None = None,
) -> Fort14Mesh:
    """Stitch overlapping meshes via ``ocsmesh.ops.merge_overlapping_meshes``.

    The first mesh in ``meshes`` is the background; subsequent meshes
    are foregrounds that get carved into the background and re-meshed
    at the seam. Defaults follow OCSMesh's recommended settings; tune
    ``buffer_size`` and ``buffer_domain`` if the seam is too narrow or
    too wide.
    """
    from ocsmesh.ops.combine_mesh import merge_overlapping_meshes

    msht_list = [fort14_to_meshdata(m, crs=crs) for m in meshes]
    merged = merge_overlapping_meshes(
        msht_list,
        adjacent_layers=adjacent_layers,
        buffer_size=buffer_size,
        buffer_domain=buffer_domain,
        min_int_ang=min_int_ang,
        crs=crs,
        clip_final=clip_final,
    )
    if title is None:
        title = "fmesh-mesh-combine overlap: " + " + ".join(
            m.title or "(untitled)" for m in meshes
        )
    return meshdata_to_fort14(merged, title=title)


def combine_neighbor(
    meshes: list[Fort14Mesh],
    *,
    crs: str = "EPSG:4326",
    title: str | None = None,
) -> Fort14Mesh:
    """Snap meshes that share boundary nodes via ``merge_neighboring_meshes``.

    OCSMesh uses a 1e-8 (lon/lat ~= 1 mm) tolerance internally to
    identify shared vertices. Inputs whose edges do not coincide
    within that tolerance will not deduplicate cleanly - prefer
    ``combine_overlap`` for those.
    """
    from ocsmesh.ops.combine_mesh import merge_neighboring_meshes

    msht_list = [fort14_to_meshdata(m, crs=crs) for m in meshes]
    merged = merge_neighboring_meshes(*msht_list)
    if title is None:
        title = "fmesh-mesh-combine neighbor: " + " + ".join(
            m.title or "(untitled)" for m in meshes
        )
    out = meshdata_to_fort14(merged, title=title)
    if out.depths.shape[0] != out.n_nodes or float(abs(out.depths).sum()) == 0.0:
        # OCSMesh's neighbor merge does not always carry values across
        # cleanly. As a defensive fallback, populate depths by
        # interpolating each output node back to whichever input it
        # came from. The KDTree approach matches OCSMesh's own snap
        # strategy.
        import numpy as np
        from scipy.spatial import cKDTree

        all_pts = np.vstack([m.nodes for m in meshes])
        all_dep = np.concatenate([m.depths for m in meshes])
        tree = cKDTree(all_pts)
        _, idx = tree.query(out.nodes, k=1)
        out.depths = all_dep[idx]
    return out
