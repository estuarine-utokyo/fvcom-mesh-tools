"""Boundary extraction and open / land classification for triangular meshes.

The classifier here is geared at meshes produced by the minimal OCSMesh
pipeline (PoC #5), where the mesh polygon's outer ring consists of two
kinds of segments:

* Coast: where the input DEM crosses ``z = zmax`` and gmsh follows the
  zero-elevation contour.
* Open ocean: where the DEM is clipped to its own bounding box and the
  mesh edge runs along that rectangular limit.

Splitting the outer ring on "is this node within ``tol_deg`` of the DEM
bounding box?" turns out to be a robust, hyperparameter-light heuristic
for real meshes. The same module exposes the lower-level boundary-edge
walker so other classifiers can be layered on top.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

# ---------------------------------------------------------------------------
# Edge / loop extraction
# ---------------------------------------------------------------------------


def boundary_edges_from_tris(elements: np.ndarray) -> np.ndarray:
    """Return ``(K, 2)`` int array of undirected boundary edges.

    A boundary edge belongs to exactly one triangle. ``elements`` must be
    ``(NE, 3)`` of 0-indexed node indices.
    """
    if elements.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    e = np.vstack(
        [elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]]]
    ).astype(np.int64)
    e.sort(axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    return uniq[cnt == 1].copy()


def chain_edges_to_loops(edges: np.ndarray) -> list[np.ndarray]:
    """Chain undirected boundary edges into closed loops of node indices.

    Each returned loop is an array ``[n0, n1, ..., nk, n0]`` (closed by
    repeating the first node at the end). Open chains are dropped, since
    a watertight mesh's boundary edges always form closed cycles.
    """
    if edges.size == 0:
        return []

    nbrs: dict[int, list[int]] = defaultdict(list)
    for i, j in edges.tolist():
        nbrs[int(i)].append(int(j))
        nbrs[int(j)].append(int(i))

    used: set[tuple[int, int]] = set()

    def ek(i: int, j: int) -> tuple[int, int]:
        return (i, j) if i < j else (j, i)

    loops: list[np.ndarray] = []
    for start_node in sorted(nbrs.keys()):
        for first_step in sorted(nbrs[start_node]):
            if ek(start_node, first_step) in used:
                continue
            # Begin a new loop walk from start_node -> first_step.
            path = [start_node, first_step]
            used.add(ek(start_node, first_step))
            prev, cur = start_node, first_step
            while cur != start_node:
                cands = [v for v in nbrs[cur] if v != prev and ek(cur, v) not in used]
                if not cands:
                    # Dead-end (non-manifold mesh) or already-closed elsewhere.
                    break
                # Deterministic choice: smallest unused neighbour.
                nxt = min(cands)
                used.add(ek(cur, nxt))
                path.append(nxt)
                prev, cur = cur, nxt
            if len(path) >= 4 and path[0] == path[-1]:
                loops.append(np.asarray(path, dtype=np.int64))
    return loops


def _ring_signed_area(loop: np.ndarray, nodes: np.ndarray) -> float:
    xy = nodes[loop[:-1]]  # drop repeated closing node
    x = xy[:, 0]
    y = xy[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def outer_loop(loops: list[np.ndarray], nodes: np.ndarray) -> np.ndarray:
    """Pick the loop with maximum absolute signed area as the outer ring.

    Raises :class:`ValueError` if ``loops`` is empty.
    """
    if not loops:
        raise ValueError("no boundary loops to choose from")
    areas = [abs(_ring_signed_area(loop, nodes)) for loop in loops]
    return loops[int(np.argmax(areas))]


# ---------------------------------------------------------------------------
# Open / land classification by DEM bounding box proximity
# ---------------------------------------------------------------------------


def _on_bbox(
    pts: np.ndarray,
    bbox: tuple[float, float, float, float],
    tol: float,
) -> np.ndarray:
    """Per-row mask: True if the point lies within ``tol`` of the rectangle."""
    xmin, ymin, xmax, ymax = bbox
    x = pts[:, 0]
    y = pts[:, 1]
    near_x = (np.abs(x - xmin) <= tol) | (np.abs(x - xmax) <= tol)
    near_y = (np.abs(y - ymin) <= tol) | (np.abs(y - ymax) <= tol)
    in_x = (x >= xmin - tol) & (x <= xmax + tol)
    in_y = (y >= ymin - tol) & (y <= ymax + tol)
    return (near_x & in_y) | (near_y & in_x)


def _runs(mask: np.ndarray) -> list[tuple[int, int, bool]]:
    """Run-length encode a 1-D bool array.

    Returns a list of ``(start, end, value)`` with ``end`` inclusive,
    over the index range ``[0, len(mask) - 1]``.
    """
    if mask.size == 0:
        return []
    out: list[tuple[int, int, bool]] = []
    i = 0
    while i < mask.size:
        v = bool(mask[i])
        j = i
        while j + 1 < mask.size and bool(mask[j + 1]) == v:
            j += 1
        out.append((i, j, v))
        i = j + 1
    return out


def _bridge_short_coast_gaps(is_open: np.ndarray, max_gap: int) -> np.ndarray:
    """Reclassify any "land" run shorter than ``max_gap`` as "open" if it
    sits between two open runs (cyclic).

    Operates on a 1-D bool mask whose True entries mean "open". Returns a
    copy with the bridging applied. ``max_gap=0`` returns the input
    unchanged.
    """
    if max_gap <= 0 or is_open.size == 0:
        return is_open.copy()
    if bool(is_open.all()) or not bool(is_open.any()):
        return is_open.copy()
    out = is_open.copy()
    runs = _runs(out)
    n = len(runs)
    for k, (i, j, v) in enumerate(runs):
        if v:
            continue
        run_len = j - i + 1
        if run_len > max_gap:
            continue
        prev_v = runs[(k - 1) % n][2]
        next_v = runs[(k + 1) % n][2]
        if prev_v and next_v:
            out[i : j + 1] = True
    return out


def classify_outer_loop_by_bbox(
    outer: np.ndarray,
    nodes: np.ndarray,
    bbox: tuple[float, float, float, float],
    tol: float,
    *,
    open_merge_coast_gap: int = 0,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Split a closed outer loop into open- and land-boundary segments.

    A boundary node is classified as *open* if it lies within ``tol`` of
    the rectangle ``bbox = (xmin, ymin, xmax, ymax)``. Consecutive same-
    class nodes form one segment; the loop is rotated so the first
    segment starts at a class transition (so segments do not span the
    arbitrary closing index).

    With ``open_merge_coast_gap > 0``, any *land* run shorter than this
    threshold (in nodes) sandwiched between two open runs is reclassified
    as open. This collapses what is geometrically one open arc broken by
    a small coast intrusion into one segment, matching the convention
    used by typical FVCOM regional meshes.

    Returns
    -------
    (open_segments, land_segments)
        Lists of 0-indexed node-index arrays. Each array represents a
        polyline (not closed). Empty lists if all of one kind.
    """
    if outer.size < 4 or outer[0] != outer[-1]:
        raise ValueError("outer must be a closed loop with >= 3 distinct nodes")
    ring = outer[:-1]  # drop closing duplicate
    pts = nodes[ring]
    is_open = _on_bbox(pts, bbox, tol)

    if open_merge_coast_gap > 0:
        is_open = _bridge_short_coast_gaps(is_open, open_merge_coast_gap)

    if bool(is_open.all()):
        return [ring.copy()], []
    if not bool(is_open.any()):
        return [], [ring.copy()]

    # Rotate so index 0 is at a class boundary (transition in -> out or
    # out -> in). Otherwise a segment can wrap across the array end.
    transitions = np.where(is_open != np.roll(is_open, 1))[0]
    rot = int(transitions[0])
    ring = np.roll(ring, -rot)
    is_open = np.roll(is_open, -rot)

    open_segments: list[np.ndarray] = []
    land_segments: list[np.ndarray] = []
    for i, j, v in _runs(is_open):
        seg = ring[i : j + 1].copy()
        # Append the first node of the next run so adjacent segments share
        # one node; this matches FVCOM/ADCIRC boundary-list convention
        # where land and open segments meet at a shared endpoint.
        next_idx = (j + 1) % ring.size
        seg = np.concatenate([seg, ring[next_idx : next_idx + 1]])
        if v:
            open_segments.append(seg)
        else:
            land_segments.append(seg)
    return open_segments, land_segments


# ---------------------------------------------------------------------------
# High-level: classify a Fort14Mesh against a DEM bounding box
# ---------------------------------------------------------------------------


def classify_boundaries_by_bbox(
    mesh: Fort14Mesh,
    bbox: tuple[float, float, float, float],
    tol: float,
    land_ibtype: int = 0,
    *,
    open_merge_coast_gap: int = 0,
) -> tuple[list[np.ndarray], list[tuple[int, np.ndarray]]]:
    """Compute open- and land-boundary segments for ``mesh``.

    Parameters
    ----------
    mesh:
        Mesh whose boundaries are not yet populated (typically straight
        out of OCSMesh).
    bbox:
        ``(xmin, ymin, xmax, ymax)`` in the same coordinate system as
        ``mesh.nodes``. Boundary nodes within ``tol`` of this rectangle
        are flagged as open.
    tol:
        Distance tolerance, in the units of ``mesh.nodes`` (degrees for
        lon/lat). A few times the local edge length usually works.
    land_ibtype:
        Integer ibtype written to fort.14 for every land segment.
        ADCIRC convention is 0 for normal flow / coast; FVCOM regional
        meshes sometimes use 20 / 21.

    Returns
    -------
    (open_segments, land_boundaries)
        Lists shaped to drop straight into
        ``Fort14Mesh.open_boundaries`` / ``Fort14Mesh.land_boundaries``.
        Inner-loop holes (islands) are emitted as land segments only.
    """
    edges = boundary_edges_from_tris(mesh.elements)
    loops = chain_edges_to_loops(edges)
    if not loops:
        return [], []

    outer = outer_loop(loops, mesh.nodes)
    open_segs, land_segs = classify_outer_loop_by_bbox(
        outer, mesh.nodes, bbox, tol,
        open_merge_coast_gap=open_merge_coast_gap,
    )

    # Hole loops (islands) are always coast.
    for loop in loops:
        if loop is outer:
            continue
        # Drop the closing duplicate; islands are stored as polygons
        # whose first and last node ids match in fort.14 too, but the
        # convention in this repo is to store polylines without the
        # repeated closing node.
        land_segs.append(loop[:-1].copy())

    land_boundaries = [(int(land_ibtype), seg) for seg in land_segs]
    return open_segs, land_boundaries
