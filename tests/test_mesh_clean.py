"""Tests for ``fvcom_mesh_tools.mesh_clean``."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
    clean_mesh,
    keep_components,
    rebuild_boundaries,
    remove_elements,
    repair_thin_chains,
    trim_dead_ends,
    widen_thin_elements_at_centroid,
)


def _mesh(
    nodes: list[list[float]] | np.ndarray,
    elements: list[list[int]] | np.ndarray,
    *,
    open_boundaries: list[np.ndarray] | None = None,
    land_boundaries: list[tuple[int, np.ndarray]] | None = None,
) -> Fort14Mesh:
    return Fort14Mesh(
        title="test",
        nodes=np.asarray(nodes, dtype=np.float64),
        depths=np.zeros(len(nodes), dtype=np.float64),
        elements=np.asarray(elements, dtype=np.int64),
        open_boundaries=open_boundaries or [],
        land_boundaries=land_boundaries or [],
    )


def _square_around_centre() -> Fort14Mesh:
    """Unit square + central node, 4 triangles, OB on the bottom edge."""
    nodes = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
    elements = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
    return _mesh(
        nodes, elements,
        open_boundaries=[np.array([0, 1])],
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )


def _square_plus_disjoint_triangle() -> Fort14Mesh:
    """``_square_around_centre`` plus an isolated triangle far to the east.

    The OB is only on the square; the isolated triangle has no OB node.
    """
    base = _square_around_centre()
    extra_nodes = np.array([[10.0, 0.0], [11.0, 0.0], [10.5, 1.0]])
    nodes = np.vstack([base.nodes, extra_nodes])
    elements = np.vstack([base.elements, np.array([[5, 6, 7]])])
    return Fort14Mesh(
        title=base.title,
        nodes=nodes,
        depths=np.zeros(len(nodes), dtype=np.float64),
        elements=elements,
        open_boundaries=base.open_boundaries,
        land_boundaries=base.land_boundaries,
    )


# ---------------------------------------------------------------------------
# remove_elements
# ---------------------------------------------------------------------------


def test_remove_elements_drops_unused_nodes_and_renumbers() -> None:
    """Two adjacent triangles. Drop the second; its private vertex
    (node 3) must disappear, and the surviving element must reference
    the renumbered nodes 0..2."""
    nodes = [[0, 0], [1, 0], [0.5, 1], [1.5, 1]]
    elements = [[0, 1, 2], [1, 3, 2]]
    mesh = _mesh(nodes, elements)
    out = remove_elements(mesh, np.array([True, False]))
    assert out.n_elements == 1
    assert out.n_nodes == 3
    np.testing.assert_array_equal(out.elements, [[0, 1, 2]])
    np.testing.assert_array_equal(out.nodes, [[0, 0], [1, 0], [0.5, 1]])
    assert out.open_boundaries == []
    assert out.land_boundaries == []


def test_remove_elements_empty_result_returns_empty_mesh() -> None:
    mesh = _square_around_centre()
    out = remove_elements(mesh, np.zeros(mesh.n_elements, dtype=bool))
    assert out.n_elements == 0
    assert out.n_nodes == 0


# ---------------------------------------------------------------------------
# keep_components
# ---------------------------------------------------------------------------


def test_keep_components_default_keeps_only_largest() -> None:
    mesh = _square_plus_disjoint_triangle()
    cleaned, info = keep_components(mesh)
    assert cleaned.n_elements == 4
    assert cleaned.n_nodes == 5
    assert info["n_components_before"] == 2
    assert info["n_components_kept"] == 1
    assert info["kept_component_sizes"] == [4]
    assert info["n_elements_removed"] == 1
    assert info["n_nodes_removed"] == 3


def test_keep_components_min_elements_drops_small_pools() -> None:
    """Three components: 4, 1, 1. ``min_elements=2`` keeps only the 4-elem one."""
    sq = _square_around_centre()
    extra = np.array([
        [10, 0], [11, 0], [10.5, 1],   # second component
        [20, 0], [21, 0], [20.5, 1],   # third component
    ], dtype=np.float64)
    nodes = np.vstack([sq.nodes, extra])
    elements = np.vstack([
        sq.elements, np.array([[5, 6, 7], [8, 9, 10]], dtype=np.int64),
    ])
    mesh = _mesh(nodes, elements, open_boundaries=sq.open_boundaries)
    cleaned, info = keep_components(mesh, min_elements=2)
    assert cleaned.n_elements == 4
    assert info["n_components_before"] == 3
    assert info["kept_component_sizes"] == [4]


def test_keep_components_require_open_boundary() -> None:
    """Two components; the ONLY OB node is on the smaller one. Requiring
    OB-touch keeps the smaller component, not the larger.
    """
    sq = _square_around_centre()
    extra_nodes = np.array([[10, 0], [11, 0], [10.5, 1]])
    nodes = np.vstack([sq.nodes, extra_nodes])
    elements = np.vstack([sq.elements, np.array([[5, 6, 7]])])
    # Move the OB onto the small component only.
    mesh = _mesh(nodes, elements, open_boundaries=[np.array([5, 6])])
    cleaned, info = keep_components(mesh, require_open_boundary=True)
    assert cleaned.n_elements == 1
    assert info["kept_component_sizes"] == [1]


def test_keep_components_falls_back_to_largest_when_filter_empties() -> None:
    """``min_elements`` larger than every component would empty the mesh;
    safety net keeps the largest component instead.
    """
    mesh = _square_plus_disjoint_triangle()
    cleaned, info = keep_components(mesh, min_elements=999)
    assert cleaned.n_elements == 4
    assert info["n_components_kept"] == 1


# ---------------------------------------------------------------------------
# rebuild_boundaries
# ---------------------------------------------------------------------------


def test_rebuild_boundaries_reclassifies_outer_ring() -> None:
    """Provide a bbox that picks the bottom side as open. The classifier
    must produce one open segment along the bottom and a land segment
    covering the rest of the ring.
    """
    mesh = _square_around_centre()
    # Strip boundaries first, then rebuild from scratch.
    blank = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes,
        depths=mesh.depths,
        elements=mesh.elements,
        open_boundaries=[],
        land_boundaries=[],
    )
    rebuilt = rebuild_boundaries(
        blank, bbox=(0.0, 0.0, 1.0, 0.0), tol_deg=1e-6, land_ibtype=20,
    )
    assert len(rebuilt.open_boundaries) == 1
    assert len(rebuilt.land_boundaries) == 1
    open_seg = rebuilt.open_boundaries[0]
    # Open segment lives along y=0; both nodes 0 and 1 must appear.
    assert {0, 1}.issubset(set(int(i) for i in open_seg))


# ---------------------------------------------------------------------------
# trim_dead_ends
# ---------------------------------------------------------------------------


def test_trim_dead_ends_no_op_on_clean_mesh() -> None:
    mesh = _square_around_centre()
    out, info = trim_dead_ends(
        mesh, max_iters=5,
        bbox=(0.0, 0.0, 1.0, 0.0), tol_deg=1e-6, land_ibtype=20,
    )
    assert out.n_elements == 4
    assert info["iterations_run"] in (0, 1)
    assert info["total_elements_removed"] == 0
    assert info["converged"]


def test_trim_dead_ends_removes_dangling_pair() -> None:
    """Two adjacent triangles with no OB; both are degree-1 and have no
    OB edge, so trimming should empty the mesh in one round."""
    nodes = np.array(
        [[0, 0], [1, 0], [0.5, 1], [1.5, 1]],
        dtype=np.float64,
    )
    elements = np.array(
        [[0, 1, 2], [1, 3, 2]],
        dtype=np.int64,
    )
    mesh = _mesh(nodes, elements,
                 land_boundaries=[(0, np.array([0, 1, 3, 2, 0]))])
    out, info = trim_dead_ends(
        mesh, max_iters=5,
        bbox=(0.0, 0.0, 0.5, 0.0), tol_deg=1e-6, land_ibtype=20,
    )
    # Both starting elements were degree-1 and had no OB; after one round
    # both are removed and the mesh is empty.
    assert out.n_elements == 0
    assert info["per_iter_dead_end_count"][0] == 2


# ---------------------------------------------------------------------------
# clean_mesh integration
# ---------------------------------------------------------------------------


def test_clean_mesh_removes_disjoint_and_rebuilds_boundaries() -> None:
    mesh = _square_plus_disjoint_triangle()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 1.0, 0.0),  # bottom edge of the big square
        bbox_tol_m=10.0,            # small tol; main square's y=0 nodes are close
        land_ibtype=20,
        trim_dead_ends_iters=0,
    )
    assert info["input"]["n_elements"] == 5
    assert info["output"]["n_elements"] == 4
    assert info["output"]["n_nodes"] == 5
    # rebuild produced at least one boundary segment.
    assert (
        info["output"]["n_open_boundaries"]
        + info["output"]["n_land_boundaries"]
    ) >= 1
    # First (and only) phase recorded.
    assert info["phases"][0]["name"] == "keep_components"
    assert info["phases"][0]["n_elements_removed"] == 1


def test_clean_mesh_pass_through_when_all_phases_disabled() -> None:
    mesh = _square_around_centre()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 1.0, 1.0),
        bbox_tol_m=1.0,
        remove_disjoint=False,
        trim_dead_ends_iters=0,
        thin_chain_mode="none",
    )
    # No element deletion, but boundaries are normalised.
    assert info["output"]["n_elements"] == info["input"]["n_elements"]
    assert info["phases"] == []


# ---------------------------------------------------------------------------
# Phase C: widen / delete thin chains
# ---------------------------------------------------------------------------


def _thin_chain_strip(chain_length: int) -> Fort14Mesh:
    """``chain_length`` quads forming a 1-cell-wide channel; every node
    on the top or bottom polyline is on a land-boundary segment, so all
    triangles in the strip are thin.
    """
    n = chain_length + 1
    bot = np.column_stack([np.arange(n, dtype=np.float64), np.zeros(n)])
    top = np.column_stack([np.arange(n, dtype=np.float64), np.ones(n)])
    nodes = np.vstack([bot, top])
    elems = []
    for i in range(chain_length):
        elems.append([i, i + 1, n + i])
        elems.append([i + 1, n + i + 1, n + i])
    elements = np.asarray(elems, dtype=np.int64)
    bot_seg = np.arange(n, dtype=np.int64)
    top_seg = np.arange(n, dtype=np.int64) + n
    return _mesh(
        nodes, elements,
        land_boundaries=[(0, bot_seg), (0, top_seg)],
    )


def test_widen_thin_elements_inserts_centroid_per_target() -> None:
    """A single thin triangle: widening adds 1 node + 2 net elements."""
    mesh = _thin_chain_strip(chain_length=1)   # 1 quad, 2 thin triangles
    flag = np.array([True, False])             # widen only the first
    out, info = widen_thin_elements_at_centroid(mesh, flag)
    assert info["n_widened"] == 1
    assert info["n_new_nodes"] == 1
    assert info["n_new_elements"] == 2
    assert out.n_nodes == mesh.n_nodes + 1
    assert out.n_elements == mesh.n_elements + 2
    # The new node is the centroid of the original first triangle.
    centroid = mesh.nodes[mesh.elements[0]].mean(axis=0)
    np.testing.assert_allclose(out.nodes[-1], centroid, atol=1e-12)
    # Boundary lists preserved (centroid is interior).
    assert len(out.land_boundaries) == len(mesh.land_boundaries)


def test_repair_thin_chains_widen_replaces_every_chain_element() -> None:
    """Strip of 8 thin elements (chain length=4 quads), all in one chain."""
    mesh = _thin_chain_strip(chain_length=4)   # 8 thin triangles
    out, info = repair_thin_chains(mesh, mode="widen", min_chain_length=3)
    assert info["mode"] == "widen"
    assert info["n_chain_elements"] == 8
    assert info["n_widened"] == 8
    assert out.n_nodes == mesh.n_nodes + 8
    # 8 thin triangles -> 24 sub-triangles + 0 originals = NE+16
    assert out.n_elements == mesh.n_elements + 16


def test_repair_thin_chains_delete_removes_every_chain_element() -> None:
    mesh = _thin_chain_strip(chain_length=4)
    out, info = repair_thin_chains(
        mesh, mode="delete", min_chain_length=3,
        bbox=(0.0, 0.0, 4.0, 0.0), tol_deg=1e-6, land_ibtype=0,
    )
    assert info["mode"] == "delete"
    assert info["n_chain_elements"] == 8
    assert info["n_elements_removed"] == 8
    assert out.n_elements == 0       # nothing left


def test_repair_thin_chains_none_is_noop() -> None:
    mesh = _thin_chain_strip(chain_length=4)
    out, info = repair_thin_chains(mesh, mode="none")
    assert info["mode"] == "none"
    assert out.n_elements == mesh.n_elements
    assert out.n_nodes == mesh.n_nodes


def test_repair_thin_chains_skips_short_chains() -> None:
    """A chain of length 2 (4 thin triangles) below threshold 5 → noop."""
    mesh = _thin_chain_strip(chain_length=2)   # 4 thin triangles
    out, info = repair_thin_chains(mesh, mode="widen", min_chain_length=5)
    assert info["n_chain_elements"] == 0
    assert info["skipped"] is True
    assert out.n_elements == mesh.n_elements


def test_clean_mesh_phase_c_widen_default() -> None:
    """Full pipeline default mode is 'widen'; thin chain in input gets widened."""
    mesh = _thin_chain_strip(chain_length=4)
    out, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 4.0, 0.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        # thin_chain_mode defaults to widen
        min_thin_chain=3,
    )
    phase_c = next(p for p in info["phases"] if p["name"] == "repair_thin_chains")
    assert phase_c["mode"] == "widen"
    assert phase_c["n_widened"] == 8
    assert info["output"]["n_nodes"] > info["input"]["n_nodes"]


def test_clean_mesh_phase_c_delete_explicit() -> None:
    mesh = _thin_chain_strip(chain_length=4)
    out, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 4.0, 0.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="delete", min_thin_chain=3,
    )
    phase_c = next(p for p in info["phases"] if p["name"] == "repair_thin_chains")
    assert phase_c["mode"] == "delete"
    assert phase_c["n_elements_removed"] == 8
    assert info["output"]["n_elements"] == 0
