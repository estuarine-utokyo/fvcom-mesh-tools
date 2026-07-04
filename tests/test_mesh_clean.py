"""Tests for ``fvcom_mesh_tools.mesh_clean``."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
    analyze_under_resolved_channels,
    clean_mesh,
    keep_components,
    rebuild_boundaries,
    remove_elements,
    repair_overconnected_nodes,
    repair_skewed_elements,
    repair_thin_chains,
    repair_under_resolved_channels,
    smooth_mesh_laplacian,
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


# ---------------------------------------------------------------------------
# Phase D: over-connected node repair
# ---------------------------------------------------------------------------


def _fan_mesh(n_wedges: int) -> Fort14Mesh:
    """Regular ``n_wedges``-spoke fan around an interior centre node.

    Centre node valence = ``n_wedges``; every rim node has valence 2.
    The rim is a single closed land segment.
    """
    centre = np.array([[0.0, 0.0]])
    angles = np.linspace(0.0, 2.0 * np.pi, n_wedges + 1)[:-1]
    rim = np.column_stack([np.cos(angles), np.sin(angles)])
    nodes = np.vstack([centre, rim])
    elements = np.array(
        [[0, 1 + i, 1 + (i + 1) % n_wedges] for i in range(n_wedges)],
        dtype=np.int64,
    )
    rim_chain = np.concatenate([np.arange(1, n_wedges + 1), [1]])
    return _mesh(nodes, elements, land_boundaries=[(0, rim_chain.astype(np.int64))])


def test_repair_overconnected_nodes_floor20_rejects_fan_flips() -> None:
    """The 12-wedge fan has a centre at valence 12; with the FVCOM-safe
    20° quality floor every flip would create a sliver chord triangle
    so the algorithm is expected to reject every candidate.
    """
    mesh = _fan_mesh(12)
    out, info = repair_overconnected_nodes(
        mesh, max_nbr_elem=8, max_iters=10, min_angle_floor_deg=20.0,
    )
    assert info["max_valence_before"] == 12
    assert info["max_valence_after"] == 12
    assert info["total_swaps"] == 0


def test_repair_overconnected_nodes_floor0_balances_fan() -> None:
    """With floor=0 the algorithm is allowed to introduce slivers and
    can drive the centre valence below the cap.
    """
    mesh = _fan_mesh(12)
    out, info = repair_overconnected_nodes(
        mesh, max_nbr_elem=8, max_iters=10, min_angle_floor_deg=0.0,
    )
    assert info["max_valence_before"] == 12
    assert info["max_valence_after"] <= 8
    assert info["total_swaps"] >= 1
    assert info["n_overconn_after"] == 0


def test_repair_overconnected_nodes_noop_when_cap_already_met() -> None:
    """6-wedge fan: centre valence 6, already <= 8, so the algorithm is
    a no-op irrespective of floor.
    """
    mesh = _fan_mesh(6)
    out, info = repair_overconnected_nodes(
        mesh, max_nbr_elem=8, max_iters=10, min_angle_floor_deg=20.0,
    )
    assert info["total_swaps"] == 0
    assert info["max_valence_before"] == 6
    assert info["max_valence_after"] == 6
    assert info["n_overconn_before"] == 0


def test_clean_mesh_phase_d_default_off() -> None:
    """Phase D is off by default — the 12-wedge fan stays at valence 12
    after a default mesh-clean run.
    """
    mesh = _fan_mesh(12)
    cleaned, info = clean_mesh(
        mesh,
        bbox=(-1.0, -1.0, 1.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        # Default: repair_overconnected_iters=0
    )
    assert all(p["name"] != "repair_overconnected_nodes" for p in info["phases"])


def test_clean_mesh_phase_d_explicit_floor0_balances_fan() -> None:
    """Enabling Phase D with floor=0 reduces fan centre valence."""
    mesh = _fan_mesh(12)
    cleaned, info = clean_mesh(
        mesh,
        bbox=(-1.0, -1.0, 1.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        repair_overconnected_iters=10,
        max_nbr_elem=8,
        overconn_min_angle_floor_deg=0.0,
    )
    phase_d = next(p for p in info["phases"]
                   if p["name"] == "repair_overconnected_nodes")
    assert phase_d["max_valence_after"] <= 8
    assert phase_d["n_overconn_after"] == 0


# ---------------------------------------------------------------------------
# Phase E: under-resolved channel repair (detector 6)
# ---------------------------------------------------------------------------


def _strip_mesh_n_rows(
    *, length_deg: float, width_deg: float, n_x: int, n_rows: int,
) -> Fort14Mesh:
    """Strip in lon/lat: top + bottom land, left + right open boundary.

    Same fixture as in ``test_diagnostics`` — duplicated here so the
    two test files do not depend on each other.
    """
    n_pts_x = n_x + 1
    n_pts_y = n_rows + 1
    nodes_rows: list[np.ndarray] = []
    for j in range(n_pts_y):
        y = (j / (n_pts_y - 1)) * width_deg
        nodes_rows.append(np.column_stack([
            np.linspace(0.0, length_deg, n_pts_x),
            np.full(n_pts_x, y),
        ]))
    nodes = np.vstack(nodes_rows)
    elems: list[list[int]] = []
    for j in range(n_rows):
        for i in range(n_x):
            i00 = j * n_pts_x + i
            i01 = j * n_pts_x + (i + 1)
            i10 = (j + 1) * n_pts_x + i
            i11 = (j + 1) * n_pts_x + (i + 1)
            elems.append([i00, i01, i11])
            elems.append([i00, i11, i10])
    elements = np.asarray(elems, dtype=np.int64)
    bot_seg = np.arange(n_pts_x, dtype=np.int64)
    top_seg = np.arange(n_pts_x, dtype=np.int64) + n_rows * n_pts_x
    left_seg = np.array(
        [j * n_pts_x for j in range(n_pts_y)], dtype=np.int64,
    )
    right_seg = np.array(
        [j * n_pts_x + (n_pts_x - 1) for j in range(n_pts_y)], dtype=np.int64,
    )
    return _mesh(
        nodes, elements,
        open_boundaries=[left_seg, right_seg],
        land_boundaries=[(0, bot_seg), (0, top_seg)],
    )


def test_repair_under_resolved_channels_widen_inserts_centroids() -> None:
    """A 1-row strip is fully flagged by detector 6 (w/h ≈ 1); widening
    inserts one centroid per element."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    out, info = repair_under_resolved_channels(
        mesh, mode="widen", min_w_h=3.0,
    )
    assert info["mode"] == "widen"
    assert info["n_flagged"] == mesh.n_elements
    assert info["n_widened"] == mesh.n_elements
    assert out.n_nodes == mesh.n_nodes + mesh.n_elements
    # Each flagged element becomes 3 sub-triangles → net +2 per flagged.
    assert out.n_elements == 3 * mesh.n_elements


def test_repair_under_resolved_channels_delete_removes_flagged() -> None:
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    out, info = repair_under_resolved_channels(
        mesh, mode="delete", min_w_h=3.0,
        bbox=(0.0, 0.0, 0.08, 0.01), tol_deg=1e-3, land_ibtype=0,
    )
    assert info["mode"] == "delete"
    assert info["n_flagged"] == mesh.n_elements
    assert info["n_elements_removed"] == mesh.n_elements
    assert out.n_elements == 0


def test_repair_under_resolved_channels_none_is_noop() -> None:
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    out, info = repair_under_resolved_channels(mesh, mode="none")
    assert info["mode"] == "none"
    assert out.n_elements == mesh.n_elements
    assert out.n_nodes == mesh.n_nodes


def test_repair_under_resolved_channels_skips_when_unflagged() -> None:
    """A wide 6-row strip has at least some elements above min_w_h=3 —
    relaxing the threshold below 1 marks zero elements (no flags)."""
    mesh = _strip_mesh_n_rows(
        length_deg=0.5, width_deg=0.10, n_x=20, n_rows=6,
    )
    out, info = repair_under_resolved_channels(
        mesh, mode="widen", min_w_h=0.1,
    )
    assert info["n_flagged"] == 0
    assert info.get("skipped") is True
    assert out.n_nodes == mesh.n_nodes
    assert out.n_elements == mesh.n_elements


def test_repair_under_resolved_channels_delete_requires_bbox() -> None:
    """``mode='delete'`` without bbox/tol_deg raises a clear error."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    import pytest

    with pytest.raises(ValueError, match="delete mode requires bbox"):
        repair_under_resolved_channels(mesh, mode="delete", min_w_h=3.0)


def test_repair_under_resolved_channels_medial_one_row_strip() -> None:
    """Stage 2: a 1-row strip's single channel is replaced by a Delaunay
    of (rim ∪ spine). The new mesh has more nodes (spine inserts) and
    more elements (the patch is now 2 cells across at most points)."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    out, info = repair_under_resolved_channels(
        mesh, mode="medial", min_w_h=3.0, min_channel_elements=1,
    )
    assert info["mode"] == "medial"
    assert info["n_components"] == 1
    assert info["n_components_replaced"] == 1
    assert info["n_components_skipped"] == 0
    assert info["n_nodes_inserted"] > 0
    assert info["n_elements_removed"] == mesh.n_elements
    assert info["n_elements_inserted"] >= mesh.n_elements
    # Original rim node IDs stay at the same positions.
    np.testing.assert_array_equal(out.nodes[: mesh.n_nodes], mesh.nodes)
    # Boundary lists are carried forward by node ID.
    assert len(out.open_boundaries) == 2
    assert len(out.land_boundaries) == 2
    # No flipped triangles in the output.
    n0, n1, n2 = (
        out.elements[:, 0], out.elements[:, 1], out.elements[:, 2],
    )
    p0, p1, p2 = out.nodes[n0], out.nodes[n1], out.nodes[n2]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all(), "Stage 2 produced a flipped or degenerate element"


def test_repair_under_resolved_channels_medial_skips_short_components() -> None:
    """``min_channel_elements`` larger than the strip's one component
    drops every flag, leaving the mesh unchanged."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    out, info = repair_under_resolved_channels(
        mesh, mode="medial", min_w_h=3.0, min_channel_elements=999,
    )
    assert info["n_flagged"] == 0
    assert info.get("skipped") is True
    assert out.n_nodes == mesh.n_nodes
    assert out.n_elements == mesh.n_elements


def test_repair_under_resolved_channels_medial_empty_flag_is_noop() -> None:
    """No flagged elements (high min_w_h floor) → mesh unchanged."""
    mesh = _strip_mesh_n_rows(
        length_deg=0.5, width_deg=0.10, n_x=20, n_rows=6,
    )
    out, info = repair_under_resolved_channels(
        mesh, mode="medial", min_w_h=0.1,
    )
    assert info["n_flagged"] == 0
    assert info.get("skipped") is True
    assert out.n_nodes == mesh.n_nodes
    assert out.n_elements == mesh.n_elements


def test_repair_under_resolved_channels_medial_idempotent_on_repair() -> None:
    """Running medial mode twice is a no-op on the second pass: the
    first pass widens the channel enough that detector 6 no longer
    flags it (or flags only short residual clusters)."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    once, _ = repair_under_resolved_channels(
        mesh, mode="medial", min_w_h=3.0, min_channel_elements=1,
    )
    twice, info2 = repair_under_resolved_channels(
        once, mode="medial", min_w_h=3.0, min_channel_elements=1,
    )
    # Either no further flagging, or the second pass is a much smaller
    # change than the first (no runaway insertion).
    if info2.get("skipped"):
        assert twice.n_nodes == once.n_nodes
        assert twice.n_elements == once.n_elements
    else:
        assert info2["n_nodes_inserted"] < info2.get(
            "n_components_replaced", 0,
        ) * once.n_elements + 1


def test_clean_mesh_phase_e_default_off() -> None:
    """Phase E is off by default — a fully-flagged 1-row strip is left
    untouched after a default ``clean_mesh`` run with all other phases
    disabled."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 0.08, 0.01), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        # under_resolved_mode defaults to "none"
    )
    assert all(p["name"] != "repair_under_resolved_channels"
               for p in info["phases"])
    assert cleaned.n_elements == mesh.n_elements


def test_clean_mesh_phase_e_widen_explicit() -> None:
    """Enabling ``under_resolved_mode='widen'`` widens every flagged
    element in the 1-row strip."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 0.08, 0.01), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        under_resolved_mode="widen",
        under_resolved_min_w_h=3.0,
    )
    phase_e = next(p for p in info["phases"]
                   if p["name"] == "repair_under_resolved_channels")
    assert phase_e["mode"] == "widen"
    assert phase_e["n_flagged"] == mesh.n_elements
    assert phase_e["n_widened"] == mesh.n_elements
    assert info["output"]["n_nodes"] == mesh.n_nodes + mesh.n_elements


def test_clean_mesh_phase_e_min_channel_elements_filter_propagates() -> None:
    """``under_resolved_min_channel_elements`` larger than the strip's
    one component drops every flag, so Phase E reports n_flagged=0."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 0.08, 0.01), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        under_resolved_mode="widen",
        under_resolved_min_w_h=3.0,
        under_resolved_min_channel_elements=999,
    )
    phase_e = next(p for p in info["phases"]
                   if p["name"] == "repair_under_resolved_channels")
    assert phase_e["min_channel_elements"] == 999
    assert phase_e["n_flagged"] == 0
    assert phase_e.get("skipped") is True
    assert info["output"]["n_nodes"] == mesh.n_nodes
    assert info["output"]["n_elements"] == mesh.n_elements


def test_clean_mesh_phase_e_invalid_mode_raises() -> None:
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    import pytest

    with pytest.raises(ValueError, match="under_resolved_mode"):
        clean_mesh(
            mesh,
            bbox=(0.0, 0.0, 0.08, 0.01), bbox_tol_m=1.0,
            under_resolved_mode="bogus",  # type: ignore[arg-type]
        )


# ---------------------------------------------------------------------------
# Phase F: angle-based skewed-element removal (wraps ocsmesh)
# ---------------------------------------------------------------------------


def _mesh_with_one_sliver() -> Fort14Mesh:
    """Two unit triangles + one extreme sliver attached to the strip.

    The sliver has interior angles of approximately (179.7°, 0.15°,
    0.15°) — well outside the default Phase F bounds [1°, 175°]. The
    other two triangles are 45-45-90 right triangles, well inside the
    default bounds.
    """
    nodes = np.array([
        [0.0, 0.0],     # 0
        [1.0, 0.0],     # 1
        [1.0, 1.0],     # 2
        [0.0, 1.0],     # 3
        [2.0, 0.0],     # 4 — sliver tip far off
        [2.0, 0.001],   # 5 — sliver tip very close to 4
    ], dtype=float)
    elements = np.array([
        [0, 1, 2],     # right tri (good)
        [0, 2, 3],     # right tri (good)
        [1, 4, 5],     # sliver: max angle ~179.7°
    ], dtype=np.int64)
    return _mesh(
        nodes, elements,
        open_boundaries=[np.array([0, 3])],
        land_boundaries=[(0, np.array([3, 2, 1, 4, 5]))],
    )


def test_repair_skewed_elements_removes_sliver() -> None:
    mesh = _mesh_with_one_sliver()
    out, info = repair_skewed_elements(
        mesh,
        min_angle_deg=1.0,
        max_angle_deg=175.0,
        bbox=(0.0, 0.0, 2.0, 1.0),
        tol_deg=1e-3,
        land_ibtype=0,
    )
    assert info["n_elements_removed"] == 1
    assert out.n_elements == mesh.n_elements - 1


def test_repair_skewed_elements_noop_preserves_boundaries() -> None:
    """A clean mesh: no element removed, original boundary lists kept."""
    mesh = _mesh_with_one_sliver()
    # Drop the sliver up front so the input has none.
    mesh = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes,
        depths=mesh.depths,
        elements=mesh.elements[:2].copy(),
        open_boundaries=mesh.open_boundaries,
        land_boundaries=mesh.land_boundaries,
    )
    out, info = repair_skewed_elements(
        mesh, min_angle_deg=1.0, max_angle_deg=175.0,
    )
    assert info["n_elements_removed"] == 0
    assert info.get("skipped") is True
    assert len(out.open_boundaries) == len(mesh.open_boundaries)
    assert len(out.land_boundaries) == len(mesh.land_boundaries)


def test_repair_skewed_elements_requires_bbox_when_deletes() -> None:
    """If a deletion happens, bbox/tol_deg must have been supplied."""
    mesh = _mesh_with_one_sliver()
    import pytest

    with pytest.raises(ValueError, match="Phase F removed elements"):
        repair_skewed_elements(
            mesh, min_angle_deg=1.0, max_angle_deg=175.0,
        )


def test_repair_skewed_elements_rejects_inverted_thresholds() -> None:
    mesh = _mesh_with_one_sliver()
    import pytest

    with pytest.raises(ValueError, match="min_angle_deg"):
        repair_skewed_elements(mesh, min_angle_deg=10.0, max_angle_deg=5.0)


def test_clean_mesh_phase_f_default_off() -> None:
    """Phase F is off by default; the sliver triangle survives."""
    mesh = _mesh_with_one_sliver()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 2.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
    )
    assert all(p["name"] != "repair_skewed_elements" for p in info["phases"])
    assert info["output"]["n_elements"] == mesh.n_elements


def test_clean_mesh_phase_f_explicit_removes_sliver() -> None:
    mesh = _mesh_with_one_sliver()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 2.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        repair_skewed=True,
    )
    phase_f = next(p for p in info["phases"]
                   if p["name"] == "repair_skewed_elements")
    assert phase_f["n_elements_removed"] == 1
    assert info["output"]["n_elements"] == mesh.n_elements - 1


# ---------------------------------------------------------------------------
# Phase G: Laplacian smoothing (wraps oceanmesh.laplacian2)
# ---------------------------------------------------------------------------


def _square_with_offcenter_node() -> Fort14Mesh:
    """Unit square with a centre node *displaced* from the geometric
    centre. Laplacian smoothing should pull it back toward (0.5, 0.5).
    """
    nodes = [
        [0.0, 0.0],          # 0
        [1.0, 0.0],          # 1
        [1.0, 1.0],          # 2
        [0.0, 1.0],          # 3
        [0.7, 0.7],          # 4 — interior, off-centre
    ]
    elements = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
    return _mesh(
        nodes, elements,
        open_boundaries=[np.array([0, 1])],
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )


def test_smooth_mesh_laplacian_pulls_interior_node_to_centre() -> None:
    """The off-centre interior node should move toward (0.5, 0.5);
    the four boundary corners must not move at all.
    """
    mesh = _square_with_offcenter_node()
    pre_centre = mesh.nodes[4].copy()
    out, info = smooth_mesh_laplacian(mesh, max_iter=20, tol=1e-6)

    # Interior moves toward the geometric centre.
    new_centre = out.nodes[4]
    assert np.linalg.norm(new_centre - np.array([0.5, 0.5])) < (
        np.linalg.norm(pre_centre - np.array([0.5, 0.5]))
    )
    # Boundary corners pinned.
    np.testing.assert_array_equal(out.nodes[:4], mesh.nodes[:4])
    # Connectivity, depths, boundaries preserved.
    np.testing.assert_array_equal(out.elements, mesh.elements)
    np.testing.assert_array_equal(out.depths, mesh.depths)
    assert len(out.open_boundaries) == 1
    assert len(out.land_boundaries) == 1
    assert info["max_iter"] == 20
    assert info["n_nodes_moved"] >= 1


def test_smooth_mesh_laplacian_empty_mesh_is_noop() -> None:
    mesh = Fort14Mesh(
        title="empty",
        nodes=np.empty((0, 2), dtype=float),
        depths=np.empty((0,), dtype=float),
        elements=np.empty((0, 3), dtype=np.int64),
        open_boundaries=[], land_boundaries=[],
    )
    out, info = smooth_mesh_laplacian(mesh)
    assert info.get("skipped") is True
    assert out.n_elements == 0


def test_smooth_mesh_laplacian_rejects_invalid_args() -> None:
    mesh = _square_with_offcenter_node()
    import pytest

    with pytest.raises(ValueError, match="max_iter"):
        smooth_mesh_laplacian(mesh, max_iter=0)
    with pytest.raises(ValueError, match="tol"):
        smooth_mesh_laplacian(mesh, tol=0.0)


def test_clean_mesh_phase_g_default_off() -> None:
    mesh = _square_with_offcenter_node()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 1.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
    )
    assert all(p["name"] != "smooth_mesh_laplacian" for p in info["phases"])
    np.testing.assert_array_equal(cleaned.nodes[4], mesh.nodes[4])


def test_clean_mesh_phase_g_explicit_smooths_interior() -> None:
    mesh = _square_with_offcenter_node()
    cleaned, info = clean_mesh(
        mesh,
        bbox=(0.0, 0.0, 1.0, 1.0), bbox_tol_m=1.0,
        remove_disjoint=False, trim_dead_ends_iters=0,
        thin_chain_mode="none",
        smooth_laplacian=True, smooth_laplacian_iters=20,
        smooth_laplacian_tol=1e-6,
    )
    phase_g = next(p for p in info["phases"]
                   if p["name"] == "smooth_mesh_laplacian")
    assert phase_g["n_nodes_moved"] >= 1
    # Interior node moved toward (0.5, 0.5).
    moved = cleaned.nodes[4]
    assert np.linalg.norm(moved - np.array([0.5, 0.5])) < (
        np.linalg.norm(mesh.nodes[4] - np.array([0.5, 0.5]))
    )


# ---------------------------------------------------------------------------
# Phase G repair: flip detection + rollback safety net
# ---------------------------------------------------------------------------


def test_repair_flipped_elements_no_flip_is_noop() -> None:
    """If post is identical to pre (no flips), repair returns post
    unchanged with zeroed counts."""
    from fvcom_mesh_tools.mesh_clean import _repair_flipped_elements

    pre = np.array([[0, 0], [1, 0], [0, 1]], dtype=float)
    post = pre.copy()
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    repaired, info = _repair_flipped_elements(pre, post, elements)
    np.testing.assert_array_equal(repaired, pre)
    assert info["n_flipped_post_smooth"] == 0
    assert info["n_flipped_after_repair"] == 0
    assert info["n_nodes_rolled_back"] == 0
    assert info["n_rollback_passes"] == 0
    assert info["full_rollback"] is False


def test_repair_flipped_elements_rolls_back_offending_nodes() -> None:
    """Construct a hand-crafted post-smoothing array where the centre
    node has been displaced outside the square, flipping every
    triangle. The repair must restore it."""
    from fvcom_mesh_tools.mesh_clean import _repair_flipped_elements

    pre = np.array([
        [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5],
    ], dtype=float)
    post = pre.copy()
    post[4] = [2.0, 2.0]   # outside the square — flips every triangle
    elements = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
                        dtype=np.int64)

    repaired, info = _repair_flipped_elements(pre, post, elements)
    # Repair restored the centre node.
    np.testing.assert_array_equal(repaired[4], pre[4])
    assert info["n_flipped_post_smooth"] >= 1
    assert info["n_flipped_after_repair"] == 0
    assert info["n_nodes_rolled_back"] >= 1
    assert info["full_rollback"] is False


def test_repair_flipped_elements_full_rollback_safety_net() -> None:
    """If max_passes=0 leaves flips in place, the safety net fully
    reverts to ``pre``."""
    from fvcom_mesh_tools.mesh_clean import _repair_flipped_elements

    pre = np.array([
        [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5],
    ], dtype=float)
    post = pre.copy()
    post[4] = [2.0, 2.0]
    elements = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
                        dtype=np.int64)

    repaired, info = _repair_flipped_elements(
        pre, post, elements, max_passes=0,
    )
    # Full rollback fired.
    np.testing.assert_array_equal(repaired, pre)
    assert info["full_rollback"] is True
    assert info["n_flipped_after_repair"] == 0


def test_smooth_mesh_laplacian_emits_repair_info_keys() -> None:
    """Even on a non-flipping fixture, the new info keys are present."""
    mesh = _square_with_offcenter_node()
    _, info = smooth_mesh_laplacian(mesh, max_iter=20, tol=1e-6)
    assert "n_flipped_post_smooth" in info
    assert "n_flipped_after_repair" in info
    assert "n_nodes_rolled_back" in info
    assert "n_rollback_passes" in info
    assert "full_rollback" in info
    # Square test fixture must not flip.
    assert info["n_flipped_post_smooth"] == 0
    assert info["n_flipped_after_repair"] == 0


def test_smooth_mesh_laplacian_repair_off_keeps_flips() -> None:
    """Synthetic test: monkey-patch oceanmesh.laplacian2 to return
    flipped vertices, and verify ``repair_flipped=False`` surfaces
    the raw output (n_flipped_after_repair > 0)."""
    import fvcom_mesh_tools.mesh_clean as mc

    mesh = _square_with_offcenter_node()
    flipped_post = mesh.nodes.copy()
    flipped_post[4] = [2.0, 2.0]

    # Stand-in laplacian2 that ignores its inputs and returns the flip.
    def fake_laplacian2(vertices, entities, **kwargs):  # noqa: ARG001
        return flipped_post.copy(), entities

    import oceanmesh
    real = oceanmesh.laplacian2
    oceanmesh.laplacian2 = fake_laplacian2
    try:
        out_unsafe, info_unsafe = mc.smooth_mesh_laplacian(
            mesh, max_iter=1, tol=1.0, repair_flipped=False,
        )
        out_safe, info_safe = mc.smooth_mesh_laplacian(
            mesh, max_iter=1, tol=1.0, repair_flipped=True,
        )
    finally:
        oceanmesh.laplacian2 = real

    # repair_flipped=False propagates the flips.
    assert info_unsafe["n_flipped_post_smooth"] >= 1
    assert info_unsafe["n_flipped_after_repair"] == info_unsafe[
        "n_flipped_post_smooth"
    ]
    np.testing.assert_array_equal(out_unsafe.nodes[4], flipped_post[4])

    # repair_flipped=True (default) repairs the flips.
    assert info_safe["n_flipped_post_smooth"] >= 1
    assert info_safe["n_flipped_after_repair"] == 0
    # The centre node was reverted.
    np.testing.assert_array_equal(out_safe.nodes[4], mesh.nodes[4])


def test_smooth_mesh_laplacian_rejects_negative_max_repair_passes() -> None:
    mesh = _square_with_offcenter_node()
    import pytest

    with pytest.raises(ValueError, match="max_repair_passes"):
        smooth_mesh_laplacian(mesh, max_repair_passes=-1)


# ---------------------------------------------------------------------------
# Phase E Stage 1: medial-axis potential analysis (no re-meshing)
# ---------------------------------------------------------------------------


def test_analyze_under_resolved_channels_clean_mesh_returns_zero() -> None:
    """A wide 6-row strip has no flagged elements at min_w_h=3, so the
    analysis returns an empty report."""
    mesh = _strip_mesh_n_rows(
        length_deg=0.5, width_deg=0.10, n_x=20, n_rows=6,
    )
    report = analyze_under_resolved_channels(mesh, min_w_h=0.1)
    assert report["total_flagged_elements"] == 0
    assert report["n_components"] == 0
    assert report["components"] == []
    assert report["current_phase_e_new_nodes"] == 0
    assert report["medial_axis_new_nodes_estimate"] == 0


def test_analyze_under_resolved_channels_one_row_strip_one_component() -> None:
    """A 1-row strip is fully flagged at min_w_h=3 and forms a single
    face-face-adjacent component."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    report = analyze_under_resolved_channels(
        mesh, min_w_h=3.0, target_cells_across=3,
    )
    assert report["total_flagged_elements"] == mesh.n_elements
    assert report["n_components"] == 1
    comp = report["components"][0]
    assert comp["n_elements"] == mesh.n_elements
    # All 18 nodes (2 rows x 9 columns) belong to the strip.
    assert comp["n_nodes"] == mesh.n_nodes
    # h_local_median should be in metres (lon/lat strip ≈ 1 km in 0.01°).
    assert 100.0 < comp["h_local_median_m"] < 1.5e4
    assert comp["long_axis_m"] > comp["h_local_median_m"]
    # Current Phase E cost = one centroid per flagged element.
    assert comp["current_phase_e_new_nodes"] == mesh.n_elements


def test_analyze_under_resolved_channels_aggregates() -> None:
    """The aggregate dict equals the sum of per-component contributions."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    report = analyze_under_resolved_channels(mesh, min_w_h=3.0)
    assert report["current_phase_e_new_nodes"] == sum(
        c["current_phase_e_new_nodes"] for c in report["components"]
    )
    assert report["medial_axis_new_nodes_estimate"] == sum(
        c["medial_axis_new_nodes_estimate"] for c in report["components"]
    )
    assert report["delta_nodes_vs_current"] == (
        report["medial_axis_new_nodes_estimate"]
        - report["current_phase_e_new_nodes"]
    )


def test_analyze_under_resolved_channels_target_cells_must_be_at_least_2() -> None:
    """The medial-axis interior rows = target_cells_across - 1 must be
    at least 1."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    import pytest

    with pytest.raises(ValueError, match="target_cells_across"):
        analyze_under_resolved_channels(mesh, target_cells_across=1)


def test_analyze_under_resolved_channels_target_4_cells_doubles_estimate() -> None:
    """target_cells_across=4 gives 3 interior rows; the estimate should
    be 3/2 = 1.5× the target_cells_across=3 case (target_cells - 1
    rows, same nodes_per_row)."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    r3 = analyze_under_resolved_channels(mesh, min_w_h=3.0, target_cells_across=3)
    r4 = analyze_under_resolved_channels(mesh, min_w_h=3.0, target_cells_across=4)
    # target_cells - 1 rows of nodes_per_row each.
    assert r3["medial_axis_new_nodes_estimate"] > 0
    assert r4["medial_axis_new_nodes_estimate"] > r3["medial_axis_new_nodes_estimate"]
    # The ratio should be 3/2 (3 rows / 2 rows).
    np.testing.assert_allclose(
        r4["medial_axis_new_nodes_estimate"]
        / r3["medial_axis_new_nodes_estimate"],
        3.0 / 2.0, atol=0.5,  # rounding in nodes-per-row
    )


def test_analyze_under_resolved_channels_min_channel_elements_filters() -> None:
    """``min_channel_elements`` larger than every component's size drops
    all flagged elements before the per-component tally is computed."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    base = analyze_under_resolved_channels(mesh, min_w_h=3.0)
    assert base["total_flagged_elements"] > 0
    assert base["n_components"] >= 1

    filtered = analyze_under_resolved_channels(
        mesh, min_w_h=3.0, min_channel_elements=base["total_flagged_elements"] + 1,
    )
    assert filtered["min_channel_elements"] == base["total_flagged_elements"] + 1
    assert filtered["total_flagged_elements"] == 0
    assert filtered["n_components"] == 0
    assert filtered["components"] == []
    assert filtered["current_phase_e_new_nodes"] == 0
    assert filtered["medial_axis_new_nodes_estimate"] == 0


# ---------------------------------------------------------------------------
# compact_nodes
# ---------------------------------------------------------------------------


def test_compact_nodes_drops_orphans_and_remaps_boundaries():
    from fvcom_mesh_tools.mesh_clean import compact_nodes

    # Node 2 is an orphan in the middle of the index range; elements
    # and boundaries reference nodes on either side of it.
    nodes = np.array([
        [0.0, 0.0], [1000.0, 0.0], [9999.0, 9999.0],
        [1000.0, 1000.0], [0.0, 1000.0],
    ])
    mesh = Fort14Mesh(
        title="orphan",
        nodes=nodes,
        depths=np.array([2.0, 3.0, 99.0, 4.0, 5.0]),
        elements=np.array([[0, 1, 3], [0, 3, 4]]),
        open_boundaries=[np.array([1, 3])],
        land_boundaries=[(20, np.array([3, 4, 0, 1]))],
    )
    out, info = compact_nodes(mesh)
    assert info["n_orphans_removed"] == 1
    assert info["n_boundary_refs_dropped"] == 0
    assert out.n_nodes == 4
    # Old node 3 -> new node 2, old 4 -> 3; coordinates/depths follow.
    assert np.allclose(out.nodes[2], [1000.0, 1000.0])
    assert np.allclose(out.depths, [2.0, 3.0, 4.0, 5.0])
    assert out.elements.max() == 3
    assert list(out.open_boundaries[0]) == [1, 2]
    assert list(out.land_boundaries[0][1]) == [2, 3, 0, 1]


def test_compact_nodes_noop_when_dense():
    from fvcom_mesh_tools.mesh_clean import compact_nodes

    mesh = Fort14Mesh(
        title="dense",
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        depths=np.full(3, 2.0),
        elements=np.array([[0, 1, 2]]),
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2]))],
    )
    out, info = compact_nodes(mesh)
    assert info["n_orphans_removed"] == 0
    assert out is mesh


def test_weld_close_nodes_merges_and_drops_degenerates():
    from fvcom_mesh_tools.mesh_clean import weld_close_nodes

    # Nodes 1 and 2 nearly coincide; the triangle using both becomes
    # degenerate and must vanish; boundaries remap without repeats.
    nodes = np.array([
        [0.0, 0.0], [1000.0, 0.0], [1000.0, 0.4], [2000.0, 0.0],
        [1000.0, 1000.0],
    ])
    mesh = Fort14Mesh(
        title="weld",
        nodes=nodes,
        depths=np.array([2.0, 3.0, 4.0, 5.0, 6.0]),
        elements=np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4]]),
        open_boundaries=[np.array([0, 1, 2, 3])],
        land_boundaries=[(20, np.array([3, 4, 0]))],
    )
    out, info = weld_close_nodes(mesh, tol=1.0)
    assert info["n_welded"] == 1
    assert info["n_elements_dropped"] == 1
    assert out.n_elements == 2
    assert out.n_nodes == 4
    seg = out.open_boundaries[0]
    assert len(seg) == 3 and len(np.unique(seg)) == 3
    from fvcom_mesh_tools.algorithms import signed_areas

    assert (signed_areas(out) > 0).all()
