"""Tests for ``fvcom_mesh_tools.mesh_clean``."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
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
