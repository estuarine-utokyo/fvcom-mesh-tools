"""Tests for ``fvcom_mesh_tools.diagnostics``.

Each test builds a small hand-crafted Fort14Mesh that exercises one
detector in isolation; the integration test at the bottom drives
:func:`run_diagnostics` end-to-end on a clean mesh.
"""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.diagnostics import (
    channel_width_metric,
    dead_end_elements_flag,
    disjoint_components_flag,
    face_face_adjacency,
    node_valence,
    overconnected_nodes_flag,
    run_diagnostics,
    thin_chain_elements_flag,
    thin_elements_flag,
    under_resolved_channels_flag,
    unreachable_elements_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _mesh(
    nodes: list[list[float]] | np.ndarray,
    elements: list[list[int]] | np.ndarray,
    *,
    open_boundaries: list[np.ndarray] | None = None,
    land_boundaries: list[tuple[int, np.ndarray]] | None = None,
) -> Fort14Mesh:
    nodes = np.asarray(nodes, dtype=np.float64)
    elements = np.asarray(elements, dtype=np.int64)
    return Fort14Mesh(
        title="test",
        nodes=nodes,
        depths=np.zeros(len(nodes), dtype=np.float64),
        elements=elements,
        open_boundaries=open_boundaries or [],
        land_boundaries=land_boundaries or [],
    )


def _square_around_centre() -> Fort14Mesh:
    """Unit square with a central node, 4 triangles. Open boundary on the
    bottom edge, land boundary on the remaining 3 sides — every triangle
    has exactly one interior vertex (node 4) so the thin detector should
    not fire.
    """
    nodes = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
    elements = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
    return _mesh(
        nodes, elements,
        open_boundaries=[np.array([0, 1])],
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )


# ---------------------------------------------------------------------------
# Mesh primitives
# ---------------------------------------------------------------------------


def test_face_face_adjacency_counts_neighbours_correctly() -> None:
    mesh = _square_around_centre()
    adj = face_face_adjacency(mesh.elements)
    # The four "wedge" triangles each share one edge with each adjacent
    # wedge (cyclic), so degree=2 for every element.
    deg = np.asarray(adj.sum(axis=1)).ravel()
    np.testing.assert_array_equal(deg, [2, 2, 2, 2])


def test_node_valence_counts_incident_triangles() -> None:
    mesh = _square_around_centre()
    val = node_valence(mesh.elements, mesh.n_nodes)
    # Outer 4 nodes appear in 2 triangles each; the centre node appears in 4.
    np.testing.assert_array_equal(val, [2, 2, 2, 2, 4])


# ---------------------------------------------------------------------------
# Detector: disjoint components
# ---------------------------------------------------------------------------


def test_disjoint_components_flags_smaller_disjoint_subgraph() -> None:
    """Two independent triangles — one larger group, one solitary triangle.
    The smaller group must be flagged.
    """
    base = _square_around_centre()  # 4 elements
    # Add an isolated triangle with a fresh set of 3 nodes far away.
    extra_nodes = np.array([[10.0, 0.0], [11.0, 0.0], [10.5, 1.0]])
    nodes = np.vstack([base.nodes, extra_nodes])
    extra_elements = np.array([[5, 6, 7]])
    elements = np.vstack([base.elements, extra_elements])
    mesh = _mesh(
        nodes, elements,
        open_boundaries=base.open_boundaries,
        land_boundaries=base.land_boundaries,
    )
    adj = face_face_adjacency(mesh.elements)
    flag, labels, sizes = disjoint_components_flag(adj)
    # Expect exactly one flagged element (the lone extra triangle).
    assert flag.sum() == 1
    assert flag[4]            # the isolated triangle is at index 4
    assert not flag[:4].any()
    assert int(sizes.max()) == 4


def test_disjoint_single_component_returns_no_flag() -> None:
    mesh = _square_around_centre()
    adj = face_face_adjacency(mesh.elements)
    flag, _labels, sizes = disjoint_components_flag(adj)
    assert flag.sum() == 0
    np.testing.assert_array_equal(sizes, [4])


# ---------------------------------------------------------------------------
# Detector: thin elements / thin chains
# ---------------------------------------------------------------------------


def _thin_chain_mesh(chain_length: int) -> Fort14Mesh:
    """Strip of ``chain_length`` triangles between two parallel lines, so
    every node lies on a (land) boundary. The whole strip is therefore a
    chain of thin elements of the requested length.
    """
    n = chain_length + 1     # nodes per row
    top = np.column_stack([np.arange(n, dtype=np.float64), np.ones(n)])
    bot = np.column_stack([np.arange(n, dtype=np.float64), np.zeros(n)])
    nodes = np.vstack([bot, top])
    elements = []
    for i in range(chain_length):
        # i_bot, i+1_bot, i_top  + i+1_bot, i+1_top, i_top
        elements.append([i, i + 1, n + i])
        elements.append([i + 1, n + i + 1, n + i])
    elements = np.asarray(elements, dtype=np.int64)
    bot_chain = np.arange(n, dtype=np.int64)
    top_chain = np.arange(n, dtype=np.int64) + n
    return _mesh(
        nodes, elements,
        land_boundaries=[(0, bot_chain), (0, top_chain)],
    )


def test_thin_elements_fire_on_strip_with_all_boundary_nodes() -> None:
    mesh = _thin_chain_mesh(chain_length=4)
    flag = thin_elements_flag(mesh)
    # Every element in this strip has 3 boundary vertices.
    assert flag.all()


def test_thin_chain_filters_short_chains_below_threshold() -> None:
    """A chain of 8 thin elements with min_chain_length=10 should not fire."""
    mesh = _thin_chain_mesh(chain_length=4)   # 8 thin elements
    adj = face_face_adjacency(mesh.elements)
    thin = thin_elements_flag(mesh)
    flag = thin_chain_elements_flag(adj, thin, min_chain_length=10)
    assert flag.sum() == 0


def test_thin_chain_fires_when_chain_length_meets_threshold() -> None:
    mesh = _thin_chain_mesh(chain_length=5)   # 10 thin elements
    adj = face_face_adjacency(mesh.elements)
    thin = thin_elements_flag(mesh)
    flag = thin_chain_elements_flag(adj, thin, min_chain_length=10)
    assert flag.sum() == 10


def test_thin_chain_does_not_fire_on_clean_mesh() -> None:
    mesh = _square_around_centre()
    adj = face_face_adjacency(mesh.elements)
    thin = thin_elements_flag(mesh)
    flag = thin_chain_elements_flag(adj, thin, min_chain_length=3)
    assert flag.sum() == 0


# ---------------------------------------------------------------------------
# Detector: over-connected nodes
# ---------------------------------------------------------------------------


def _fan_mesh(n_wedges: int) -> Fort14Mesh:
    """Regular fan of ``n_wedges`` triangles around a central node. The
    centre node has valence ``n_wedges``; outer nodes have valence 2.
    """
    centre = np.array([[0.0, 0.0]])
    angles = np.linspace(0, 2 * np.pi, n_wedges + 1)[:-1]
    rim = np.column_stack([np.cos(angles), np.sin(angles)])
    nodes = np.vstack([centre, rim])
    elements = np.array(
        [[0, 1 + i, 1 + (i + 1) % n_wedges] for i in range(n_wedges)],
        dtype=np.int64,
    )
    rim_chain = np.concatenate([np.arange(1, n_wedges + 1), [1]])
    return _mesh(nodes, elements, land_boundaries=[(0, rim_chain)])


def test_overconnected_node_detected_above_threshold() -> None:
    mesh = _fan_mesh(n_wedges=12)
    flag, val = overconnected_nodes_flag(mesh.elements, mesh.n_nodes, max_nbr=8)
    assert flag[0]              # centre node (valence 12) is over-connected
    assert not flag[1:].any()   # rim nodes have valence 2
    assert int(val.max()) == 12


def test_overconnected_node_quiet_below_threshold() -> None:
    mesh = _fan_mesh(n_wedges=6)
    flag, val = overconnected_nodes_flag(mesh.elements, mesh.n_nodes, max_nbr=8)
    assert not flag.any()
    assert int(val.max()) == 6


# ---------------------------------------------------------------------------
# Detector: dead-end elements
# ---------------------------------------------------------------------------


def test_dead_end_excludes_open_boundary_corner() -> None:
    """Wedge mesh with the entire rim declared as the OPEN boundary. Every
    triangle is degree-2 in the dual; a corner triangle in particular is
    degree-2 and shares an OB edge, so the dead-end detector must stay
    silent.
    """
    n = 6
    centre = np.array([[0.0, 0.0]])
    angles = np.linspace(0, 2 * np.pi, n + 1)[:-1]
    rim = np.column_stack([np.cos(angles), np.sin(angles)])
    nodes = np.vstack([centre, rim])
    elements = np.array(
        [[0, 1 + i, 1 + (i + 1) % n] for i in range(n)], dtype=np.int64,
    )
    rim_chain = np.concatenate([np.arange(1, n + 1), [1]])
    mesh = _mesh(
        nodes, elements,
        open_boundaries=[rim_chain],
    )
    adj = face_face_adjacency(mesh.elements)
    flag = dead_end_elements_flag(adj, mesh)
    assert not flag.any()


def test_dead_end_fires_on_dangling_triangle() -> None:
    """A dangling triangle attached by a single edge to the rest, where
    both of its other edges are land boundary (no open boundary) — the
    dead-end detector must fire.
    """
    nodes = np.array([
        [0, 0], [1, 0], [0, 1],     # main triangle 0,1,2
        [1, 0], [2, 0], [1.5, 1],   # dangling triangle attached via 1-2 edge
    ], dtype=np.float64)
    # Re-use shared edge: triangle B = (1, 4, 5) shares edge (1, ?) — let's
    # instead set up a clean configuration: two triangles sharing edge 1-2,
    # plus a third triangle dangling off via edge 4-5 with one shared node.
    nodes = np.array([
        [0, 0],     # 0
        [1, 0],     # 1
        [0.5, 1],   # 2
        [2, 0],     # 3 (extends right of 1)
        [1.5, 1],   # 4
    ], dtype=np.float64)
    elements = np.array([
        [0, 1, 2],       # main triangle
        [1, 3, 4],       # neighbour, shares only node 1 — wait that's not adjacent
    ], dtype=np.int64)
    # The above does NOT share an edge. We need to share an edge for the
    # dual graph to connect them at all. Build properly:
    nodes = np.array([
        [0, 0],   # 0
        [1, 0],   # 1
        [0.5, 1], # 2
        [1.5, 1], # 3
    ], dtype=np.float64)
    elements = np.array([
        [0, 1, 2],   # element 0
        [1, 3, 2],   # element 1, shares edge 1-2 with element 0
    ], dtype=np.int64)
    # Element 0 has neighbours: only element 1 (degree=1)
    # Element 1 has neighbours: only element 0 (degree=1)
    # No open boundary -> both flagged as dead-end.
    mesh = _mesh(
        nodes, elements,
        land_boundaries=[(0, np.array([0, 1, 3, 2, 0]))],
    )
    adj = face_face_adjacency(mesh.elements)
    flag = dead_end_elements_flag(adj, mesh)
    assert flag.all()


# ---------------------------------------------------------------------------
# Detector: unreachable
# ---------------------------------------------------------------------------


def test_unreachable_silent_when_no_open_boundary() -> None:
    mesh = _square_around_centre()
    mesh.open_boundaries.clear()  # no open boundary
    adj = face_face_adjacency(mesh.elements)
    flag = unreachable_elements_flag(adj, mesh)
    assert not flag.any()


def test_unreachable_fires_on_disjoint_component_without_ob() -> None:
    """Two disjoint components, the second contains no open-boundary node.
    Detector 5 must flag the elements of the second component.
    """
    base = _square_around_centre()
    extra_nodes = np.array([[10.0, 0.0], [11.0, 0.0], [10.5, 1.0]])
    nodes = np.vstack([base.nodes, extra_nodes])
    extra_elements = np.array([[5, 6, 7]])
    elements = np.vstack([base.elements, extra_elements])
    mesh = _mesh(
        nodes, elements,
        open_boundaries=[np.array([0, 1])],            # only on the first component
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )
    adj = face_face_adjacency(mesh.elements)
    flag = unreachable_elements_flag(adj, mesh)
    assert flag[4]                  # the isolated triangle is unreachable
    assert not flag[:4].any()       # the OB-touching component is reachable


def test_unreachable_quiet_when_disjoint_component_touches_ob() -> None:
    """Same disjoint pair, but the second component carries its own open
    boundary segment. Both components are reachable from an OB node, so
    detector 5 stays silent (detector 1 still flags the smaller one).
    """
    base = _square_around_centre()
    extra_nodes = np.array([[10.0, 0.0], [11.0, 0.0], [10.5, 1.0]])
    nodes = np.vstack([base.nodes, extra_nodes])
    extra_elements = np.array([[5, 6, 7]])
    elements = np.vstack([base.elements, extra_elements])
    mesh = _mesh(
        nodes, elements,
        open_boundaries=[
            np.array([0, 1]),            # on base
            np.array([5, 6]),            # on extra triangle
        ],
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )
    adj = face_face_adjacency(mesh.elements)
    flag = unreachable_elements_flag(adj, mesh)
    assert not flag.any()


# ---------------------------------------------------------------------------
# Integration: run_diagnostics on a clean mesh
# ---------------------------------------------------------------------------


def test_run_diagnostics_clean_mesh_reports_no_flags() -> None:
    """The 4-element square mesh is too coarse for the medial-axis
    detector (its single 'channel' is one cell across), so we disable
    detector 7 by setting min_w_h=0 in this test of the other six.
    """
    mesh = _square_around_centre()
    report = run_diagnostics(mesh, name="clean", min_w_h=0.0)
    assert not report.any_flagged()
    assert report.n_nodes == 5
    assert report.n_elements == 4
    assert report.max_nbr_elem == 8
    assert int(report.valence.max()) == 4


def test_run_diagnostics_with_overconnected_threshold() -> None:
    """A 12-wedge fan mesh — over-connected with default cap, clean if cap
    is raised to 12.
    """
    mesh = _fan_mesh(n_wedges=12)
    rep_default = run_diagnostics(mesh, max_nbr_elem=8)
    rep_relaxed = run_diagnostics(mesh, max_nbr_elem=12)
    assert rep_default.overconnected_flag.any()
    assert not rep_relaxed.overconnected_flag.any()


# ---------------------------------------------------------------------------
# Detector 7: channel-width / h ratio (medial-axis-style)
# ---------------------------------------------------------------------------


def _strip_mesh_n_rows(
    *, length_deg: float, width_deg: float, n_x: int, n_rows: int,
) -> Fort14Mesh:
    """A rectangular ``length_deg × width_deg`` strip in lon/lat space,
    triangulated as ``n_rows`` rows × ``n_x`` columns of quads (each quad
    split into two triangles). Top and bottom long edges are land; the
    two short edges are open.

    With more rows the cross-channel resolution increases and the local
    h shrinks, so the channel-width / h ratio scales roughly linearly
    with ``n_rows``.
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


def test_channel_width_metric_one_row_strip_flagged() -> None:
    """A 1-row strip is one cell across the channel; w/h ≈ 1 < 3."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    flag, info = under_resolved_channels_flag(mesh, min_w_h=3.0)
    assert flag.all()
    assert (info["w_h_ratio"] < 3.0).all()


def test_channel_width_metric_six_row_strip_partially_unflagged() -> None:
    """A 6-row strip has ~6 cells across the channel; w/h is well above
    3 in the interior even though boundary-adjacent rows can be lower.
    """
    mesh = _strip_mesh_n_rows(
        length_deg=0.5, width_deg=0.10, n_x=20, n_rows=6,
    )
    flag, info = under_resolved_channels_flag(mesh, min_w_h=3.0)
    # The middle-row elements are away from both banks, so most of them
    # should NOT be flagged.
    assert (~flag).any()
    # Median ratio should be comfortably above 3.
    assert float(np.median(info["w_h_ratio"])) > 3.0


def test_channel_width_metric_handles_no_boundary() -> None:
    """A mesh with no boundary lists returns inf-ratio (no flags)."""
    mesh = _square_around_centre()
    mesh.open_boundaries.clear()
    mesh.land_boundaries.clear()
    flag, info = under_resolved_channels_flag(mesh, min_w_h=3.0)
    assert not flag.any()
    assert np.isinf(info["w_h_ratio"]).all()


def test_run_diagnostics_includes_under_resolved_channels() -> None:
    """The 7th detector is wired into the high-level driver and the
    summary text mentions the channel ratio."""
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    rep = run_diagnostics(mesh, min_w_h=3.0)
    assert rep.under_resolved_channels_flag.any()
    assert rep.min_w_h == 3.0
    # The ratio array is populated.
    assert rep.w_h_ratio.shape == (mesh.n_elements,)


def test_channel_width_metric_metric_keys_present() -> None:
    mesh = _strip_mesh_n_rows(length_deg=0.08, width_deg=0.01, n_x=8, n_rows=1)
    info = channel_width_metric(mesh)
    for key in ("channel_width_m", "h_local_m", "w_h_ratio", "sample_count"):
        assert key in info
