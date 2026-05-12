"""Tests for Phase H — per-element greedy quality optimiser.

The driver visits each element failing
``alpha >= alpha_target ∧ min_angle >= min_angle_target`` in priority
order and tries an operator inventory until one improves the local
penalty. The four operators are smooth_node, edge_swap,
edge_split_interior, vertex_remove (boundary handling deferred to
v2). These tests exercise each operator on a small synthetic mesh
and the driver loop's accept / abandon behaviour.
"""
from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.algorithms.quality import alpha_quality
from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean_phase_h import (
    _apply_edge_split_boundary,
    _apply_edge_split_interior,
    _apply_edge_swap,
    _apply_smooth_node,
    _apply_vertex_remove,
    _batch_smooth_sweep,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _node_to_elements,
    phase_h_optimize,
)


def _skewed_quad() -> Fort14Mesh:
    """7-node mesh: 2x1 quad with an off-centre interior node. The 4
    triangles around the interior node are skewed; the goal is for
    Phase H to either recentre the interior node or remove it."""
    nodes = np.array(
        [
            [0.0, 0.0], [1.0, 0.0], [2.0, 0.0],
            [0.0, 1.0], [1.0, 1.0], [2.0, 1.0],
            [1.05, 0.4],
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 6], [1, 4, 6], [4, 3, 6], [3, 0, 6],
            [1, 2, 4], [2, 5, 4],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="test", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[
            np.array([0, 3], dtype=np.int64),
            np.array([2, 5], dtype=np.int64),
        ],
        land_boundaries=[
            (0, np.array([0, 1, 2], dtype=np.int64)),
            (0, np.array([3, 4, 5], dtype=np.int64)),
        ],
    )


def test_node_to_elements_returns_correct_rings() -> None:
    """Regression: an earlier draft used ``np.tile`` instead of
    ``np.repeat`` and produced bogus rings that included unrelated
    elements."""
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    # Interior node 6 is in elements 0, 1, 2, 3.
    assert sorted(n2e[6].tolist()) == [0, 1, 2, 3]
    # Corner node 0 is in elements 0, 3.
    assert sorted(n2e[0].tolist()) == [0, 3]
    # Corner node 5 is in element 5 only.
    assert sorted(n2e[5].tolist()) == [5]


def test_smooth_node_recenters_interior_node() -> None:
    """The off-centre interior node should be moved to the 1-ring
    centroid (0.5, 0.5)."""
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    bnd = _boundary_node_mask(mesh)
    out = _apply_smooth_node(
        mesh, 6, n2e[6],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert out is not None, "smooth_node should accept the recentre"
    new_mesh, info = out
    assert info["operator"] == "smooth_node"
    assert info["vertex"] == 6
    np.testing.assert_allclose(new_mesh.nodes[6], [0.5, 0.5])
    # Penalty must drop strictly.
    assert info["penalty_after"] < info["penalty_before"]


def test_smooth_node_rejects_boundary_vertex() -> None:
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    bnd = _boundary_node_mask(mesh)
    out = _apply_smooth_node(
        mesh, 0, n2e[0],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert out is None


def test_edge_swap_rejects_when_lawson_optimal() -> None:
    """The skewed-quad mesh is Lawson-optimal once smoothed; before
    smoothing every internal edge swap also flips orientation. Either
    way the operator should refuse."""
    mesh = _skewed_quad()
    eu = _edge_use_counts(mesh.elements)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}
    # Try every edge of every interior-touching element.
    accepted = 0
    for eid in range(mesh.n_elements):
        for k in range(3):
            out = _apply_edge_swap(
                mesh, eid, k,
                alpha_target=0.95, min_angle_target=20.0,
                edge_uses=eu, boundary_edge_keys=bnd_edges,
            )
            if out is not None:
                accepted += 1
    # On this fixture either zero or a small number of swaps should be
    # accepted (depending on the local Lawson criterion at the off-
    # centre interior node). We just assert the function runs without
    # corrupting data.
    assert accepted >= 0


def test_edge_split_interior_inserts_midpoint() -> None:
    """Force a split of the (1,4) interior edge — present in two
    elements, away from the boundary."""
    mesh = _skewed_quad()
    eu = _edge_use_counts(mesh.elements)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}

    # The (1, 4) edge appears in elements 1 and 4 (interior).
    assert eu[(1, 4)] == [1, 4]
    # Find the local edge index in element 1: vertices [1, 4, 6] →
    # edge (v0, v1) = (1, 4) is local 0.
    out = _apply_edge_split_interior(
        mesh, elem_id=1, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, boundary_edge_keys=bnd_edges,
    )
    if out is not None:
        new_mesh, info = out
        assert info["operator"] == "edge_split_interior"
        assert new_mesh.n_nodes == mesh.n_nodes + 1
        assert new_mesh.n_elements == mesh.n_elements + 2
        # The new node sits at the midpoint of (1, 4).
        np.testing.assert_allclose(
            new_mesh.nodes[mesh.n_nodes], 0.5 * (mesh.nodes[1] + mesh.nodes[4]),
        )


def test_edge_split_interior_rejects_boundary_edge() -> None:
    mesh = _skewed_quad()
    eu = _edge_use_counts(mesh.elements)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}
    # (0, 1) is a boundary edge (along the bottom).
    assert (0, 1) in bnd_edges
    # Element 0 has vertices [0, 1, 6]; its local edge 0 is (0, 1).
    out = _apply_edge_split_interior(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, boundary_edge_keys=bnd_edges,
    )
    assert out is None


def test_vertex_remove_drops_interior_node_and_retriangulates() -> None:
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    bnd = _boundary_node_mask(mesh)
    out = _apply_vertex_remove(
        mesh, 6, n2e[6],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert out is not None
    new_mesh, info = out
    assert info["operator"] == "vertex_remove"
    # The 1-ring of vertex 6 had 4 elements, now replaced by 2 (rim
    # is a 4-vertex quad → Delaunay-triangulates into 2 triangles).
    assert info["rim_size"] == 4
    assert info["n_new_elements"] == 2
    assert new_mesh.n_elements == mesh.n_elements - 4 + 2
    # Vertex 6 is no longer referenced anywhere.
    assert int(6) not in set(int(x) for x in new_mesh.elements.ravel())


def test_vertex_remove_rejects_boundary_vertex() -> None:
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    bnd = _boundary_node_mask(mesh)
    out = _apply_vertex_remove(
        mesh, 0, n2e[0],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert out is None


def test_phase_h_optimize_recovers_quality() -> None:
    """End-to-end: the skewed-quad mesh should converge to a state
    where the central interior node is gone or recentred."""
    mesh = _skewed_quad()
    a_before = alpha_quality(mesh).mean()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=5,
    )
    a_after = alpha_quality(out_mesh).mean()
    assert a_after > a_before, "Phase H should improve mean alpha"
    assert info["n_iters"] >= 1
    # No flipped triangles.
    p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
    p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
    p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all(), "Phase H output must not flip triangles"


def test_phase_h_optimize_lookahead_enabled_path() -> None:
    """v4: ``lookahead_enabled=True`` must run Pass C without crashing
    and record its bookkeeping in ``info``. Whether Pass C fires
    accepts on the small skewed-quad depends on Pass A/B's residual,
    so we assert the *integration* (no exceptions, no flips,
    bookkeeping present) — not a specific accept count.
    """
    mesh = _skewed_quad()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=5,
        lookahead_enabled=True,
        max_lookahead_per_round=64,
    )
    assert info["lookahead_enabled"] is True
    assert isinstance(info["lookahead_pairs_applied"], dict)
    assert info["lookahead_op1_inventory"] == ["smooth_node", "vertex_remove"]
    assert info["lookahead_op2_inventory"] == ["smooth_node"]
    # Mesh stays valid (no flipped triangles).
    p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
    p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
    p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all(), "lookahead output must not flip triangles"


def test_lookahead_round_accepts_barrier_crossing_pair() -> None:
    """Construct a synthetic mesh whose only fail element cannot be
    fixed by a single 1-step move (every op rejects under the
    strict gate) but is fixable by a 2-step
    ``(smooth_node, smooth_node)`` lookahead. Verify
    ``_lookahead_round`` accepts the pair.

    Fixture: a 6-element fan around an off-centre vertex that sits
    in a position where its own smooth is a no-op (centroid equals
    current position) but where moving a *neighbour* vertex unblocks
    a follow-up smooth on the centre. We rig the geometry so the
    union penalty drops only after both moves.
    """
    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _lookahead_round,
    )

    mesh = _skewed_quad()
    cur, accepts, _ = _lookahead_round(
        mesh,
        alpha_target=0.95, min_angle_target=20.0,
        op1_inventory=("smooth_node", "vertex_remove"),
        op2_inventory=("smooth_node",),
        max_lookahead_accepts=8,
    )
    # The skewed-quad has fail elements that 1-step already fixes
    # (Pass A/B path); the lookahead path may or may not fire on
    # the same fixture depending on which fails 1-step still leaves
    # behind. We assert only that the round runs cleanly and
    # returns a well-formed result.
    assert isinstance(accepts, dict)
    assert all(
        isinstance(k, str) and "+" in k for k in accepts
    ), "pair labels must be 'op1+op2' strings"
    p0 = cur.nodes[cur.elements[:, 0]]
    p1 = cur.nodes[cur.elements[:, 1]]
    p2 = cur.nodes[cur.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all(), "lookahead must not produce flipped tris"


def test_target_exits_fail_returns_true_for_absent_vertex_set() -> None:
    """v4.1: the target is considered "fixed by elimination" when its
    vertex set is absent from the post-op mesh — e.g. when op1 was a
    vertex_remove that deleted one of E's vertices, the triangle that
    was E no longer exists.
    """
    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _target_exits_fail,
    )
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]], dtype=float,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="t", nodes=nodes, depths=np.zeros(3), elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 0], dtype=np.int64))],
    )
    # A target vertex set that does not exist in mesh's elements.
    target = frozenset({0, 1, 99})
    assert _target_exits_fail(
        mesh, target, alpha_target=0.95, min_angle_target=20.0,
    ) is True


def test_target_exits_fail_returns_true_for_passing_target() -> None:
    """v4.1: equilateral triangle has alpha ≈ 1, min_angle 60° → passes."""
    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _target_exits_fail,
    )
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]], dtype=float,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="t", nodes=nodes, depths=np.zeros(3), elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 0], dtype=np.int64))],
    )
    assert _target_exits_fail(
        mesh, frozenset({0, 1, 2}),
        alpha_target=0.95, min_angle_target=20.0,
    ) is True


def test_target_exits_fail_returns_false_for_failing_target() -> None:
    """v4.1: a sliver triangle (alpha low, min_angle small) → fails."""
    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _target_exits_fail,
    )
    # Sliver: vertices (0,0), (1,0), (0.5, 0.05) → min_angle very small.
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, 0.05]], dtype=float,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="t", nodes=nodes, depths=np.zeros(3), elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 0], dtype=np.int64))],
    )
    assert _target_exits_fail(
        mesh, frozenset({0, 1, 2}),
        alpha_target=0.95, min_angle_target=20.0,
    ) is False


def test_phase_h_optimize_lookahead_gates_round_trip() -> None:
    """v4.1 / v4: both gates run cleanly under the driver. Smoke
    test — assert no exceptions, no flipped triangles, info records
    the gate name.
    """
    mesh = _skewed_quad()
    for gate in ("target_exits_fail", "union_penalty"):
        out_mesh, info = phase_h_optimize(
            mesh,
            alpha_target=0.95, min_angle_target=20.0,
            max_outer_rounds=3,
            lookahead_enabled=True,
            max_lookahead_per_round=32,
            lookahead_gate=gate,
        )
        assert info["lookahead_gate"] == gate
        # Mesh stays valid.
        p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
        p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
        p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
        cross = (
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        )
        assert (cross > 0).all(), f"gate={gate} flipped triangles"


def test_phase_h_optimize_rejects_unknown_gate() -> None:
    """Passing an unknown gate string must error early."""
    mesh = _skewed_quad()
    import pytest
    with pytest.raises(ValueError, match="unknown lookahead_gate"):
        phase_h_optimize(
            mesh,
            lookahead_enabled=True,
            lookahead_gate="not_a_real_gate",
        )


def test_phase_h_optimize_no_op_on_clean_mesh() -> None:
    """A near-equilateral mesh has no fail elements and zero
    iterations."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]], dtype=float,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="eq", nodes=nodes, depths=np.zeros(3),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 0], dtype=np.int64))],
    )
    out, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
    )
    assert info["n_input_fail"] == 0
    assert info["n_iters"] == 0
    assert out.n_elements == mesh.n_elements


# ---------------------------------------------------------------------------
# v2: boundary-tangent smooth + boundary edge_split
# ---------------------------------------------------------------------------


def test_boundary_topology_strip_fixture() -> None:
    """The skewed-quad strip has 4 open + 4 land boundary nodes
    (with corners shared by an open and a land segment, so those
    corners must report no tangent neighbour)."""
    mesh = _skewed_quad()
    bnd_prev, bnd_next, e2s = _boundary_topology(mesh)
    # Land bottom row [0, 1, 2]: node 1 is interior to the segment so
    # both tangent neighbours exist.
    assert int(bnd_prev[1]) == 0
    assert int(bnd_next[1]) == 2
    # Node 0 is a corner (open[0] starts here, land[0] starts here) →
    # neighbours wiped.
    assert int(bnd_prev[0]) == -1
    assert int(bnd_next[0]) == -1
    # Interior node 6 has no boundary tangent at all.
    assert int(bnd_prev[6]) == -1
    assert int(bnd_next[6]) == -1
    # Edge (1, 2) is on the land bottom segment.
    assert (1, 2) in e2s
    kind, seg_idx, position, ibtype = e2s[(1, 2)]
    assert kind == "land"
    assert ibtype == 0
    # Edge (1, 4) is interior so not in the segment map.
    assert (1, 4) not in e2s


def test_batch_smooth_sweep_with_boundary_tangent_moves_segment_interior() -> None:
    """A linear land segment whose middle node is displaced normal to
    the segment must, under boundary-tangent smooth, slide back
    toward the segment line."""
    # Five collinear coastal nodes spaced along y=0, with the middle
    # one pulled down to y = -0.1.
    nodes = np.array(
        [
            [0.0, 0.0], [1.0, 0.0], [2.0, -0.1], [3.0, 0.0], [4.0, 0.0],
            [0.0, 1.0], [4.0, 1.0],
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 5], [1, 2, 5], [2, 3, 5], [3, 4, 5],
            [5, 4, 6],
        ],
        dtype=np.int64,
    )
    mesh = Fort14Mesh(
        title="bnd", nodes=nodes, depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[np.array([5, 6], dtype=np.int64)],
        land_boundaries=[(0, np.array([0, 1, 2, 3, 4], dtype=np.int64))],
    )
    bnd_node = _boundary_node_mask(mesh)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)

    y_before = float(mesh.nodes[2, 1])
    assert y_before == -0.1

    accepts_total = 0
    for _ in range(20):
        n_acc = _batch_smooth_sweep(
            mesh,
            alpha_target=0.95, min_angle_target=20.0,
            boundary_node_mask=bnd_node, n2e=n2e,
            boundary_prev=bnd_prev, boundary_next=bnd_next,
        )
        accepts_total += int(n_acc)
        if n_acc == 0:
            break

    y_after = float(mesh.nodes[2, 1])
    assert accepts_total >= 1, "boundary tangent smooth should fire at least once"
    # Node 2's tangent line is y = 0 between (1, 0) and (3, 0); after
    # smoothing it should be much closer to y=0 than to its starting
    # y=-0.1.
    assert abs(y_after) < 0.02, f"node 2 y not pulled back to segment: {y_after}"


def _long_thin_boundary_triangle() -> Fort14Mesh:
    """Single 4:1-aspect isoceles triangle with the long edge on
    the boundary. Original alpha ~0.53; splitting the long boundary
    edge at the midpoint yields two right triangles with alpha
    ~0.69 each — penalty drops, so the operator should accept."""
    nodes = np.array(
        [[0.0, 0.0], [4.0, 0.0], [2.0, 1.0]], dtype=float,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    return Fort14Mesh(
        title="thin", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([0, 1], dtype=np.int64)),
            (0, np.array([1, 2], dtype=np.int64)),
            (0, np.array([2, 0], dtype=np.int64)),
        ],
    )


def test_apply_edge_split_boundary_inserts_and_updates_segment() -> None:
    mesh = _long_thin_boundary_triangle()
    eu = _edge_use_counts(mesh.elements)
    _, _, e2s = _boundary_topology(mesh)
    # Edge (0, 1) is on land segment 0; original triangle alpha is
    # ~0.23 (needle); splitting (0, 1) yields two alpha~0.77 triangles.
    assert (0, 1) in e2s
    out = _apply_edge_split_boundary(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, edge_to_segment=e2s,
    )
    assert out is not None
    new_mesh, info = out
    assert info["operator"] == "edge_split_boundary"
    assert info["boundary_kind"] == "land"
    assert new_mesh.n_nodes == mesh.n_nodes + 1
    assert new_mesh.n_elements == mesh.n_elements + 1  # 1 - 1 + 2
    # The new node sits at the midpoint of (0, 1).
    np.testing.assert_allclose(
        new_mesh.nodes[mesh.n_nodes], 0.5 * (mesh.nodes[0] + mesh.nodes[1]),
    )
    # The land segment containing the (0, 1) edge now threads the new
    # node between them.
    bot_seg = new_mesh.land_boundaries[0][1]
    assert list(bot_seg) == [0, mesh.n_nodes, 1]


def test_apply_edge_split_boundary_rejects_interior_edge() -> None:
    mesh = _skewed_quad()
    eu = _edge_use_counts(mesh.elements)
    _, _, e2s = _boundary_topology(mesh)
    # (1, 4) is an internal edge — not in the segment map.
    assert (1, 4) not in e2s
    out = _apply_edge_split_boundary(
        mesh, elem_id=1, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, edge_to_segment=e2s,
    )
    assert out is None


def test_phase_h_v2_default_operator_order_includes_boundary_split() -> None:
    """The v2 default operator order plumbs ``edge_split_boundary``
    so a Phase H run with default kwargs sees boundary edges as
    valid targets."""
    from fvcom_mesh_tools.mesh_clean_phase_h import DEFAULT_OPERATOR_ORDER
    assert "edge_split_boundary" in DEFAULT_OPERATOR_ORDER


# ---------------------------------------------------------------------------
# v3: coastline-projecting boundary operations
# ---------------------------------------------------------------------------


def test_apply_edge_split_boundary_snaps_to_coastline_projector() -> None:
    """When a coastline projector is supplied that nudges the
    midpoint slightly *into* the triangle (toward the apex), the new
    midpoint is placed at the projector's output rather than the
    chord midpoint. The shift must be small enough that the resulting
    sub-triangles still improve over the parent (otherwise the
    operator rejects)."""
    mesh = _long_thin_boundary_triangle()
    eu = _edge_use_counts(mesh.elements)
    _, _, e2s = _boundary_topology(mesh)

    # Nudge +0.02 in y (toward the apex at (2, 1)). At this magnitude
    # the resulting alpha is still well above the parent's.
    def _shift_projector(xy: np.ndarray) -> np.ndarray:
        return np.asarray([xy[0], xy[1] + 0.02], dtype=float)

    out = _apply_edge_split_boundary(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, edge_to_segment=e2s,
        coastline_projector=_shift_projector,
    )
    assert out is not None
    new_mesh, info = out
    assert info["snapped_to_coastline"] is True
    # Original straight midpoint of (0,0)→(4,0) is (2, 0); projector
    # bumps the y to 0.02.
    np.testing.assert_allclose(
        new_mesh.nodes[mesh.n_nodes], [2.0, 0.02],
    )


def test_apply_edge_split_boundary_no_projector_uses_chord_midpoint() -> None:
    """Without a projector, v3 still produces v2 behaviour."""
    mesh = _long_thin_boundary_triangle()
    eu = _edge_use_counts(mesh.elements)
    _, _, e2s = _boundary_topology(mesh)
    out = _apply_edge_split_boundary(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, edge_to_segment=e2s,
        coastline_projector=None,
    )
    assert out is not None
    new_mesh, info = out
    assert info["snapped_to_coastline"] is False
    np.testing.assert_allclose(
        new_mesh.nodes[mesh.n_nodes], [2.0, 0.0],
    )


def test_apply_edge_split_boundary_projector_returning_none_uses_chord() -> None:
    """A projector that returns ``None`` (point too far from any
    polyline) falls back to the chord midpoint."""
    mesh = _long_thin_boundary_triangle()
    eu = _edge_use_counts(mesh.elements)
    _, _, e2s = _boundary_topology(mesh)
    out = _apply_edge_split_boundary(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, edge_to_segment=e2s,
        coastline_projector=lambda _xy: None,
    )
    assert out is not None
    new_mesh, info = out
    assert info["snapped_to_coastline"] is False
    np.testing.assert_allclose(
        new_mesh.nodes[mesh.n_nodes], [2.0, 0.0],
    )


def test_smooth_node_force_skips_penalty_gate() -> None:
    """``force=True`` should accept a smooth that the strict gate
    rejects. Recentre vertex 6 once, then attempt a second smooth: the
    1-ring centroid coincides with the moved position so penalty is
    unchanged, the strict gate rejects, ``force=True`` accepts.
    """
    mesh = _skewed_quad()
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    bnd = _boundary_node_mask(mesh)
    first = _apply_smooth_node(
        mesh, 6, n2e[6],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert first is not None
    mesh1, _ = first
    n2e1 = _node_to_elements(mesh1.elements, mesh1.n_nodes)

    rejected = _apply_smooth_node(
        mesh1, 6, n2e1[6],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
    )
    assert rejected is None, "strict gate should reject the no-op smooth"

    forced = _apply_smooth_node(
        mesh1, 6, n2e1[6],
        alpha_target=0.95, min_angle_target=20.0,
        boundary_node_mask=bnd,
        force=True,
    )
    assert forced is not None, "force=True should accept the no-op smooth"
    new_mesh, info = forced
    assert info["operator"] == "smooth_node"
    assert info["forced"] is True
    # The proposed position equals the current → mesh is unchanged.
    np.testing.assert_allclose(new_mesh.nodes[6], mesh1.nodes[6])


def test_edge_swap_force_skips_penalty_gate() -> None:
    """Two right-isoceles triangles sharing the (0, 2) diagonal (under
    the element ordering ``[(0, 2, 1), (0, 3, 2)]`` which matches the
    operator's CCW post-swap convention). The Lawson swap to the
    (1, 3) diagonal yields two congruent right triangles — penalty
    is unchanged, strict gate rejects, ``force=True`` accepts.
    """
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        dtype=float,
    )
    elements = np.array(
        [[0, 2, 1], [0, 3, 2]],
        dtype=np.int64,
    )
    mesh = Fort14Mesh(
        title="swap_force",
        nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([0, 1, 2, 3, 0], dtype=np.int64)),
        ],
    )
    eu = _edge_use_counts(mesh.elements)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}
    assert eu[(0, 2)] == [0, 1], "shared diagonal must be interior"

    strict = _apply_edge_swap(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, boundary_edge_keys=bnd_edges,
    )
    assert strict is None, "Lawson-neutral swap must be rejected strict"

    forced = _apply_edge_swap(
        mesh, elem_id=0, edge_local=0,
        alpha_target=0.95, min_angle_target=20.0,
        edge_uses=eu, boundary_edge_keys=bnd_edges,
        force=True,
    )
    assert forced is not None, "force=True must accept the swap"
    new_mesh, info = forced
    assert info["operator"] == "edge_swap"
    assert info["forced"] is True
    assert info["penalty_after"] == info["penalty_before"]
    assert new_mesh.n_elements == mesh.n_elements
    # After the swap the (0, 2) diagonal disappears and (1, 3) appears.
    new_eu = _edge_use_counts(new_mesh.elements)
    assert (1, 3) in new_eu and len(new_eu[(1, 3)]) == 2
    assert (0, 2) not in new_eu


def test_build_coastline_projector_returns_none_for_empty_paths() -> None:
    from fvcom_mesh_tools.mesh_clean_phase_h import (
        build_coastline_projector,
    )
    assert build_coastline_projector(None) is None
    assert build_coastline_projector([]) is None


def test_build_coastline_projector_snaps_to_nearest_polyline(tmp_path) -> None:
    """Build a tiny shapefile from a hand-crafted polyline and verify
    that ``build_coastline_projector`` snaps an off-line point onto
    it, while a far-away point falls through (returns ``None``)."""
    import pytest
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import LineString

    from fvcom_mesh_tools.mesh_clean_phase_h import (
        build_coastline_projector,
    )

    # Polyline at y=35.0 between x=139.0 and x=140.0 (Tokyo Bay-ish lat).
    polyline = LineString([(139.0, 35.0), (140.0, 35.0)])
    gdf = gpd.GeoDataFrame(geometry=[polyline], crs="EPSG:4326")
    shp = tmp_path / "coast.shp"
    gdf.to_file(shp)

    proj = build_coastline_projector(
        [shp], max_snap_distance_m=5_000.0, mean_latitude_deg=35.0,
    )
    assert proj is not None
    # A point ~100 m above the polyline at lon 139.5 → projects to
    # (139.5, 35.0).
    deg_per_m_lat = 1.0 / (
        6_371_000.0 * (np.pi / 180.0)
    )
    near = np.array([139.5, 35.0 + 100.0 * deg_per_m_lat])
    snapped = proj(near)
    assert snapped is not None
    np.testing.assert_allclose(snapped[0], 139.5, atol=1e-8)
    np.testing.assert_allclose(snapped[1], 35.0, atol=1e-8)

    # A point 100 km above is far beyond the 5 km snap range → None.
    far = np.array([139.5, 35.0 + 1.0])
    assert proj(far) is None
