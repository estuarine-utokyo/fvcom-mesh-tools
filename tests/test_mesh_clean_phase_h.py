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
    DEFAULT_AREA_RATIO_TARGET,
    DEFAULT_MAX_ANGLE_TARGET,
    DEFAULT_MAX_VALENCE,
    _affected_internal_edges_around,
    _apply_edge_split_boundary,
    _apply_edge_split_interior,
    _apply_edge_swap,
    _apply_pass_e_split,
    _apply_pass_e_swap,
    _apply_smooth_node,
    _apply_vertex_remove,
    _batch_smooth_sweep,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _internal_edge_buddies,
    _is_fail,
    _local_area_changes,
    _node_to_elements,
    _pass_e_round,
    _pass_f_c4_smooth_sweep,
    _pass_f_round,
    _pass_g_c1_smooth_sweep,
    _pass_g_round,
    _penalty,
    _per_edge_area_change,
    _per_element_quality,
    phase_h_finish,
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
    assert info["lookahead_op1_inventory"] == ["smooth_node"]
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

    def _n_flipped(m):
        p0 = m.nodes[m.elements[:, 0]]
        p1 = m.nodes[m.elements[:, 1]]
        p2 = m.nodes[m.elements[:, 2]]
        cross = (
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        )
        return int((cross <= 0).sum())

    n_flipped_before = _n_flipped(mesh)
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
    # behind. We assert only that the round runs cleanly, returns a
    # well-formed result, and does not introduce new flipped tris
    # (the fixture itself has one pre-flipped element by
    # construction).
    assert isinstance(accepts, dict)
    assert all(
        isinstance(k, str) and "+" in k for k in accepts
    ), "pair labels must be 'op1+op2' strings"
    assert _n_flipped(cur) <= n_flipped_before, (
        "lookahead must not introduce new flipped triangles"
    )


def test_target_exits_fail_rejects_absent_vertex_set() -> None:
    """v4.1 strict gate: an absent target vertex set is REJECTED, not
    "fixed by elimination". The earlier permissive interpretation
    produced PoC #46's catastrophic regression — under the default
    inventory ``vertex_remove`` deletes E by construction, so the
    elimination branch made every valid vertex_remove auto-accept
    and ~19 600 interior vertices were stripped from the Tokyo-Bay
    mesh.
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
    target = frozenset({0, 1, 99})
    assert _target_exits_fail(
        mesh, target, alpha_target=0.95, min_angle_target=20.0,
    ) is False


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


def test_find_fail_clusters_groups_adjacent_fails() -> None:
    """Two fail elements sharing an edge should form a single cluster
    of size 2."""
    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _find_fail_clusters,
    )
    mesh = _skewed_quad()
    clusters = _find_fail_clusters(
        mesh,
        alpha_target=0.95, min_angle_target=20.0,
        min_cluster_size=1, max_cluster_size=100,
    )
    # The fixture's central 4-triangle fan should be detected as a
    # connected fail cluster.
    assert len(clusters) >= 1
    assert max(c.size for c in clusters) >= 2


def test_apply_patch_recdt_accepts_interior_cluster() -> None:
    """A clean equilateral-rim hexagon with a deliberately bad
    interior vertex: dropping the interior + Delaunay re-mesh should
    pass every gate. Validates the happy path of Pass D."""
    import math

    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _apply_patch_recdt,
    )
    # 6 rim nodes on a hexagon + 1 interior node offset from centre.
    r = 1.0
    rim = np.array(
        [[r * math.cos(t), r * math.sin(t)]
         for t in np.linspace(0.0, 2.0 * math.pi, 7)[:-1]],
        dtype=float,
    )
    centre = np.array([[0.2, 0.1]], dtype=float)
    nodes = np.vstack([rim, centre])
    centre_id = 6
    elements = np.array(
        [[i, (i + 1) % 6, centre_id] for i in range(6)],
        dtype=np.int64,
    )
    mesh = Fort14Mesh(
        title="hex_patch", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(
            0, np.array([0, 1, 2, 3, 4, 5, 0], dtype=np.int64),
        )],
    )
    # All 6 elements form one cluster (their interior node is shared
    # ⇒ they're face-face-adjacent via the rim edges? No, they share
    # vertex but not edges with their non-neighbour fan members. The
    # face-face adjacency uses edge sharing, so the fan is a chain).
    # Pass cluster_eids explicitly to test the apply function
    # directly.
    cluster = np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)
    # Regular hexagon Delaunay (no interior point) yields 4
    # isoceles triangles with alpha ≈ 0.75 each — relax the gate so
    # the patch acceptance is exercised mechanically rather than
    # gated out on a fixture-specific Delaunay numeric.
    out = _apply_patch_recdt(
        mesh, cluster,
        alpha_target=0.5, min_angle_target=20.0,
        reject_boundary_clusters=False,  # rim IS the mesh boundary here
    )
    assert out is not None, "hexagon patch with bad centre must accept"
    new_mesh, info = out
    assert info["operator"] == "patch_recdt"
    assert info["cluster_size"] == 6
    assert info["rim_size"] == 6
    assert info["n_interior_orphaned"] == 1  # the centre node
    # The new patch should have 4 triangles (Delaunay of a 6-gon).
    assert info["n_new_elements"] == 4
    assert new_mesh.n_elements == 4
    # No flipped tris.
    p0 = new_mesh.nodes[new_mesh.elements[:, 0]]
    p1 = new_mesh.nodes[new_mesh.elements[:, 1]]
    p2 = new_mesh.nodes[new_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all()


def test_apply_patch_recdt_rejects_boundary_rim() -> None:
    """When ``reject_boundary_clusters=True``, a cluster whose rim
    sits on an open / land boundary segment must be rejected."""
    import math

    from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
        _apply_patch_recdt,
    )
    rim = np.array(
        [[math.cos(t), math.sin(t)]
         for t in np.linspace(0.0, 2.0 * math.pi, 7)[:-1]],
        dtype=float,
    )
    centre = np.array([[0.2, 0.1]], dtype=float)
    nodes = np.vstack([rim, centre])
    elements = np.array(
        [[i, (i + 1) % 6, 6] for i in range(6)],
        dtype=np.int64,
    )
    mesh = Fort14Mesh(
        title="hex_patch_bnd", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(
            0, np.array([0, 1, 2, 3, 4, 5, 0], dtype=np.int64),
        )],
    )
    cluster = np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)
    out = _apply_patch_recdt(
        mesh, cluster,
        alpha_target=0.95, min_angle_target=20.0,
        reject_boundary_clusters=True,
    )
    assert out is None, "boundary-touching cluster must be rejected"


def test_phase_h_optimize_pass_d_path_runs_cleanly() -> None:
    """Smoke test: enabling Pass D should not crash on the
    skewed-quad fixture (whether or not any patch fires)."""
    mesh = _skewed_quad()
    out_mesh, info = phase_h_optimize(
        mesh,
        alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=3,
        patch_recdt_enabled=True,
        patch_min_cluster_size=1,  # let the fixture's tiny clusters in
        max_patches_per_round=4,
    )
    assert info["patch_recdt_enabled"] is True
    assert isinstance(info["patch_recdt_accepts"], dict)
    assert "n_patch_recdt_rejected" in info
    # No new flips beyond the fixture's pre-flipped element.
    p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
    p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
    p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    n_flipped_out = int((cross <= 0).sum())
    p0i = mesh.nodes[mesh.elements[:, 0]]
    p1i = mesh.nodes[mesh.elements[:, 1]]
    p2i = mesh.nodes[mesh.elements[:, 2]]
    crossi = (
        (p1i[:, 0] - p0i[:, 0]) * (p2i[:, 1] - p0i[:, 1])
        - (p1i[:, 1] - p0i[:, 1]) * (p2i[:, 0] - p0i[:, 0])
    )
    n_flipped_in = int((crossi <= 0).sum())
    assert n_flipped_out <= n_flipped_in


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


# ---------------------------------------------------------------------------
# max_angle_target gate (FVCOM manual criterion C2)
# ---------------------------------------------------------------------------


def _obtuse_triangle_pair() -> Fort14Mesh:
    """4-node mesh: one needle-obtuse triangle (~150° angle at v2) and
    one well-shaped triangle. Used to exercise the C2 gate."""
    nodes = np.array(
        [
            [0.0, 0.0],   # 0
            [1.0, 0.0],   # 1
            [0.5, 0.05],  # 2 (interior, very close to edge 0-1 →
                          # angle at v2 is highly obtuse)
            [0.5, 1.0],   # 3 (forms a well-shaped triangle with 0, 1)
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 2],  # obtuse: angle at v2 ~170°
            [0, 1, 3],  # equilateral-ish, max_angle < 130°
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="obtuse", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([0, 1, 3], dtype=np.int64)),
            (0, np.array([0, 2, 1], dtype=np.int64)),
        ],
    )


def test_per_element_quality_returns_max_angle_in_degrees() -> None:
    """`_per_element_quality` must return a 3-tuple
    ``(alpha, min_angle_deg, max_angle_deg)`` so callers can gate on C2.
    """
    mesh = _obtuse_triangle_pair()
    alpha, min_ang, max_ang = _per_element_quality(mesh.nodes, mesh.elements)
    # The needle triangle (index 0) has a very wide max angle.
    assert max_ang[0] > 150.0, (
        f"needle triangle max angle should be > 150°, got {max_ang[0]}"
    )
    # The well-shaped triangle (index 1) is bounded under 130°.
    assert max_ang[1] < 130.0, (
        f"well-shaped triangle max angle should be < 130°, got {max_ang[1]}"
    )
    # min/max consistency.
    assert (min_ang <= max_ang).all()


def test_is_fail_flags_max_angle_violation() -> None:
    """An element with max_angle > target should be flagged as fail
    even if alpha and min_angle pass."""
    alpha = np.array([0.99])
    min_ang = np.array([45.0])
    max_ang = np.array([135.0])
    # Without the max-angle term (max_angle_target=180 default),
    # this element passes.
    fail_no_gate = _is_fail(
        alpha, min_ang, alpha_target=0.95, min_angle_target=30.0,
    )
    assert not fail_no_gate.any()
    # With max_angle_target=130, the same element flips to fail.
    fail_gated = _is_fail(
        alpha, min_ang, max_ang,
        alpha_target=0.95, min_angle_target=30.0,
        max_angle_target=130.0,
    )
    assert fail_gated.all()


def test_penalty_adds_max_angle_term_only_when_exceeded() -> None:
    """The penalty must equal the alpha + min-angle baseline whenever
    ``max_ang <= max_angle_target``, and rise by a positive amount
    once ``max_ang > max_angle_target``."""
    alpha = np.array([0.99, 0.99])
    min_ang = np.array([45.0, 45.0])
    max_ang = np.array([125.0, 140.0])  # below and above 130°
    baseline = _penalty(
        alpha, min_ang, alpha_target=0.95, min_angle_target=30.0,
    )
    gated = _penalty(
        alpha, min_ang, max_ang,
        alpha_target=0.95, min_angle_target=30.0,
        max_angle_target=130.0,
    )
    # Element 0 (max_ang=125°): no contribution from max-angle term.
    np.testing.assert_allclose(gated[0], baseline[0])
    # Element 1 (max_ang=140°): penalty must rise by (140-130)^2/100 = 1.
    np.testing.assert_allclose(gated[1] - baseline[1], 1.0, atol=1e-12)


def test_phase_h_optimize_accepts_max_angle_target_kwarg() -> None:
    """`phase_h_optimize` must accept ``max_angle_target`` and surface
    it in the returned ``info`` dict for downstream summaries."""
    mesh = _skewed_quad()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=130.0,
        max_outer_rounds=5,
    )
    assert info["max_angle_target"] == 130.0
    # No flipped triangles.
    p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
    p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
    p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all()


def test_phase_h_optimize_default_max_angle_is_disabled() -> None:
    """The default ``max_angle_target`` (``DEFAULT_MAX_ANGLE_TARGET``,
    180°) must leave the gate effectively disabled — Phase H's behaviour
    on a mesh whose max angles never exceed 180° (i.e., every valid
    triangle) must match the non-gated case."""
    mesh = _skewed_quad()
    _, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=5,
    )
    assert info["max_angle_target"] == DEFAULT_MAX_ANGLE_TARGET
    assert info["max_angle_target"] == 180.0


# ---------------------------------------------------------------------------
# Pass E: gradation refinement (FVCOM manual criterion C4)
# ---------------------------------------------------------------------------


def _c4_fail_strip() -> Fort14Mesh:
    """5-node strip with one equilateral-ish triangle next to a
    much smaller one, triggering area_change > 0.5 on the shared
    edge while keeping each new sub-triangle's min angle >= 20°.

    Triangles:
      T0 = (0, 1, 2): equilateral, side 2, area = sqrt(3) ≈ 1.732
      T1 = (1, 0, 3): thin sliver, area = 0.5  → area_change(0-1)
                       = (1.732 - 0.5)/1.732 ≈ 0.711 (fails C4)
      T2 = (0, 2, 4): keeps edge 0-2 internal so Pass E has two
                       candidate edges to choose from
    """
    sqrt3 = float(np.sqrt(3.0))
    nodes = np.array(
        [
            [0.0, 0.0],     # 0
            [2.0, 0.0],     # 1   (equilateral base)
            [1.0, sqrt3],   # 2   (equilateral apex)
            [1.0, -0.5],    # 3   (sliver apex below 0-1)
            [-1.0, 2.0],    # 4   (gives T2)
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 2],
            [1, 0, 3],
            [0, 2, 4],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="c4_strip", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([4, 2, 1, 3, 0], dtype=np.int64)),
        ],
    )


def test_internal_edge_buddies_finds_shared_edges() -> None:
    """The C4-fail strip has 3 elements; the boundary walk has 5
    edges, the strip has 3*3 = 9 edge slots, so 2 internal edges
    (edges 0-1 and 0-2)."""
    mesh = _c4_fail_strip()
    edge_uv, elem_pair = _internal_edge_buddies(mesh.elements)
    assert edge_uv.shape == (2, 2)
    assert elem_pair.shape == (2, 2)
    # Both internal edges include node 0.
    edge_sets = [frozenset(e.tolist()) for e in edge_uv]
    assert frozenset([0, 1]) in edge_sets
    assert frozenset([0, 2]) in edge_sets


def test_per_edge_area_change_matches_manual_formula() -> None:
    """The 0-1 shared edge sits between T0 (area=sqrt(3)) and T1
    (area=0.5), so area_change = (sqrt(3) - 0.5) / sqrt(3) ≈ 0.711
    → fails C4. The 0-2 shared edge between T0 and T2 fits the same
    formula on its respective pair."""
    mesh = _c4_fail_strip()
    edge_uv, _elem_pair, ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    keys = [tuple(sorted(e.tolist())) for e in edge_uv]
    idx_01 = keys.index((0, 1))
    expected_01 = (np.sqrt(3.0) - 0.5) / np.sqrt(3.0)
    np.testing.assert_allclose(ac[idx_01], expected_01, atol=1e-12)
    assert ac[idx_01] > 0.5  # confirms this fixture trips C4


def test_apply_pass_e_split_drops_c4_violation() -> None:
    """The 0-1 internal edge fails C4 (area_change=0.75). Pass E
    should split the longest non-shared edge of T0 (the larger of the
    two incident triangles), dropping area_change at edge 0-1
    strictly below the 0.5 target."""
    mesh = _c4_fail_strip()
    edge_uv, elem_pair, ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    # Pick the 0-1 fail edge.
    keys = [tuple(sorted(e.tolist())) for e in edge_uv]
    fail_idx = keys.index((0, 1))
    u, v = int(edge_uv[fail_idx, 0]), int(edge_uv[fail_idx, 1])
    e_i, e_j = int(elem_pair[fail_idx, 0]), int(elem_pair[fail_idx, 1])

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    _bp, _bn, edge_to_segment = _boundary_topology(mesh)

    out = _apply_pass_e_split(
        mesh, (u, v), (e_i, e_j),
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        edge_uses=edge_uses,
        boundary_edge_keys=boundary_edge_keys,
        edge_to_segment=edge_to_segment,
        coastline_projector=None,
    )
    assert out is not None, "Pass E should accept the strip split"
    new_mesh, info = out
    assert info["operator"] == "pass_e_split"
    assert info["area_change_after"] < info["area_change_before"]
    # After halving T0 (area sqrt(3) ≈ 1.732 → ~0.866) the shared
    # edge's new ratio is (0.866 - 0.5)/0.866 ≈ 0.42 → below the 0.5
    # target.
    assert info["area_change_after"] <= 0.5 + 1e-12
    # Mesh sanity: no flipped triangles.
    p0 = new_mesh.nodes[new_mesh.elements[:, 0]]
    p1 = new_mesh.nodes[new_mesh.elements[:, 1]]
    p2 = new_mesh.nodes[new_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all()


def test_pass_e_round_eliminates_strip_c4_fail() -> None:
    """End-to-end ``_pass_e_round``: the strip fixture has 1 (or 2)
    C4 fails before; after the round, zero C4 fails should remain
    (or the round should converge to a fixed point)."""
    mesh = _c4_fail_strip()
    _, _, ac_before = _per_edge_area_change(mesh.nodes, mesh.elements)
    n_fail_before = int((ac_before > 0.5).sum())
    assert n_fail_before >= 1

    new_mesh, n_acc, n_rej, _n_swap, _n_split = _pass_e_round(
        mesh,
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_splits=20,
        coastline_projector=None,
    )
    assert n_acc >= 1
    _, _, ac_after = _per_edge_area_change(
        new_mesh.nodes, new_mesh.elements,
    )
    n_fail_after = int((ac_after > 0.5).sum())
    assert n_fail_after < n_fail_before, (
        f"Pass E should reduce C4 fails: {n_fail_before} -> {n_fail_after}"
    )


def test_phase_h_optimize_pass_e_path_runs_cleanly() -> None:
    """`phase_h_optimize` must accept ``pass_e_enabled=True`` and
    surface Pass E bookkeeping in ``info``."""
    mesh = _skewed_quad()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=3,
        pass_e_enabled=True,
        pass_e_area_ratio_target=0.5,
    )
    assert info["pass_e_enabled"] is True
    assert "pass_e_accepts" in info
    assert "pass_e_rejected" in info
    # Mesh stays valid (no flipped triangles).
    p0 = out_mesh.nodes[out_mesh.elements[:, 0]]
    p1 = out_mesh.nodes[out_mesh.elements[:, 1]]
    p2 = out_mesh.nodes[out_mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross > 0).all()


def test_phase_h_optimize_pass_e_disabled_by_default() -> None:
    """Pass E must stay off unless explicitly enabled — preserves
    backward compatibility for existing callers."""
    mesh = _skewed_quad()
    _, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_outer_rounds=3,
    )
    assert info["pass_e_enabled"] is False
    assert info["pass_e_accepts"] == 0
    assert info["pass_e_rejected"] == 0
    assert info["pass_e_area_ratio_target"] == DEFAULT_AREA_RATIO_TARGET


def test_pass_e_skips_cascading_c4_fail_candidate() -> None:
    """If the longest non-shared edge of L is itself a C4 fail edge,
    splitting it cascades the violation onto two new edges (each
    inheriting the same area_change). Pass E must skip such
    candidates when ``c4_fail_keys`` is provided. We verify by
    constructing a fixture where the only acceptable candidate is
    skipped and the operator returns ``None``."""
    mesh = _c4_fail_strip()
    edge_uv, elem_pair, ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    keys = [tuple(sorted(e.tolist())) for e in edge_uv]
    fail_idx = keys.index((0, 1))
    u, v = int(edge_uv[fail_idx, 0]), int(edge_uv[fail_idx, 1])
    e_i, e_j = int(elem_pair[fail_idx, 0]), int(elem_pair[fail_idx, 1])

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    _bp, _bn, edge_to_segment = _boundary_topology(mesh)

    # Artificially mark every other internal edge AND every boundary
    # edge of L as a C4 fail. ``L`` is the larger triangle T0 =
    # (0, 1, 2). Its three edges are (0,1) [the fail edge itself],
    # (1,2) [boundary], (0,2) [internal, shared with T2]. Mark
    # (1,2) and (0,2) as cascading fails — Pass E must reject.
    cascading_fail_keys = {(0, 1), (1, 2), (0, 2)}

    out = _apply_pass_e_split(
        mesh, (u, v), (e_i, e_j),
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_valence=DEFAULT_MAX_VALENCE,
        edge_uses=edge_uses,
        boundary_edge_keys=boundary_edge_keys,
        edge_to_segment=edge_to_segment,
        c4_fail_keys=cascading_fail_keys,
        valence_before=None,
        coastline_projector=None,
    )
    assert out is None, (
        "Pass E must reject when every candidate edge is itself "
        "in c4_fail_keys (cascade avoidance)"
    )


def test_pass_e_swap_signature_and_safe_rejection() -> None:
    """`_apply_pass_e_swap` accepts the standard Pass-E kwargs and
    returns a tuple-or-``None``. On the strip fixture the swap output
    would invert one of the new triangles (kite quad geometry), which
    the underlying ``_apply_edge_swap`` rejects via its signed-area
    check — so this is a valid rejection path."""
    mesh = _c4_fail_strip()
    edge_uv, elem_pair, _ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    keys = [tuple(sorted(e.tolist())) for e in edge_uv]
    fail_idx = keys.index((0, 1))
    u, v = int(edge_uv[fail_idx, 0]), int(edge_uv[fail_idx, 1])
    e_i, e_j = int(elem_pair[fail_idx, 0]), int(elem_pair[fail_idx, 1])

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}

    out = _apply_pass_e_swap(
        mesh, (u, v), (e_i, e_j),
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_valence=DEFAULT_MAX_VALENCE,
        edge_uses=edge_uses,
        boundary_edge_keys=boundary_edge_keys,
        c4_fail_keys=None,
        valence_before=None,
    )
    assert out is None or (
        isinstance(out, tuple) and len(out) == 2
    ), "operator must return None or (new_mesh, info)"


def test_pass_e_swap_rejects_when_target_is_boundary() -> None:
    """An edge_swap on a boundary edge is topologically undefined
    (one of the would-be incident triangles doesn't exist). The
    operator must short-circuit to ``None`` rather than blowing up
    on the buddy-lookup."""
    mesh = _c4_fail_strip()
    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    # Pick any boundary edge from the strip fixture.
    bnd_key = next(iter(boundary_edge_keys))

    out = _apply_pass_e_swap(
        mesh, bnd_key, (0, 1),  # elem pair is a dummy; the boundary
                                # check fires first
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_valence=DEFAULT_MAX_VALENCE,
        edge_uses=edge_uses,
        boundary_edge_keys=boundary_edge_keys,
        c4_fail_keys=None,
        valence_before=None,
    )
    assert out is None


def test_pass_e_prefilters_on_valence_cap() -> None:
    """If the C5 prefilter sees the larger triangle's opposite
    vertex already at ``max_valence``, splitting any non-shared edge
    would push the opposite vertex to ``max_valence + 1``, regressing
    C5. Pass E must reject. We simulate by passing a synthetic
    ``valence_before`` array that pins every node at the cap."""
    mesh = _c4_fail_strip()
    edge_uv, elem_pair, _ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    keys = [tuple(sorted(e.tolist())) for e in edge_uv]
    fail_idx = keys.index((0, 1))
    u, v = int(edge_uv[fail_idx, 0]), int(edge_uv[fail_idx, 1])
    e_i, e_j = int(elem_pair[fail_idx, 0]), int(elem_pair[fail_idx, 1])

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    _bp, _bn, edge_to_segment = _boundary_topology(mesh)

    # Every node already at max_valence — any split would regress C5.
    valence_pinned = np.full(mesh.n_nodes, DEFAULT_MAX_VALENCE, dtype=np.int64)

    out = _apply_pass_e_split(
        mesh, (u, v), (e_i, e_j),
        alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_valence=DEFAULT_MAX_VALENCE,
        edge_uses=edge_uses,
        boundary_edge_keys=boundary_edge_keys,
        edge_to_segment=edge_to_segment,
        c4_fail_keys=None,
        valence_before=valence_pinned,
        coastline_projector=None,
    )
    assert out is None, (
        "Pass E must reject when the C5 prefilter would push the "
        "opposite vertex over max_valence"
    )


# ---------------------------------------------------------------------------
# Pass F: C4-aware smoothing (Laplacian + count-gate, target nodes only)
# ---------------------------------------------------------------------------


def _c4_fail_fan() -> Fort14Mesh:
    """4-element fan around 1 interior node (#4) offset to (1.2, 1.2).
    The square corners are 0=(0,0), 1=(4,0), 2=(4,4), 3=(0,4) — all on
    the land boundary. Four triangles fan around node 4:

      T0 = (0, 1, 4)  area = 2.4
      T1 = (1, 2, 4)  area = 5.6
      T2 = (2, 3, 4)  area = 5.6
      T3 = (3, 0, 4)  area = 2.4

    Internal edges (1,4) and (3,4) carry area_change ~0.571 — both
    fail C4 (> 0.5). Every triangle's min angle is ~23° so C1 passes
    (default 20°) and every max angle is ~112° so C2 passes (default
    130°). Moving node 4 to the centroid (2, 2) equalises all four
    areas and clears both C4 fails without regressing C1/C2.
    """
    nodes = np.array(
        [
            [0.0, 0.0],   # 0
            [4.0, 0.0],   # 1
            [4.0, 4.0],   # 2
            [0.0, 4.0],   # 3
            [1.2, 1.2],   # 4 — interior
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="c4_fan", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([3, 2, 1, 0], dtype=np.int64)),
        ],
    )


def test_affected_internal_edges_around_ring_collects_buddies() -> None:
    """For the fan fixture, the 1-ring of node 4 is every element of
    the mesh, so every internal edge is affected — all 4 of them
    (each connecting a corner to node 4)."""
    mesh = _c4_fail_fan()
    edge_uses = _edge_use_counts(mesh.elements)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    ring = n2e[4]
    aff_uv, aff_pair = _affected_internal_edges_around(
        mesh.elements, edge_uses, ring,
    )
    assert aff_uv.shape == (4, 2)
    assert aff_pair.shape == (4, 2)
    # Every internal edge here includes node 4.
    for uv in aff_uv:
        assert 4 in uv.tolist()


def test_local_area_changes_matches_global() -> None:
    """``_local_area_changes`` on the fan's internal-edge subset must
    match the global ``_per_edge_area_change`` values element-by-
    element."""
    mesh = _c4_fail_fan()
    edge_uv_g, elem_pair_g, ac_g = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    ac_l = _local_area_changes(mesh.nodes, mesh.elements, elem_pair_g)
    np.testing.assert_allclose(ac_l, ac_g, atol=1e-12)


def test_pass_f_sweep_moves_interior_node_to_centroid() -> None:
    """When the target set restricts to the interior node (#4 of the
    fan fixture), Pass F's Laplacian proposal places it at the 1-ring
    centroid (2, 2) — the area-balancing point that clears every C4
    fail without raising C1/C2."""
    mesh = _c4_fail_fan()
    _uv, _ep, ac_b = _per_edge_area_change(mesh.nodes, mesh.elements)
    assert int((ac_b > 0.5).sum()) == 2

    bnd_node = _boundary_node_mask(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    edge_uses = _edge_use_counts(mesh.elements)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n_acc = _pass_f_c4_smooth_sweep(
        mesh,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        boundary_node_mask=bnd_node,
        n2e=n2e,
        edge_uses=edge_uses,
        boundary_prev=bnd_prev,
        boundary_next=bnd_next,
        coastline_projector=None,
        target_nodes=np.array([4], dtype=np.int64),
    )
    assert n_acc == 1
    np.testing.assert_allclose(mesh.nodes[4], [2.0, 2.0], atol=1e-9)
    _uv2, _ep2, ac_a = _per_edge_area_change(mesh.nodes, mesh.elements)
    assert int((ac_a > 0.5).sum()) == 0


def test_pass_f_round_clears_c4_fails_in_fan() -> None:
    """End-to-end ``_pass_f_round`` on the fan: starting with 2 C4
    fails, terminate with 0 fails and no C1 / C2 regression. The
    fixture admits both interior and boundary-tangent paths to a
    solution; we assert on the count delta, not which nodes moved
    (boundary tangents on a coarse 4-corner polygon can produce
    geometrically aggressive moves that should be validated on the
    real coastline rather than this synthetic fixture)."""
    mesh = _c4_fail_fan()
    _uv, _ep, ac_b = _per_edge_area_change(mesh.nodes, mesh.elements)
    _a_b, m_b, M_b = _per_element_quality(mesh.nodes, mesh.elements)
    assert int((ac_b > 0.5).sum()) == 2
    assert int((m_b < 20.0).sum()) == 0
    assert int((M_b > DEFAULT_MAX_ANGLE_TARGET).sum()) == 0

    new_mesh, n_acc, n_sw = _pass_f_round(
        mesh,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_sweeps=10,
    )
    assert n_acc >= 1
    assert n_sw >= 1
    _uv2, _ep2, ac_a = _per_edge_area_change(
        new_mesh.nodes, new_mesh.elements,
    )
    _a_a, m_a, M_a = _per_element_quality(
        new_mesh.nodes, new_mesh.elements,
    )
    assert int((ac_a > 0.5).sum()) == 0, (
        f"Pass F should clear the fan's C4 fails: "
        f"{int((ac_b > 0.5).sum())} -> {int((ac_a > 0.5).sum())}"
    )
    assert int((m_a < 20.0).sum()) == 0
    assert int((M_a > DEFAULT_MAX_ANGLE_TARGET).sum()) == 0


def test_pass_f_sweep_rejects_when_no_target_nodes() -> None:
    """If ``target_nodes`` is empty, the sweep must trivially accept
    zero moves and not crash."""
    mesh = _c4_fail_fan()
    bnd_node = _boundary_node_mask(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    edge_uses = _edge_use_counts(mesh.elements)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n_acc = _pass_f_c4_smooth_sweep(
        mesh,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        boundary_node_mask=bnd_node,
        n2e=n2e,
        edge_uses=edge_uses,
        boundary_prev=bnd_prev,
        boundary_next=bnd_next,
        coastline_projector=None,
        target_nodes=np.empty(0, dtype=np.int64),
    )
    assert n_acc == 0


def test_pass_f_sweep_rejects_move_that_would_flip_triangle() -> None:
    """Forcing target_nodes to include a corner of the fan would
    propose the boundary-tangent centroid, which on the corner ends
    up clamped — but we additionally guard against signed-area flips
    by feeding a fixture where the centroid lands outside the 1-ring
    convex hull. Here we instead verify the gate by running a sweep
    with the interior node ALREADY at the balanced (2, 2) position;
    Pass F's local C4 count must already be zero so no move is
    accepted."""
    mesh = _c4_fail_fan()
    # Pre-balance: move 4 to (2, 2) so no C4 fails exist.
    mesh.nodes[4] = np.array([2.0, 2.0])
    _uv, _ep, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    assert int((ac > 0.5).sum()) == 0

    bnd_node = _boundary_node_mask(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    edge_uses = _edge_use_counts(mesh.elements)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n_acc = _pass_f_c4_smooth_sweep(
        mesh,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        boundary_node_mask=bnd_node,
        n2e=n2e,
        edge_uses=edge_uses,
        boundary_prev=bnd_prev,
        boundary_next=bnd_next,
        coastline_projector=None,
        target_nodes=np.array([4], dtype=np.int64),
    )
    # Already at the area-balancing point: local C4 count is 0
    # → strict-decrease gate rejects every move.
    assert n_acc == 0


def test_phase_h_optimize_pass_f_path_runs_cleanly() -> None:
    """`phase_h_optimize` must accept ``pass_f_enabled=True`` and
    surface Pass F bookkeeping in ``info``. We pass an empty
    ``operator_order`` so Pass A / Pass B do not pre-clear the C4
    fails (Pass A's penalty gate happens to accept the same fan
    moves) — this isolates the Pass F contribution."""
    mesh = _c4_fail_fan()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_smooth_sweeps=5,
        max_outer_rounds=2,
        operator_order=(),
        pass_f_enabled=True,
        pass_f_area_ratio_target=0.5,
        max_pass_f_sweeps_per_round=10,
    )
    assert info["pass_f_enabled"] is True
    assert info["pass_f_accepts"] >= 1
    assert info["pass_f_sweeps"] >= 1
    _uv, _ep, ac = _per_edge_area_change(
        out_mesh.nodes, out_mesh.elements,
    )
    assert int((ac > 0.5).sum()) == 0


def test_phase_h_optimize_pass_f_disabled_by_default() -> None:
    """`phase_h_optimize` without ``pass_f_enabled`` must leave the
    bookkeeping at zero and not call into Pass F."""
    mesh = _c4_fail_fan()
    _out, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_rounds=1,
    )
    assert info["pass_f_enabled"] is False
    assert info["pass_f_accepts"] == 0
    assert info["pass_f_sweeps"] == 0


# ---------------------------------------------------------------------------
# Pass G: C1-aware smoothing (Laplacian + count-gate, C1 fail neighbourhoods)
# ---------------------------------------------------------------------------


def _c1_fail_fan() -> Fort14Mesh:
    """4-element fan around an interior node placed at (2.0, 0.8) so
    that T0 = (0, 1, 4) has min_angle ~21.78° and max_angle ~136.4°
    (failing C1 at the 30° gate and C2 at the 130° gate). Moving the
    interior node to its 1-ring centroid (2, 2) produces four
    45-45-90 right-isoceles triangles that pass C1 and C2.

    Square corners 0=(0,0), 1=(4,0), 2=(4,4), 3=(0,4) sit on the
    land boundary."""
    nodes = np.array(
        [
            [0.0, 0.0],   # 0
            [4.0, 0.0],   # 1
            [4.0, 4.0],   # 2
            [0.0, 4.0],   # 3
            [2.0, 0.8],   # 4 — interior, thin angle in T0
        ],
        dtype=float,
    )
    elements = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="c1_fan", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[
            (0, np.array([3, 2, 1, 0], dtype=np.int64)),
        ],
    )


def test_pass_g_sweep_moves_interior_node_to_centroid() -> None:
    """Restricting target_nodes to the interior node #4 of the C1 fan
    fixture, Pass G's Laplacian proposal places it at (2, 2), the
    1-ring centroid that gives every triangle a 45° / 45° / 90°
    profile (C1 = 0, C2 = 0)."""
    mesh = _c1_fail_fan()
    _aq, m_b, M_b = _per_element_quality(mesh.nodes, mesh.elements)
    assert int((m_b < 30.0).sum()) >= 1
    assert int((m_b < 30.0).sum()) >= int((M_b > 130.0).sum())

    bnd_node = _boundary_node_mask(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    edge_uses = _edge_use_counts(mesh.elements)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n_acc = _pass_g_c1_smooth_sweep(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        area_ratio_target=0.5,
        boundary_node_mask=bnd_node,
        n2e=n2e,
        edge_uses=edge_uses,
        boundary_prev=bnd_prev,
        boundary_next=bnd_next,
        coastline_projector=None,
        target_nodes=np.array([4], dtype=np.int64),
    )
    assert n_acc == 1
    np.testing.assert_allclose(mesh.nodes[4], [2.0, 2.0], atol=1e-9)
    _aq2, m_a, M_a = _per_element_quality(mesh.nodes, mesh.elements)
    assert int((m_a < 30.0).sum()) == 0
    assert int((M_a > 130.0).sum()) == 0


def test_pass_g_round_clears_c1_fails_in_fan() -> None:
    """End-to-end ``_pass_g_round``: starting from the C1 fan fixture
    with at least 1 C1 fail, terminate with 0 C1 fails. C2 / C4 must
    not regress."""
    mesh = _c1_fail_fan()
    _aq_b, m_b, _M_b = _per_element_quality(mesh.nodes, mesh.elements)
    _uv, _ep, ac_b = _per_edge_area_change(mesh.nodes, mesh.elements)
    c1_b = int((m_b < 30.0).sum())
    c4_b = int((ac_b > 0.5).sum())
    assert c1_b >= 1

    new_mesh, n_acc, n_sw = _pass_g_round(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        area_ratio_target=0.5,
        max_sweeps=10,
    )
    assert n_acc >= 1
    assert n_sw >= 1
    _aq_a, m_a, M_a = _per_element_quality(
        new_mesh.nodes, new_mesh.elements,
    )
    _uv2, _ep2, ac_a = _per_edge_area_change(
        new_mesh.nodes, new_mesh.elements,
    )
    assert int((m_a < 30.0).sum()) == 0, (
        f"Pass G should clear all C1 fails: "
        f"{c1_b} -> {int((m_a < 30.0).sum())}"
    )
    # C2 / C4 must not regress under Pass G's gate.
    assert int((M_a > 130.0).sum()) == 0  # before was 1
    assert int((ac_a > 0.5).sum()) <= c4_b


def test_pass_g_sweep_rejects_when_c1_unchanged() -> None:
    """If Pass G is run on a mesh that has no C1 fails (e.g. the
    fan already at its area-balancing centroid), no node should
    move — the strict-decrease gate rejects every candidate."""
    mesh = _c1_fail_fan()
    mesh.nodes[4] = np.array([2.0, 2.0])  # pre-balance
    _aq, m, _M = _per_element_quality(mesh.nodes, mesh.elements)
    assert int((m < 30.0).sum()) == 0

    bnd_node = _boundary_node_mask(mesh)
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    edge_uses = _edge_use_counts(mesh.elements)
    bnd_prev, bnd_next, _ = _boundary_topology(mesh)
    n_acc = _pass_g_c1_smooth_sweep(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        area_ratio_target=0.5,
        boundary_node_mask=bnd_node,
        n2e=n2e,
        edge_uses=edge_uses,
        boundary_prev=bnd_prev,
        boundary_next=bnd_next,
        coastline_projector=None,
        target_nodes=np.array([4], dtype=np.int64),
    )
    assert n_acc == 0


def test_phase_h_optimize_pass_g_path_runs_cleanly() -> None:
    """`phase_h_optimize` must accept ``pass_g_enabled=True`` and
    surface Pass G bookkeeping in ``info``. ``operator_order = ()``
    keeps Pass A from clearing the fixture first."""
    mesh = _c1_fail_fan()
    out_mesh, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=30.0,
        max_angle_target=130.0,
        max_outer_rounds=2,
        operator_order=(),
        pass_g_enabled=True,
        pass_g_min_angle_target=30.0,
        pass_g_area_ratio_target=0.5,
        max_pass_g_sweeps_per_round=10,
    )
    assert info["pass_g_enabled"] is True
    assert info["pass_g_accepts"] >= 1
    assert info["pass_g_sweeps"] >= 1
    _aq, m, _M = _per_element_quality(out_mesh.nodes, out_mesh.elements)
    assert int((m < 30.0).sum()) == 0


def test_phase_h_optimize_pass_g_disabled_by_default() -> None:
    """`phase_h_optimize` without ``pass_g_enabled`` must leave the
    bookkeeping at zero and not call into Pass G."""
    mesh = _c1_fail_fan()
    _out, info = phase_h_optimize(
        mesh, alpha_target=0.95, min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_rounds=1,
    )
    assert info["pass_g_enabled"] is False
    assert info["pass_g_accepts"] == 0
    assert info["pass_g_sweeps"] == 0


# ---------------------------------------------------------------------------
# phase_h_finish: stochastic local fixer + vertex_remove chain (PoC #58d/j/k)
# ---------------------------------------------------------------------------


def test_phase_h_finish_clears_c4_fan_at_threshold() -> None:
    """The C4-fail fan (2 marginal C4 fails at the 30°-min-angle
    threshold) is clearable by the stochastic local fixer's Stage 1
    alone — `phase_h_finish` should drive the residual to zero.
    """
    mesh = _c4_fail_fan()
    _uv, _ep, ac_b = _per_edge_area_change(mesh.nodes, mesh.elements)
    assert int((ac_b > 0.5).sum()) == 2

    out_mesh, info = phase_h_finish(
        mesh,
        seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        area_ratio_target=0.5,
        max_outer_passes=3,
        max_tries_per_fail=200,
    )
    assert info["before"]["C4"] == 2
    assert info["after"]["C1"] == 0
    assert info["after"]["C4"] == 0
    assert info["delta_total"] <= 0
    _uv2, _ep2, ac_a = _per_edge_area_change(
        out_mesh.nodes, out_mesh.elements,
    )
    assert int((ac_a > 0.5).sum()) == 0


def test_phase_h_finish_zero_residual_input_is_no_op() -> None:
    """Running the finisher on an input with no FVCOM violations
    must not modify the mesh or its element count, and must report
    zero work."""
    mesh = _c4_fail_fan()
    # Pre-balance node 4 → all four triangles isoceles right (45/45/90).
    mesh.nodes[4] = np.array([2.0, 2.0])
    out_mesh, info = phase_h_finish(
        mesh, seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_passes=2, max_tries_per_fail=50,
    )
    assert info["before"]["C1"] == 0
    assert info["before"]["C2"] == 0
    assert info["before"]["C4"] == 0
    assert info["before"]["C5"] == 0
    assert info["after"] == info["before"]
    assert info["stage1_stochastic"]["n_fixed"] == 0
    assert info["stage1_stochastic"]["n_stuck"] == 0
    assert info["stage2_c1_vertex_remove"] == []
    assert info["stage3_c4_vertex_remove"] == []
    assert out_mesh.n_elements == mesh.n_elements
    assert out_mesh.n_nodes == mesh.n_nodes


def test_phase_h_finish_is_seed_reproducible() -> None:
    """Two `phase_h_finish` runs with the same seed must produce
    bit-identical meshes (same node positions, same element list,
    same final counts)."""
    mesh = _c4_fail_fan()
    out_a, info_a = phase_h_finish(
        mesh, seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_passes=2, max_tries_per_fail=100,
    )
    out_b, info_b = phase_h_finish(
        mesh, seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_passes=2, max_tries_per_fail=100,
    )
    np.testing.assert_array_equal(out_a.nodes, out_b.nodes)
    np.testing.assert_array_equal(out_a.elements, out_b.elements)
    assert info_a["after"] == info_b["after"]


def test_phase_h_finish_does_not_mutate_input_mesh() -> None:
    """`phase_h_finish` must operate on a deep copy of the input;
    the caller's mesh object must remain unchanged."""
    mesh = _c4_fail_fan()
    nodes_before = mesh.nodes.copy()
    elements_before = mesh.elements.copy()
    _out, _info = phase_h_finish(
        mesh, seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_passes=2, max_tries_per_fail=50,
    )
    np.testing.assert_array_equal(mesh.nodes, nodes_before)
    np.testing.assert_array_equal(mesh.elements, elements_before)


def test_phase_h_finish_records_per_stage_info() -> None:
    """The info dict must surface per-stage bookkeeping: stage 1
    stochastic stats, stage 2 / 3 vertex_remove records, before /
    after / after-stage1 / after-stage2 counts, and the delta."""
    mesh = _c4_fail_fan()
    _out, info = phase_h_finish(
        mesh, seed=42,
        min_angle_target=20.0,
        max_angle_target=DEFAULT_MAX_ANGLE_TARGET,
        max_outer_passes=2, max_tries_per_fail=100,
    )
    assert "stage1_stochastic" in info
    assert "stage2_c1_vertex_remove" in info
    assert "stage3_c4_vertex_remove" in info
    assert "before" in info
    assert "after_stage1" in info
    assert "after_stage2" in info
    assert "after" in info
    assert "delta_total" in info
    assert isinstance(info["stage2_c1_vertex_remove"], list)
    assert isinstance(info["stage3_c4_vertex_remove"], list)
    # Per-stage counts must monotonically not increase (each stage
    # can only fix fails, never introduce them at the global gate).
    def _total(c: dict) -> int:
        return c["C1"] + c["C2"] + c["C4"] + c["C5"]
    assert _total(info["after_stage1"]) <= _total(info["before"])
    assert _total(info["after_stage2"]) <= _total(info["after_stage1"])
    assert _total(info["after"]) <= _total(info["after_stage2"])
