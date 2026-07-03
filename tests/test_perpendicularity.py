from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
    unique_edges,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

REFERENCE_MESH = (
    Path(__file__).parent.parent
    / "data"
    / "mesh"
    / "reference"
    / "tokyo_bay"
    / "tb_futtsu20220311.14"
)


def _mixed_parent_mesh() -> Fort14Mesh:
    """6-node strip with one movable single-parent interior (node 2, only
    parent 4) and one movable multi-parent interior (node 3, parents 4 and
    5). Top edge (4,5) is the open boundary; the diagonal triangulation
    isolates node 2 from open node 5."""
    nodes = np.array(
        [
            [0.0, 0.0],   # 0 land
            [1.0, 0.0],   # 1 land
            [0.2, 0.5],   # 2 interior, single-parent
            [0.5, 0.5],   # 3 interior, multi-parent
            [0.0, 1.0],   # 4 open
            [1.0, 1.0],   # 5 open
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            [0, 1, 3],
            [0, 3, 2],
            [2, 3, 4],
            [3, 5, 4],
            [1, 5, 3],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="mixed",
        nodes=nodes,
        depths=np.zeros(6),
        elements=elements,
        open_boundaries=[np.array([4, 5], dtype=np.int64)],
        land_boundaries=[(0, np.array([0, 1], dtype=np.int64))],
    )


def test_unique_edges_count() -> None:
    elements = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    edges = unique_edges(elements)
    assert edges.shape == (5, 2)
    pairs = {tuple(e) for e in edges.tolist()}
    assert pairs == {(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)}


def test_signed_areas_positive_for_ccw() -> None:
    mesh = _mixed_parent_mesh()
    a = signed_areas(mesh)
    assert (a > 0).all()


def test_align_perfects_single_parent_edge() -> None:
    """The single-parent movable interior edge must be driven to 0-deg
    deviation; the multi-parent edge only improves."""
    before = _mixed_parent_mesh()
    perp_before = open_bdy_perpendicularity(before)
    after, info = align_open_boundary_first_ring(before, alpha=1.0, n_iters=1)
    perp_after = open_bdy_perpendicularity(after)

    movable = info["movable_first_ring_by_parent_count"]
    assert movable.get(1, 0) == 1   # node 2
    assert movable.get(2, 0) == 1   # node 3
    assert perp_after.min() <= 1e-9
    assert perp_after.mean() < perp_before.mean()


def test_align_keeps_boundary_nodes_fixed() -> None:
    before = _mixed_parent_mesh()
    open_ids = before.open_boundaries[0]
    land_ids = before.land_boundaries[0][1]
    after, _ = align_open_boundary_first_ring(before, alpha=1.0, n_iters=1)
    np.testing.assert_array_equal(after.nodes[open_ids], before.nodes[open_ids])
    np.testing.assert_array_equal(after.nodes[land_ids], before.nodes[land_ids])


def test_damped_multi_iter_converges_to_single_alpha_one() -> None:
    """alpha=1 / n_iters=1 reaches the fixed point in one shot. Damped
    iteration (alpha<1, n_iters large) must converge to the same result
    because the targets are precomputed from the original geometry."""
    before = _mixed_parent_mesh()
    one_shot, _ = align_open_boundary_first_ring(before, alpha=1.0, n_iters=1)
    damped, _ = align_open_boundary_first_ring(before, alpha=0.5, n_iters=80)
    # alpha=0.5 / 80 iters: residual = 0.5**80 ~ 1e-24, far below float noise.
    np.testing.assert_allclose(damped.nodes, one_shot.nodes, atol=1e-12)


def test_align_no_op_when_no_movable_first_ring() -> None:
    """If every interior endpoint of every incident edge happens to be a
    boundary node (open or land), nothing is movable and the algorithm is
    a no-op."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        dtype=np.float64,
    )
    elements = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="all-boundary",
        nodes=nodes,
        depths=np.zeros(4),
        elements=elements,
        open_boundaries=[np.array([3, 2], dtype=np.int64)],
        land_boundaries=[(0, np.array([0, 1], dtype=np.int64))],
    )
    after, info = align_open_boundary_first_ring(mesh, alpha=1.0, n_iters=1)
    np.testing.assert_array_equal(after.nodes, mesh.nodes)
    assert info["moved"] == 0


@pytest.mark.skipif(not REFERENCE_MESH.exists(), reason="reference mesh symlink not in place")
def test_align_improves_reference_mesh() -> None:
    """End-to-end check on the live Tokyo Bay mesh: the algorithm must
    reduce mean perpendicularity deviation and not flip any triangle."""
    before = read_fort14(REFERENCE_MESH)
    perp_before = open_bdy_perpendicularity(before)
    after, info = align_open_boundary_first_ring(before, alpha=1.0, n_iters=1)
    perp_after = open_bdy_perpendicularity(after)

    assert info["moved"] > 0
    assert perp_after.mean() < perp_before.mean()
    assert (signed_areas(after) > 0).all()


# ---------------------------------------------------------------------------
# align_open_boundary_local (promoted from PoC #59e)
# ---------------------------------------------------------------------------


def _offset_grid_mesh(offset: float = 400.0):
    """4x4-node right-triangle grid; OBC = bottom-row nodes [1, 2].

    The first-ring interior nodes 5 and 6 are shifted +``offset`` in x,
    so each OBC node's best interior edge deviates
    ``atan(offset/1000)`` (~21.8 deg at 400) from perpendicular while
    the surrounding triangles keep enough C1 headroom for a local fix.
    (A uniformly SHEARED grid is the wrong fixture: every triangle
    sits ~2.7 deg above the 30 deg gate, so any perpendicularizing
    move is infeasible — the fixer correctly refuses it.)
    """
    n = 4
    sp = 1000.0
    nodes = np.array(
        [[i * sp, j * sp] for j in range(n) for i in range(n)],
        dtype=np.float64,
    )
    nodes[5, 0] += offset
    nodes[6, 0] += offset
    elements = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b = j * n + i, j * n + i + 1
            c, d = (j + 1) * n + i + 1, (j + 1) * n + i
            elements.append([a, b, c])
            elements.append([a, c, d])
    land = [2, 3, 7, 11, 15, 14, 13, 12, 8, 4, 0, 1]
    return Fort14Mesh(
        title="perp-local-test",
        nodes=nodes,
        depths=np.full(n * n, 5.0),
        elements=np.asarray(elements),
        open_boundaries=[np.array([1, 2])],
        land_boundaries=[(20, np.asarray(land))],
    )


def test_align_open_boundary_local_fixes_without_breaking_quality():
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.qa import run_qa

    mesh = _offset_grid_mesh()
    before = run_qa(mesh, channel_check=False)
    perp_before = [c for c in before.checks
                   if c.check_id == "obc_perpendicularity"][0]
    assert not perp_before.passed

    fixed, info = align_open_boundary_local(mesh, seed=7)
    assert info["remaining"] == []
    assert info["accepted_total"] >= 1

    after = run_qa(fixed, channel_check=False)
    by_id = {c.check_id: c for c in after.checks}
    assert by_id["obc_perpendicularity"].passed
    for cid in ("c1_min_angle", "c2_max_angle", "c4_area_change",
                "c5_valence", "ccw_all_elements"):
        assert by_id[cid].passed, cid
    # Input mesh untouched (pure function).
    assert np.array_equal(mesh.nodes, _offset_grid_mesh().nodes)


def test_align_open_boundary_local_noop_when_already_perpendicular():
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local

    mesh = _offset_grid_mesh(offset=0.0)
    fixed, info = align_open_boundary_local(mesh, seed=7)
    assert info["accepted_total"] == 0
    assert info["passes"][0]["violations"] == 0
    assert np.array_equal(fixed.nodes, mesh.nodes)
