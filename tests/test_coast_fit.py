"""Tests for the coastline-fitting pass."""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Polygon

from fvcom_mesh_tools.coast_fit import (
    boundary_edges,
    boundary_offsets,
    fit_boundary_to_coast,
)


def _strip_mesh(nx: int = 9, ny: int = 3, h: float = 100.0):
    """A regular right-triangle strip over ``[0, nx*h] x [0, ny*h]``."""
    xs, ys = np.meshgrid(np.arange(nx + 1) * h, np.arange(ny + 1) * h, indexing="ij")
    nodes = np.column_stack([xs.ravel(), ys.ravel()])

    def nid(i, j):
        return i * (ny + 1) + j

    tri = []
    for i in range(nx):
        for j in range(ny):
            tri.append([nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)])
            tri.append([nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)])
    return nodes, np.asarray(tri, dtype=np.int64)


def _land_above(y: float, nx: int = 9, h: float = 100.0) -> Polygon:
    """Land occupying everything above ``y`` over the strip's span."""
    x0, x1 = -h, nx * h + h
    return Polygon([(x0, y), (x1, y), (x1, y + 10 * h), (x0, y + 10 * h)])


def test_boundary_edges_counts_the_rim_only():
    nodes, tri = _strip_mesh(nx=4, ny=2)
    be = boundary_edges(tri)
    # A 4x2 right-triangle strip has 2*(4+2) = 12 boundary edges.
    assert len(be) == 12
    assert be.shape[1] == 2


def test_offsets_are_signed_positive_in_water():
    nodes, tri = _strip_mesh()
    # Coast at y = 340: the top row (y = 300) is 40 m short of it, in water.
    land = _land_above(340.0)
    ids, off = boundary_offsets(nodes, tri, land)
    top = nodes[ids, 1] == 300.0
    assert np.allclose(off[top], 40.0)
    assert (off > 0).all()


def test_offsets_are_negative_when_the_node_is_on_land():
    nodes, tri = _strip_mesh()
    # Coast at y = 250: the top row has crossed 50 m onto land.
    land = _land_above(250.0)
    ids, off = boundary_offsets(nodes, tri, land)
    top = nodes[ids, 1] == 300.0
    assert np.allclose(off[top], -50.0)


def test_fit_pulls_the_boundary_onto_the_coast():
    nodes, tri = _strip_mesh()
    land = _land_above(340.0)
    res = fit_boundary_to_coast(nodes, tri, land, sweeps=8)
    assert res.n_moved > 0
    # The coast only bounds the top row; the other three sides are open water
    # in this fixture, so judge the fit on the row that faces it.
    ids, after = boundary_offsets(res.nodes, tri, land)
    top = nodes[ids, 1] == 300.0
    assert after[top].max() == pytest.approx(0.0, abs=1e-3)  # relax=0.8 ** 8 sweeps
    # The far side of the strip is open water 340 m from the coast and clamped
    # by the move budget, so the aggregate only has to improve, not vanish.
    assert res.after["median"] < res.before["median"]
    # The elements are never touched, only the coordinates.
    assert res.nodes.shape == nodes.shape


def test_fit_never_moves_a_fixed_node():
    nodes, tri = _strip_mesh()
    land = _land_above(340.0)
    ids, _ = boundary_offsets(nodes, tri, land)
    frozen = ids[nodes[ids, 1] == 300.0]
    res = fit_boundary_to_coast(nodes, tri, land, fixed=frozen, sweeps=8)
    assert np.allclose(res.nodes[frozen], nodes[frozen])


def test_fit_respects_the_move_budget():
    nodes, tri = _strip_mesh()
    land = _land_above(1000.0)  # far away: every node wants a huge move
    res = fit_boundary_to_coast(nodes, tri, land, sweeps=8, max_move_frac=0.25)
    # Local boundary edge length is 100 m, so no node may travel past 25 m.
    assert res.max_move_m <= 25.0 + 1e-6


def test_fit_leaves_no_inverted_element():
    nodes, tri = _strip_mesh()
    land = _land_above(320.0)
    res = fit_boundary_to_coast(nodes, tri, land, sweeps=8)
    p = res.nodes[tri]
    a = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
         - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    assert (a > 0).all() or (a < 0).all()


def test_dt_floor_blocks_moves_that_would_cost_the_time_step():
    nodes, tri = _strip_mesh()
    land = _land_above(305.0)
    depths = np.full(len(nodes), 10.0)
    loose = fit_boundary_to_coast(nodes, tri, land, sweeps=8,
                                  depths=depths, dt_floor_s=0.0)
    tight = fit_boundary_to_coast(nodes, tri, land, sweeps=8, depths=depths)
    # The default floor is the mesh's own dt, so it can never come out lower.
    assert tight.dt_after_s >= tight.dt_before_s - 1e-9
    assert tight.n_rejected >= loose.n_rejected


def test_fit_on_a_mesh_with_no_land_boundary_is_a_no_op():
    nodes, tri = _strip_mesh(nx=3, ny=2)
    land = _land_above(340.0)
    ids, _ = boundary_offsets(nodes, tri, land)
    res = fit_boundary_to_coast(nodes, tri, land, fixed=ids, sweeps=4)
    assert res.n_moved == 0
    assert np.allclose(res.nodes, nodes)
