"""Tests for ``fvcom_mesh_tools.mesh_compose``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_compose import combine
from fvcom_mesh_tools.mesh_compose.disjoint import combine_disjoint


def _square_mesh(
    title: str = "square",
    *,
    x0: float = 0.0,
    y0: float = 0.0,
    nx: int = 3,
    ny: int = 3,
    spacing: float = 1.0,
    depth_value: float = 5.0,
) -> Fort14Mesh:
    """Tiny structured triangular mesh for unit tests.

    Lays out an ``nx`` x ``ny`` grid of nodes, divides each cell into
    two CCW triangles, marks the four corners of the rectangle as
    open boundary points and the rest of the perimeter as a single
    land boundary. Pure numpy, no GIS deps required.
    """
    xs = x0 + spacing * np.arange(nx)
    ys = y0 + spacing * np.arange(ny)
    XX, YY = np.meshgrid(xs, ys)
    nodes = np.column_stack([XX.ravel(), YY.ravel()])
    depths = np.full(nodes.shape[0], depth_value)

    elements = []

    def idx(i: int, j: int) -> int:
        return j * nx + i

    for j in range(ny - 1):
        for i in range(nx - 1):
            elements.append([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)])
            elements.append([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)])
    elements = np.asarray(elements, dtype=np.int64)

    # Trivial open / land boundary structure: south edge open,
    # remaining perimeter as one land segment (non-closed).
    south = np.array([idx(i, 0) for i in range(nx)], dtype=np.int64)
    east = [idx(nx - 1, j) for j in range(1, ny)]
    north = [idx(i, ny - 1) for i in range(nx - 2, -1, -1)]
    west = [idx(0, j) for j in range(ny - 2, 0, -1)]
    land = np.asarray(east + north + west, dtype=np.int64)

    return Fort14Mesh(
        title=title,
        nodes=nodes,
        depths=depths,
        elements=elements,
        open_boundaries=[south],
        land_boundaries=[(20, land)],
    )


def test_disjoint_node_count_adds_up() -> None:
    a = _square_mesh("A", x0=0.0, y0=0.0, depth_value=5.0)
    b = _square_mesh("B", x0=10.0, y0=10.0, depth_value=12.0)
    out = combine_disjoint([a, b])
    assert out.n_nodes == a.n_nodes + b.n_nodes
    assert out.n_elements == a.n_elements + b.n_elements


def test_disjoint_preserves_boundaries_with_offsets() -> None:
    a = _square_mesh("A", nx=3, ny=3)
    b = _square_mesh("B", nx=4, ny=4, x0=20.0)
    out = combine([
        "disjoint" if False else "disjoint",  # via dispatcher
    ][0], [a, b])

    # Each input contributes 1 open + 1 land segment.
    assert len(out.open_boundaries) == 2
    assert len(out.land_boundaries) == 2

    # First mesh's segments are unshifted; second mesh's segments are
    # shifted by a.n_nodes.
    assert np.array_equal(out.open_boundaries[0], a.open_boundaries[0])
    assert np.array_equal(
        out.open_boundaries[1], b.open_boundaries[0] + a.n_nodes
    )
    assert out.land_boundaries[0][0] == 20
    assert np.array_equal(
        out.land_boundaries[1][1], b.land_boundaries[0][1] + a.n_nodes
    )


def test_disjoint_depths_concat() -> None:
    a = _square_mesh("A", depth_value=3.0)
    b = _square_mesh("B", depth_value=8.5)
    out = combine_disjoint([a, b])
    assert np.allclose(out.depths[: a.n_nodes], 3.0)
    assert np.allclose(out.depths[a.n_nodes:], 8.5)


def test_disjoint_round_trip_through_fort14(tmp_path: Path) -> None:
    a = _square_mesh("A", x0=0.0, depth_value=2.0)
    b = _square_mesh("B", x0=10.0, depth_value=7.0)
    out = combine_disjoint([a, b], title="A+B")
    f14 = tmp_path / "combined.14"
    write_fort14(out, f14)
    rt = read_fort14(f14)
    assert rt.n_nodes == out.n_nodes
    assert rt.n_elements == out.n_elements
    assert np.allclose(rt.nodes, out.nodes)
    assert len(rt.open_boundaries) == 2
    assert len(rt.land_boundaries) == 2


def test_combine_dispatcher_unknown_strategy() -> None:
    a = _square_mesh("A")
    b = _square_mesh("B")
    with pytest.raises(ValueError, match="unknown strategy"):
        combine("nope", [a, b])


def test_combine_dispatcher_too_few_inputs() -> None:
    with pytest.raises(ValueError, match="at least 2 meshes"):
        combine("disjoint", [_square_mesh("A")])


def test_disjoint_no_overlap_in_indices() -> None:
    """Sanity: no element in B should index into A's node range after combine."""
    a = _square_mesh("A", nx=3, ny=3)
    b = _square_mesh("B", nx=5, ny=5)
    out = combine_disjoint([a, b])
    a_elems = out.elements[: a.n_elements]
    b_elems = out.elements[a.n_elements:]
    assert int(a_elems.max()) < a.n_nodes
    assert int(b_elems.min()) >= a.n_nodes
    assert int(b_elems.max()) < a.n_nodes + b.n_nodes
