import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def _grid3():
    n = 3
    nodes = np.array([[i * 1000.0, j * 1000.0]
                      for j in range(n) for i in range(n)])
    elements = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b = j * n + i, j * n + i + 1
            c, d = (j + 1) * n + i + 1, (j + 1) * n + i
            elements.append([a, b, c])
            elements.append([a, c, d])
    return Fort14Mesh(
        title="site",
        nodes=nodes,
        depths=np.full(n * n, 5.0),
        elements=np.asarray(elements),
        open_boundaries=[np.array([5, 8])],
        land_boundaries=[(20, np.array([8, 7, 6, 3, 0, 1, 2, 5]))],
    )


def test_insert_node_on_line_splits_wide_boundary_edge():
    """The operator's real use: a boundary edge ~2x the local scale
    (carrier triangle wide and squat) splits cleanly within gates.
    Mid-splitting a 45-90-45 carrier is correctly rejected by the
    quality gate (it halves the base angles below 30 deg), so the
    happy path needs the wide-edge geometry."""
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms import signed_areas
    from fvcom_mesh_tools.algorithms.site_edits import insert_node_on_line

    mesh = Fort14Mesh(
        title="wide",
        nodes=np.array([[0.0, 0.0], [2000.0, 0.0], [1000.0, 900.0]]),
        depths=np.array([4.0, 6.0, 8.0]),
        elements=np.array([[0, 1, 2]]),
        open_boundaries=[np.array([2])],
        land_boundaries=[(20, np.array([1, 0]))],
    )
    coast = LineString([(-500.0, 50.0), (2500.0, 50.0)])
    result = insert_node_on_line(mesh, 0, 1, coast)
    assert result is not None
    out, info = result
    n_new = info["new_node"]
    assert n_new == 3
    assert np.allclose(out.nodes[n_new], [1000.0, 50.0])
    assert out.n_elements == 2
    assert (signed_areas(out) > 0).all()
    assert info["boundary_updated"]
    seg = out.land_boundaries[0][1]
    assert list(seg) == [1, 3, 0]
    assert out.depths[n_new] == 5.0


def test_insert_node_on_line_rejects_sharp_split():
    """Mid-split of a right-angle carrier halves the 45 deg base
    angles below C1 -> must be refused."""
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms.site_edits import insert_node_on_line

    mesh = _grid3()
    coast = LineString([(-500.0, 80.0), (2500.0, 80.0)])
    assert insert_node_on_line(mesh, 0, 1, coast) is None


def test_insert_node_on_line_rejects_non_boundary_edge():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms.site_edits import insert_node_on_line

    mesh = _grid3()
    # (0,4) is an interior diagonal shared by two elements.
    coast = LineString([(-500.0, 80.0), (2500.0, 80.0)])
    assert insert_node_on_line(mesh, 0, 4, coast) is None


def test_insert_node_on_line_rejects_quality_loss():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms.site_edits import insert_node_on_line

    mesh = _grid3()
    # A line far above the edge midpoint would create a spike node
    # inside the carrier triangle's territory -> flips one half.
    bad = LineString([(-500.0, 900.0), (2500.0, 900.0)])
    assert insert_node_on_line(mesh, 0, 1, bad) is None


def _strip_mesh():
    """Two-row strip: y=0 (boundary toward the coast) and y=300."""
    xs = np.arange(0.0, 3300.0, 300.0)
    bottom = np.column_stack([xs, np.zeros_like(xs)])
    top = np.column_stack([xs, np.full_like(xs, 300.0)])
    nodes = np.vstack([bottom, top])
    n = len(xs)
    elements = []
    for i in range(n - 1):
        elements.append([i, i + 1, n + i + 1])
        elements.append([i, n + i + 1, n + i])
    return Fort14Mesh(
        title="strip",
        nodes=nodes,
        depths=np.full(2 * n, 3.0),
        elements=np.asarray(elements),
        open_boundaries=[np.arange(n, 2 * n)],
        land_boundaries=[(20, np.arange(n - 1, -1, -1))],
    )


def test_extrude_boundary_strip_fills_retreat_gap():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms import signed_areas
    from fvcom_mesh_tools.algorithms.site_edits import (
        extrude_boundary_strip,
    )

    mesh = _strip_mesh()
    n0 = mesh.n_nodes
    e0 = mesh.n_elements
    coast = LineString([(-500.0, -240.0), (3800.0, -240.0)])
    out, info = extrude_boundary_strip(mesh, [coast])
    assert info["n_strips"] == 1
    assert info["n_new_nodes"] == 11
    assert info["n_new_elements"] == 20
    assert out.n_nodes == n0 + 11
    assert out.n_elements == e0 + 20
    assert np.allclose(out.nodes[n0:, 1], -240.0)
    assert (signed_areas(out) > 0).all()
    assert np.allclose(out.depths[n0:], 3.0)


def test_extrude_boundary_strip_skips_far_and_near_lines():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms.site_edits import (
        extrude_boundary_strip,
    )

    mesh = _strip_mesh()
    far = LineString([(-500.0, -600.0), (3800.0, -600.0)])
    out, info = extrude_boundary_strip(mesh, [far])
    assert info["n_strips"] == 0
    near = LineString([(-500.0, -30.0), (3800.0, -30.0)])
    out, info = extrude_boundary_strip(mesh, [near])
    assert info["n_strips"] == 0


def test_collapse_edge_removes_boundary_needle():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms import signed_areas
    from fvcom_mesh_tools.algorithms.site_edits import collapse_edge

    # Needle ear: boundary nodes 0,1 are 20 m apart on a ~600 m coast
    # edge row; apex 2 is 600 m away -> sliver with 2 boundary edges.
    nodes = np.array([
        [0.0, 0.0], [20.0, 2.0], [300.0, 600.0],
        [-600.0, 500.0], [700.0, 520.0],
    ])
    mesh = Fort14Mesh(
        title="needle",
        nodes=nodes,
        depths=np.full(5, 4.0),
        elements=np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2]]),
        open_boundaries=[np.array([3])],
        land_boundaries=[(20, np.array([3, 0, 1, 4]))],
    )
    coast = LineString([(-800.0, -10.0), (900.0, 30.0)])
    res = collapse_edge(mesh, 0, 1, lines=[coast])
    assert res is not None
    out, info = res
    assert info["n_elements_removed"] == 1
    assert out.n_elements == 2
    assert (signed_areas(out) > 0).all()
    # Survivor sits ON the coast line.
    seg = out.land_boundaries[0][1]
    assert list(seg) == [3, 0, 4]
    p = out.nodes[0]
    # distance from the line y = (x+800)*40/1700 - 10 is ~0
    from shapely.geometry import Point
    assert coast.distance(Point(p[0], p[1])) < 1e-6


def test_collapse_edge_refuses_protected_and_no_gain():
    from fvcom_mesh_tools.algorithms.site_edits import collapse_edge

    mesh = _grid3()
    # OBC-protected node.
    assert collapse_edge(mesh, 5, 8, protected=[5, 8]) is None
    # Healthy grid: collapsing any edge makes things worse -> None.
    assert collapse_edge(mesh, 0, 1) is None


def test_split_edge_pair_grades_long_chord():
    from fvcom_mesh_tools.algorithms import signed_areas
    from fvcom_mesh_tools.algorithms.site_edits import split_edge_pair

    # Two big triangles over a long shared edge, flanked by small
    # ones -> C4 jump; splitting the shared edge must reduce it.
    nodes = np.array([
        [0.0, 0.0], [2000.0, 0.0],       # long edge 0-1
        [1000.0, 1500.0], [1000.0, -1500.0],
        [-400.0, 400.0], [-400.0, -400.0],
    ])
    mesh = Fort14Mesh(
        title="chord",
        nodes=nodes,
        depths=np.full(6, 10.0),
        elements=np.array([
            [0, 1, 2], [0, 3, 1], [0, 2, 4], [0, 4, 5], [0, 5, 3],
        ]),
        open_boundaries=[np.array([2])],
        land_boundaries=[(20, np.array([2, 1, 3]))],
    )
    res = split_edge_pair(mesh, 0, 1)
    assert res is not None
    out, info = res
    assert out.n_elements == 7
    assert out.n_nodes == 7
    assert (signed_areas(out) > 0).all()
    b, a = info["violations"]
    assert a < b
