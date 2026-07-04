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
