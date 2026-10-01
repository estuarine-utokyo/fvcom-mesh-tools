

def test_c4_flip_refuses_a_concave_quadrilateral():
    """Review of the extend tools, round 5 F2: the flip folded the mesh and
    reported it fixed."""
    import numpy as np

    from fvcom_mesh_tools.algorithms.obc_finish import flip_c4_edges
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh

    nodes = np.array([[3951, 982], [2396, 1743], [926, 262], [2471, 635]], float)
    m = Fort14Mesh("t", nodes, np.ones(4), np.array([[3, 2, 0], [1, 3, 0]]), [], [])
    out = flip_c4_edges(m)
    assert out["fixed"] == [] and m.elements.tolist() == [[3, 2, 0], [1, 3, 0]]


def test_flip_pair_on_a_convex_quadrilateral():
    import numpy as np

    from fvcom_mesh_tools.algorithms.obc_finish import _area, _flip_pair

    nodes = np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 1.0], [1.0, -1.0]])
    t1, t2 = _flip_pair(nodes, 0, 1, 2, 3)
    assert {frozenset(t1), frozenset(t2)} == {frozenset((2, 3, 0)), frozenset((2, 3, 1))}
    assert _area(*nodes[t1]) > 0 and _area(*nodes[t2]) > 0
    nodes[3] = [0.5, 0.2]                       # now on the same side as node 2
    assert _flip_pair(nodes, 0, 1, 2, 3) is None
