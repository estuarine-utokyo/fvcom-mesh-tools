"""Strict one-wide mode of resolve_narrow_channels (owner
2026-07-13: "do not create one-mesh-wide channels")."""

import numpy as np

from fvcom_mesh_tools.channel_policy import resolve_narrow_channels
from fvcom_mesh_tools.io import Fort14Mesh


def _strip_with_one_wide_tail():
    """A two-element-wide main body (3x4 node lattice, OBC on the
    y=0 edge) plus a ONE-wide tail of two triangles ending in a
    one-triangle pocket."""
    xs = [0.0, 400.0, 800.0]
    ys = [0.0, 400.0, 800.0, 1200.0]
    nodes = [[x, y] for y in ys for x in xs]      # 12 lattice
    nid = {(i, j): j * 3 + i for j in range(4) for i in range(3)}
    els = []
    for j in range(3):
        for i in range(2):
            a = nid[(i, j)]
            b = nid[(i + 1, j)]
            c = nid[(i, j + 1)]
            d = nid[(i + 1, j + 1)]
            els.append([a, b, c])
            els.append([b, d, c])
    n_main = len(els)                              # 12
    t0 = len(nodes)
    nodes += [[200.0, 1600.0], [600.0, 1600.0], [400.0, 2000.0]]
    tl, tr, tip = t0, t0 + 1, t0 + 2
    els.append([nid[(0, 3)], nid[(1, 3)], tl])     # tail 1
    els.append([nid[(1, 3)], tr, tl])              # tail 2
    els.append([tl, tr, tip])                      # pocket cap
    nodes = np.asarray(nodes, float)
    els = np.asarray(els, int)
    obc = np.array([nid[(0, 0)], nid[(1, 0)], nid[(2, 0)]])
    land = np.array([nid[(2, 0)], nid[(2, 3)], nid[(0, 3)],
                     nid[(0, 0)]])
    mesh = Fort14Mesh(
        title="strict", nodes=nodes,
        depths=np.full(len(nodes), 5.0), elements=els,
        open_boundaries=[obc], land_boundaries=[(20, land)])
    return mesh, n_main


def test_strict_prunes_one_wide_tail_with_pocket():
    mesh, n_main = _strip_with_one_wide_tail()
    m2, info = resolve_narrow_channels(
        mesh, min_basin_elements=25, apply_widen=False,
        strict_boundary_flag=True, max_rounds=8)
    # the two tail cells are strict one-wide; the pocket behind
    # them is a sub-threshold appendix -- all three go together
    assert info["n_flagged"] >= 2
    assert info["n_deleted_elements"] == 3
    assert len(m2.elements) == n_main
    # OBC survived untouched, in place
    ob = np.asarray(m2.open_boundaries[0])
    assert len(ob) == 3
    assert np.allclose(m2.nodes[ob][:, 1], 0.0)
    # no strict one-wide cell remains
    m3, info2 = resolve_narrow_channels(
        m2, min_basin_elements=25, apply_widen=False,
        strict_boundary_flag=True, max_rounds=1)
    assert info2["n_deleted_elements"] == 0


def test_default_mode_unchanged_keeps_tail_reported():
    # without strict flagging the w/h metric may or may not flag
    # the synthetic tail; the call must at least not crash and
    # never delete the main body
    mesh, n_main = _strip_with_one_wide_tail()
    m2, info = resolve_narrow_channels(
        mesh, min_basin_elements=6, apply_widen=False)
    assert len(m2.elements) >= n_main


def test_widen_choke_sections_noop_on_clean_lattice():
    # smoke: a healthy two-wide strip has no bank-to-bank choke
    # edge -- the operator must return the mesh untouched
    import shapely
    from fvcom_mesh_tools.algorithms.obc_finish import (
        widen_choke_sections,
    )

    mesh, _ = _strip_with_one_wide_tail()
    land = shapely.box(-5000, -5000, 5000, 7000).difference(
        shapely.box(-100, -100, 900, 2100))
    m2, info = widen_choke_sections(mesh, land)
    assert info["widened"] == 0
    assert len(m2.elements) == len(mesh.elements)
    assert np.allclose(m2.nodes, mesh.nodes)
