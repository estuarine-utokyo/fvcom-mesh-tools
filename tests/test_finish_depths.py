"""fmesh-finish-depths: the TB-FVCOM depth product, with its three numbers as options."""

from __future__ import annotations

import numpy as np
import pytest

from fvcom_mesh_tools.cli.finish_depths import main, variant_tag
from fvcom_mesh_tools.io.fort14 import Fort14Mesh
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_dep


def _strip(n=8):
    xs = np.arange(n, dtype=float) * 100.0
    nodes = np.column_stack([np.repeat(xs, 2), np.tile([0.0, 100.0], n)])
    tri = []
    for i in range(n - 1):
        a = 2 * i
        tri += [[a, a + 3, a + 1], [a, a + 2, a + 3]]      # counter-clockwise
    return nodes, np.asarray(tri, dtype=np.int64)


def _case(tmp_path, depths, name="case"):
    nodes, tri = _strip(len(depths) // 2)
    mesh = Fort14Mesh(title=name, nodes=nodes, depths=np.asarray(depths, float),
                      elements=tri, open_boundaries=[np.array([0, 1])],
                      land_boundaries=[])
    export_fvcom_case(mesh, tmp_path, name, twodm=False, obc_depth_control=False)
    return tmp_path / f"{name}_grd.dat", tri


@pytest.mark.parametrize("hmin,hmax,r,tag", [
    (3, 300, 0.2, "min3m_rfac0p2_cap300"),
    (2.5, 50, 0.15, "min2p5m_rfac0p15_cap050"),
    (3, 300, None, "min3m_cap300"),
])
def test_variant_tag_follows_the_tb_fvcom_spelling(hmin, hmax, r, tag):
    assert variant_tag(hmin, hmax, r) == tag


def test_a_case_is_floored_smoothed_and_capped_everywhere(tmp_path):
    h = [0.5, 1.0, 4.0, 4.5, 40.0, 42.0, 45.0, 44.0, 400.0, 420.0, 30.0, 31.0,
         20.0, 22.0, 18.0, 19.0]
    grd, tri = _case(tmp_path, h)
    assert main([str(grd), "--hmin", "3", "--hmax", "300", "--rfactor", "0.2"]) == 0
    out = tmp_path / "case_dep_min3m_rfac0p2_cap300.dat"
    assert out.exists() and (tmp_path / "case_dep_min3m_rfac0p2_cap300.json").exists()
    _, d = read_dep(out)
    assert d.min() >= 3.0 - 1e-6 and d.max() <= 300.0 + 1e-6
    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]),
                          axis=1), axis=0)
    r = np.abs(d[e[:, 0]] - d[e[:, 1]]) / (d[e[:, 0]] + d[e[:, 1]])
    assert r.max() <= 0.2 + 1e-5


def test_rfactor_zero_is_floor_and_cap_only(tmp_path):
    h = [0.5, 1.0, 4.0, 4.5, 40.0, 42.0, 400.0, 44.0]
    grd, _ = _case(tmp_path, h)
    assert main([str(grd), "--hmin", "3", "--hmax", "300", "--rfactor", "0"]) == 0
    _, d = read_dep(tmp_path / "case_dep_min3m_cap300.dat")
    assert np.allclose(d, np.clip(h, 3.0, 300.0))


def test_a_starting_depth_file_can_be_named(tmp_path):
    grd, _ = _case(tmp_path, [5.0] * 8)
    alt = tmp_path / "case_dep_raw.dat"
    alt.write_text((tmp_path / "case_dep.dat").read_text().replace(" 5.000000", " 1.000000"))
    assert main([str(grd), "--dep", "case_dep_raw.dat", "--hmin", "2",
                 "--rfactor", "0", "--tag", "mine"]) == 0
    _, d = read_dep(tmp_path / "case_dep_mine.dat")
    assert np.allclose(d, 2.0)


def test_a_refinement_keeps_its_frozen_depths_unless_told(tmp_path):
    h = [1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 8.0, 8.0]
    fv = tmp_path / "ref" / "fvcom"
    fv.mkdir(parents=True)
    _case(fv, h, name="p")
    # base nodes 0..3 kept, 4..7 new
    np.save(tmp_path / "ref" / "node_map.npy", np.array([0, 1, 2, 3]))
    assert main([str(tmp_path / "ref"), "--hmin", "3", "--rfactor", "0"]) == 0
    _, d = read_dep(tmp_path / "ref" / "fvcom_finished" / "p_dep_min3m_cap300.dat")
    assert np.allclose(d[:4], 1.0) and np.allclose(d[4:6], 3.0)
    assert (tmp_path / "ref" / "fvcom_finished" / "p_grd.dat").exists()
    assert main([str(tmp_path / "ref"), "--hmin", "3", "--rfactor", "0",
                 "--whole-mesh", "--outdir", str(tmp_path / "all")]) == 0
    _, d2 = read_dep(tmp_path / "all" / "p_dep_min3m_cap300.dat")
    assert np.allclose(d2[:6], 3.0)


def test_a_wrong_source_is_refused(tmp_path):
    assert main([str(tmp_path), "--hmin", "3"]) == 2
    assert main([str(tmp_path / "nothing.dat"), "--hmin", "3"]) == 2


def test_equal_smoothing_keeps_a_pair_s_sum_and_honours_frozen_nodes():
    from fvcom_mesh_tools.refine import smooth_rfactor_equal

    tri = np.array([[0, 1, 2]])
    h = np.array([10.0, 40.0, 40.0])
    out, rep = smooth_rfactor_equal(tri, h, np.ones(3, bool), 0.2, depth_min=3.0)
    assert rep["converged"] and rep["max_r"] <= 0.2 + 1e-6
    # equal and opposite: the channel spreads into its bank, it is not filled
    assert out.sum() == pytest.approx(h.sum(), rel=1e-6)
    frozen = np.array([True, False, False])
    out2, _ = smooth_rfactor_equal(tri, h, ~frozen, 0.2, depth_min=3.0)
    assert out2[0] == 10.0
