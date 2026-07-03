"""Tests for the FVCOM native writers (io/fvcom_native.py, fmesh-export-fvcom).

The reference parsers below mimic the FVCOM 5.1 cold-start readers as
documented in ``docs/fvcom_source_constraints.md``: keyword headers,
free-format rows, and — for ``_grd.dat`` — the ``NVG = (N1, N3, N2)``
column swap followed by the clockwise check (a valid file is CCW, so
the swapped connectivity must be CW for every element).
"""

from __future__ import annotations

import numpy as np
import pytest

from fvcom_mesh_tools.cli import exportfvcom
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
from fvcom_mesh_tools.io.fvcom_native import (
    export_fvcom_case,
    write_2dm,
    write_cor,
    write_dep,
    write_grd,
    write_obc,
    write_spg,
)

N = 3  # 3x3 nodes, 4 squares, 8 CCW triangles


def _nid(i: int, j: int) -> int:
    return j * N + i


def _mesh() -> Fort14Mesh:
    nodes = np.array(
        [[i * 1000.0, j * 1000.0] for j in range(N) for i in range(N)],
    )
    elements = []
    for j in range(N - 1):
        for i in range(N - 1):
            a, b = _nid(i, j), _nid(i + 1, j)
            c, d = _nid(i + 1, j + 1), _nid(i, j + 1)
            elements.append([a, b, c])
            elements.append([a, c, d])
    return Fort14Mesh(
        title="native-test",
        nodes=nodes,
        depths=np.arange(2.0, 2.0 + N * N),
        elements=np.asarray(elements, dtype=np.int64),
        open_boundaries=[np.array([_nid(0, 1), _nid(0, 2)])],
        land_boundaries=[(20, np.array(
            [_nid(0, 2), _nid(1, 2), _nid(2, 2), _nid(2, 1), _nid(2, 0),
             _nid(1, 0), _nid(0, 0), _nid(0, 1)],
        ))],
    )


def _data_rows(path):
    return [ln.split() for ln in path.read_text().splitlines() if "=" not in ln]


def _header_int(path, keyword):
    for ln in path.read_text().splitlines():
        if keyword in ln:
            return int(ln.split("=")[1])
    raise AssertionError(f"header keyword {keyword!r} not found")


# ---------------------------------------------------------------------------
# _grd.dat
# ---------------------------------------------------------------------------


def test_write_grd_matches_fvcom_reader(tmp_path):
    mesh = _mesh()
    p = write_grd(mesh, tmp_path / "case_grd.dat")
    assert _header_int(p, "Node Number") == mesh.n_nodes
    assert _header_int(p, "Cell Number") == mesh.n_elements
    rows = _data_rows(p)
    conn = np.array([[int(v) for v in r] for r in rows[: mesh.n_elements]])
    coords = np.array([[float(v) for v in r] for r in rows[mesh.n_elements:]])
    assert conn.shape == (mesh.n_elements, 4)
    assert coords.shape == (mesh.n_nodes, 3)
    assert np.array_equal(conn[:, 1:] - 1, mesh.elements)
    assert np.allclose(coords[:, 1:], mesh.nodes)
    # FVCOM reorders NVG = (N1, N3, N2) and requires the result to be
    # clockwise (cross < 0) — i.e. the file must be CCW. Check ALL
    # elements, not just #1 as the model does.
    nvg = conn[:, [1, 3, 2]] - 1
    p0, p1, p2 = (mesh.nodes[nvg[:, k]] for k in range(3))
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    assert (cross < 0).all()


def test_write_grd_rejects_flipped(tmp_path):
    mesh = _mesh()
    mesh.elements[0] = mesh.elements[0][[0, 2, 1]]
    with pytest.raises(ValueError, match="CCW"):
        write_grd(mesh, tmp_path / "bad_grd.dat")


def test_write_grd_rejects_orphans(tmp_path):
    base = _mesh()
    mesh = Fort14Mesh(
        title=base.title,
        nodes=np.vstack([base.nodes, [[9e3, 9e3]]]),
        depths=np.append(base.depths, 2.0),
        elements=base.elements,
        open_boundaries=base.open_boundaries,
        land_boundaries=base.land_boundaries,
    )
    with pytest.raises(ValueError, match="compact_nodes"):
        write_grd(mesh, tmp_path / "bad_grd.dat")


# ---------------------------------------------------------------------------
# _dep.dat / _cor.dat / _spg.dat
# ---------------------------------------------------------------------------


def test_write_dep(tmp_path):
    mesh = _mesh()
    p = write_dep(mesh, tmp_path / "case_dep.dat")
    assert _header_int(p, "Node Number") == mesh.n_nodes
    rows = np.array([[float(v) for v in r] for r in _data_rows(p)])
    assert rows.shape == (mesh.n_nodes, 3)
    assert np.allclose(rows[:, :2], mesh.nodes)
    assert np.allclose(rows[:, 2], mesh.depths)


def test_write_cor(tmp_path):
    mesh = _mesh()
    lat = np.full(mesh.n_nodes, 35.3)
    p = write_cor(mesh, tmp_path / "case_cor.dat", lat)
    rows = np.array([[float(v) for v in r] for r in _data_rows(p)])
    assert np.allclose(rows[:, 2], 35.3)
    with pytest.raises(ValueError, match="n_nodes"):
        write_cor(mesh, tmp_path / "bad_cor.dat", lat[:-1])


def test_write_spg(tmp_path):
    mesh = _mesh()
    p = write_spg(mesh, tmp_path / "case_spg.dat")
    assert _header_int(p, "Sponge Node Number") == 0
    p = write_spg(mesh, tmp_path / "case_spg2.dat", [(0, 5000.0, 0.001)])
    rows = _data_rows(p)
    assert rows[0][0] == "1"  # 1-indexed on disk
    with pytest.raises(ValueError, match="out of range"):
        write_spg(mesh, tmp_path / "bad_spg.dat", [(99, 1.0, 1.0)])


# ---------------------------------------------------------------------------
# _obc.dat
# ---------------------------------------------------------------------------


def test_write_obc(tmp_path):
    mesh = _mesh()
    p = write_obc(mesh, tmp_path / "case_obc.dat", obc_type=2)
    assert _header_int(p, "OBC Node Number") == 2
    rows = np.array([[int(v) for v in r] for r in _data_rows(p)])
    assert list(rows[:, 0]) == [1, 2]
    assert list(rows[:, 1]) == [_nid(0, 1) + 1, _nid(0, 2) + 1]
    assert set(rows[:, 2]) == {2}


def test_write_obc_validates_types(tmp_path):
    mesh = _mesh()
    with pytest.raises(ValueError, match="1-10"):
        write_obc(mesh, tmp_path / "bad_obc.dat", obc_type=11)
    with pytest.raises(ValueError, match="entries"):
        write_obc(mesh, tmp_path / "bad_obc.dat", obc_type=[1, 2])


def test_write_obc_empty(tmp_path):
    mesh = _mesh()
    mesh.open_boundaries = []
    p = write_obc(mesh, tmp_path / "case_obc.dat")
    assert _header_int(p, "OBC Node Number") == 0


# ---------------------------------------------------------------------------
# .2dm
# ---------------------------------------------------------------------------


def test_write_2dm(tmp_path):
    mesh = _mesh()
    p = write_2dm(mesh, tmp_path / "case.2dm")
    lines = p.read_text().splitlines()
    assert lines[0] == "MESH2D"
    e3t = [ln.split() for ln in lines if ln.startswith("E3T")]
    nd = [ln.split() for ln in lines if ln.startswith("ND")]
    ns = [ln.split() for ln in lines if ln.startswith("NS")]
    assert len(e3t) == mesh.n_elements
    assert len(nd) == mesh.n_nodes
    conn = np.array([[int(v) for v in r[2:5]] for r in e3t]) - 1
    assert np.array_equal(conn, mesh.elements)
    z = np.array([float(r[4]) for r in nd])
    assert np.allclose(z, mesh.depths)  # default z_convention="depth"
    ns_ids = [int(v) for r in ns for v in r[1:]]
    assert ns_ids == [_nid(0, 1) + 1, -(_nid(0, 2) + 1)]

    p2 = write_2dm(mesh, tmp_path / "case_elev.2dm", z_convention="elevation")
    nd2 = [ln.split() for ln in p2.read_text().splitlines() if ln.startswith("ND")]
    assert np.allclose([float(r[4]) for r in nd2], -mesh.depths)


def test_write_2dm_nodestring_wraps_lines(tmp_path):
    mesh = _mesh()
    long_seg = np.array(
        [_nid(0, 1), _nid(0, 2), _nid(1, 2), _nid(2, 2), _nid(2, 1),
         _nid(2, 0), _nid(1, 0), _nid(0, 0)] + [_nid(1, 1)] * 4,
    )
    mesh.open_boundaries = [long_seg]
    p = write_2dm(mesh, tmp_path / "case.2dm")
    ns = [ln for ln in p.read_text().splitlines() if ln.startswith("NS")]
    assert len(ns) == 2  # 12 ids -> 10 + 2
    assert int(ns[-1].split()[-1]) < 0


# ---------------------------------------------------------------------------
# export_fvcom_case + CLI
# ---------------------------------------------------------------------------


def test_export_fvcom_case(tmp_path):
    mesh = _mesh()
    written = export_fvcom_case(
        mesh, tmp_path, "tokyo",
        obc_type=1,
        cor=np.full(mesh.n_nodes, 35.0),
        write_empty_spg=True,
    )
    assert set(written) == {"grd", "dep", "obc", "cor", "spg", "2dm"}
    for path in written.values():
        assert path.exists()
    assert written["grd"].name == "tokyo_grd.dat"


def test_cli_export(tmp_path, capsys):
    mesh = _mesh()
    f14 = tmp_path / "tokyo.14"
    write_fort14(mesh, f14)
    rc = exportfvcom.main([str(f14), "--cor", "y", "--write-empty-spg"])
    assert rc == 0
    for suffix in ("_grd.dat", "_dep.dat", "_obc.dat", "_cor.dat", "_spg.dat"):
        assert (tmp_path / f"tokyo{suffix}").exists()
    assert (tmp_path / "tokyo.2dm").exists()
    out = capsys.readouterr().out
    assert "grd:" in out


def test_cli_export_refuses_flipped(tmp_path):
    mesh = _mesh()
    mesh.elements[0] = mesh.elements[0][[0, 2, 1]]
    f14 = tmp_path / "bad.14"
    write_fort14(mesh, f14)
    assert exportfvcom.main([str(f14)]) == 1


def test_cli_export_missing_input(tmp_path):
    assert exportfvcom.main([str(tmp_path / "nope.14")]) == 2
