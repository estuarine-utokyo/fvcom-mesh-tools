"""Tests for the FVCOM native writers (io/fvcom_native.py, fmesh-export-fvcom).

The reference parsers below mimic the FVCOM 5.1 cold-start readers as
documented in ``docs/fvcom_source_constraints.md``: keyword headers,
free-format rows, and — for ``_grd.dat`` — the ``NVG = (N1, N3, N2)``
column swap followed by the clockwise check (a valid file is CCW, so
the swapped connectivity must be CW for every element).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.cli import exportfvcom
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
from fvcom_mesh_tools.io.fvcom_native import (
    apply_obc_depth_control,
    export_fvcom_case,
    fvcom_next_obc,
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


def test_fvcom_next_obc_picks_most_normal_interior_neighbour():
    # OBC on the x=0 edge: inward normal is +x. Node (0,1) has interior
    # neighbours (0,0) [dot 0], (1,1) [dot 1] and (1,2) [dot 0.707];
    # node (0,2)'s only interior neighbour is (1,2).
    m = _mesh()
    nxt, margin = fvcom_next_obc(m.nodes, m.elements, m.open_boundaries[0])
    assert nxt.tolist() == [_nid(1, 1), _nid(1, 2)]
    assert margin[0] == pytest.approx(1 - np.sqrt(0.5))
    assert np.isinf(margin[1])


def test_fvcom_next_obc_is_orientation_independent():
    m = _mesh()
    cw = m.elements[:, [0, 2, 1]]
    a, _ = fvcom_next_obc(m.nodes, m.elements, m.open_boundaries[0])
    b, _ = fvcom_next_obc(m.nodes, cw, m.open_boundaries[0])
    assert a.tolist() == b.tolist()


def test_fvcom_next_obc_rejects_isolated_obc_node():
    m = _mesh()
    with pytest.raises(ValueError, match="no OBC neighbour"):
        fvcom_next_obc(m.nodes, m.elements, [_nid(0, 1)])


def test_apply_obc_depth_control():
    m = _mesh()
    out, change = apply_obc_depth_control(m)
    obc = m.open_boundaries[0]
    assert out.depths[obc].tolist() == [m.depths[_nid(1, 1)], m.depths[_nid(1, 2)]]
    assert change.tolist() == (out.depths[obc] - m.depths[obc]).tolist()
    interior = np.setdiff1d(np.arange(m.n_nodes), obc)
    assert np.array_equal(out.depths[interior], m.depths[interior])
    assert np.array_equal(m.depths, np.arange(2.0, 2.0 + N * N))  # input untouched


def test_export_applies_obc_depth_control_by_default(tmp_path):
    m = _mesh()
    obc = m.open_boundaries[0]
    on = export_fvcom_case(m, tmp_path / "on", "c", twodm=False)
    off = export_fvcom_case(m, tmp_path / "off", "c", twodm=False, obc_depth_control=False)
    dep_on = np.array([float(r[2]) for r in _data_rows(on["dep"])])
    dep_off = np.array([float(r[2]) for r in _data_rows(off["dep"])])
    assert dep_on[obc].tolist() == [m.depths[_nid(1, 1)], m.depths[_nid(1, 2)]]
    assert np.allclose(dep_off, m.depths)


# --------------------------------------------------------------- readers
#
# Added because local refinement takes a finished FVCOM case as its base: the
# mesh and the depth file are inputs, not something to rebuild.


def _case(tmp_path):
    import numpy as np

    from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case

    nodes = np.array([[0.0, 0.0], [100.0, 0.0], [200.0, 0.0],
                      [0.0, 100.0], [100.0, 100.0], [200.0, 100.0]])
    elements = np.array([[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]])
    mesh = Fort14Mesh(title="case", nodes=nodes,
                      depths=np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
                      elements=elements,
                      open_boundaries=[np.array([0, 1, 2])],
                      land_boundaries=[(20, np.array([2, 5, 4, 3, 0]))])
    written = export_fvcom_case(mesh, tmp_path, "c", cor=nodes[:, 1] * 0 + 35.0,
                                obc_depth_control=False)
    return mesh, written


def test_a_written_case_reads_back_the_same(tmp_path: Path) -> None:
    import numpy as np

    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    mesh, written = _case(tmp_path)
    back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
    assert np.array_equal(back.nodes, mesh.nodes)
    assert np.array_equal(back.depths, mesh.depths)
    assert np.array_equal(np.sort(back.elements, axis=1),
                          np.sort(mesh.elements, axis=1))
    assert np.array_equal(back.open_boundaries[0], mesh.open_boundaries[0])


def test_the_depth_file_wins_over_the_grd_column(tmp_path: Path) -> None:
    """The two disagree on purpose.

    goto2023's `_grd.dat` carries a depth column from whenever it was made --
    4.312072 m at node 1 -- while the b12 baseline names
    `TokyoBay_dep_m7001tp_rfac0p2_cap300.dat`, which says 7.161207 m there.
    Reading the grd's column would silently run the wrong bathymetry.
    """

    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    _, written = _case(tmp_path)
    grd = Path(written["grd"])
    lines = grd.read_text().split("\n")
    out = []
    for ln in lines:
        parts = ln.split()
        out.append(ln if len(parts) != 3 or not ln[0].isdigit()
                   else f"{parts[0]} {parts[1]} {parts[2]} 999.0")
    grd.write_text("\n".join(out))
    back = read_fvcom_case(grd, written["dep"], written["obc"])
    assert back.depths.max() < 100.0


def test_a_depth_file_for_another_mesh_is_refused(tmp_path: Path) -> None:
    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    _, written = _case(tmp_path)
    dep = Path(written["dep"])
    lines = dep.read_text().split("\n")
    lines[1] = "9999.0 9999.0 5.0"
    dep.write_text("\n".join(lines))
    with pytest.raises(ValueError, match="does not sit on"):
        read_fvcom_case(written["grd"], dep, written["obc"])


def test_land_boundaries_are_derived_around_the_open_arc(tmp_path: Path) -> None:
    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    _, written = _case(tmp_path)
    back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
    assert [b for b, _ in back.land_boundaries] == [20]
    # the land string runs from the arc's far end back round to its start
    land = back.land_boundaries[0][1].tolist()
    assert land[0] == 2 and land[-1] == 0
    assert 1 not in land[1:-1], "the interior of the open arc is not land"


def test_an_obc_type_survives_a_round_trip(tmp_path: Path) -> None:
    """It is part of the model input -- odd is elevation only, even adds
    nonlinear flux -- and a base declaring type 2 came back as type 1."""
    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    mesh, _ = _case(tmp_path)
    written = export_fvcom_case(mesh, tmp_path / "two", "t", obc_type=2,
                                cor=mesh.nodes[:, 1] * 0 + 35.0,
                                obc_depth_control=False)
    back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
    assert back.obc_type == 2
    again = export_fvcom_case(back, tmp_path / "three", "u",
                              cor=back.nodes[:, 1] * 0 + 35.0,
                              obc_depth_control=False)
    types = [int(ln.split()[2])
             for ln in Path(again["obc"]).read_text().splitlines()[1:]]
    assert set(types) == {2}


def test_a_non_finite_depth_coordinate_is_refused(tmp_path: Path) -> None:
    """`off.max() > tol` is False for NaN, so one NaN turned the check off."""
    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    _, written = _case(tmp_path)
    dep = Path(written["dep"])
    lines = dep.read_text().split("\n")
    lines[1] = "nan 0.0 5.0"
    dep.write_text("\n".join(lines))
    with pytest.raises(ValueError, match="non-finite"):
        read_fvcom_case(written["grd"], dep, written["obc"])


def test_a_non_finite_grid_coordinate_is_refused(tmp_path: Path) -> None:
    """The dep file was checked for this and the grd file was not.

    A NaN coordinate reaches FVCOM as a mesh it cannot build, and reaches
    every geometric test here as a comparison that is quietly False.
    """
    from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

    _, written = _case(tmp_path)
    grd = Path(written["grd"])
    rows = grd.read_text().splitlines()
    # the first node row: header lines, then every cell, then the nodes
    n_cells = _header_int(grd, "Cell Number")
    fields = rows[2 + n_cells].split()
    fields[1] = "nan"
    rows[2 + n_cells] = " ".join(fields)
    grd.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="non-finite"):
        read_fvcom_case(grd, written["dep"], written["obc"])


def test_the_obc_type_survives_a_dataclass_replace(tmp_path: Path) -> None:
    """``apply_obc_depth_control`` rebuilds the mesh with ``replace``.

    ``obc_type`` rode alongside as a plain attribute, so the rebuild dropped
    it and a type 2 boundary was exported as type 1 -- a different model,
    silently. It is a field now.
    """
    import dataclasses

    from fvcom_mesh_tools.io.fvcom_native import (
        apply_obc_depth_control,
        read_obc_types,
    )

    mesh = dataclasses.replace(_mesh(), obc_type=2)
    controlled, _ = apply_obc_depth_control(mesh)
    assert controlled.obc_type == 2
    written = export_fvcom_case(controlled, tmp_path / "kept", "k",
                                cor=controlled.nodes[:, 1] * 0 + 35.0,
                                obc_depth_control=False)
    assert read_obc_types(written["obc"]) == [2]


def test_the_obc_type_survives_a_rebuild_not_only_a_replace(tmp_path: Path) -> None:
    """`compact_nodes` builds a new mesh instead of replacing one.

    Making `obc_type` a field fixed `dataclasses.replace`, but the cleanup
    helpers construct `Fort14Mesh(...)` directly and did not pass it, so a
    type 3 boundary came back as type 1 from an ordinary compaction (fifth
    review). Every rebuild carries the metadata now.
    """
    import dataclasses

    import numpy as np

    from fvcom_mesh_tools.io.fvcom_native import read_obc_types
    from fvcom_mesh_tools.mesh_clean import compact_nodes

    mesh = dataclasses.replace(_mesh(), obc_type=3)
    # one unused node, which is what compaction is for
    mesh = dataclasses.replace(
        mesh, nodes=np.vstack([mesh.nodes, [[9e4, 9e4]]]),
        depths=np.append(mesh.depths, 5.0))
    out, _ = compact_nodes(mesh)
    assert out.obc_type == 3
    written = export_fvcom_case(out, tmp_path / "compacted", "c",
                                cor=out.nodes[:, 1] * 0 + 35.0,
                                obc_depth_control=False)
    assert read_obc_types(written["obc"]) == [3]
