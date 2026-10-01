from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14

FIXTURE = Path(__file__).parent / "fixtures" / "tiny.fort14"
REFERENCE_MESH = (
    Path(__file__).parent.parent
    / "data"
    / "mesh"
    / "reference"
    / "tokyo_bay"
    / "tb_futtsu20220311.14"
)


@pytest.fixture(scope="module")
def tiny() -> Fort14Mesh:
    return read_fort14(FIXTURE)


def test_returns_fort14mesh(tiny: Fort14Mesh) -> None:
    assert isinstance(tiny, Fort14Mesh)


def test_title_and_counts(tiny: Fort14Mesh) -> None:
    assert tiny.title.strip() == "tiny"
    assert tiny.n_nodes == 4
    assert tiny.n_elements == 2


def test_nodes_and_depths(tiny: Fort14Mesh) -> None:
    assert tiny.nodes.shape == (4, 2)
    assert tiny.depths.shape == (4,)
    np.testing.assert_allclose(tiny.nodes[0], [0.0, 0.0])
    np.testing.assert_allclose(tiny.nodes[3], [2.0, 0.5])
    np.testing.assert_allclose(tiny.depths, [1.0, 1.5, 2.0, 3.0])


def test_elements_are_zero_indexed(tiny: Fort14Mesh) -> None:
    # File element 1 references nodes 1,2,3 -> 0-indexed 0,1,2.
    # File element 2 references nodes 2,4,3 -> 0-indexed 1,3,2.
    assert tiny.elements.shape == (2, 3)
    np.testing.assert_array_equal(tiny.elements[0], [0, 1, 2])
    np.testing.assert_array_equal(tiny.elements[1], [1, 3, 2])


def test_open_boundary_zero_indexed(tiny: Fort14Mesh) -> None:
    assert len(tiny.open_boundaries) == 1
    np.testing.assert_array_equal(tiny.open_boundaries[0], [0, 1])


def test_land_boundaries_preserve_ibtype(tiny: Fort14Mesh) -> None:
    assert len(tiny.land_boundaries) == 2

    ibtype0, ids0 = tiny.land_boundaries[0]
    assert ibtype0 == 0
    np.testing.assert_array_equal(ids0, [2, 0])

    ibtype1, ids1 = tiny.land_boundaries[1]
    assert ibtype1 == 21
    np.testing.assert_array_equal(ids1, [3])


def test_bbox(tiny: Fort14Mesh) -> None:
    xmin, ymin, xmax, ymax = tiny.bbox
    assert (xmin, ymin, xmax, ymax) == (0.0, 0.0, 2.0, 1.0)


def _meshes_equal(a: Fort14Mesh, b: Fort14Mesh) -> None:
    np.testing.assert_allclose(a.nodes, b.nodes, atol=0, rtol=1e-12)
    np.testing.assert_allclose(a.depths, b.depths, atol=0, rtol=1e-9)
    np.testing.assert_array_equal(a.elements, b.elements)
    assert len(a.open_boundaries) == len(b.open_boundaries)
    for x, y in zip(a.open_boundaries, b.open_boundaries):
        np.testing.assert_array_equal(x, y)
    assert len(a.land_boundaries) == len(b.land_boundaries)
    for (ia, idsa), (ib, idsb) in zip(a.land_boundaries, b.land_boundaries):
        assert ia == ib
        np.testing.assert_array_equal(idsa, idsb)


def test_round_trip_tiny(tmp_path: Path, tiny: Fort14Mesh) -> None:
    out = tmp_path / "tiny_out.fort14"
    write_fort14(tiny, out)
    reread = read_fort14(out)
    _meshes_equal(tiny, reread)


@pytest.mark.skipif(not REFERENCE_MESH.exists(), reason="reference mesh symlink not in place")
def test_round_trip_reference_mesh(tmp_path: Path) -> None:
    """Round-trip the real Tokyo Bay reference mesh (95k nodes / 183k elements)."""
    src = read_fort14(REFERENCE_MESH)
    out = tmp_path / "tb_round_trip.fort14"
    write_fort14(src, out)
    reread = read_fort14(out)
    _meshes_equal(src, reread)


def test_node_count_mismatch_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.fort14"
    bad.write_text(
        # NE=1 NP=2, but only 1 node row -> np.loadtxt will pull from the
        # element row instead and the shape check should fail.
        "broken\n"
        "1 2\n"
        "1 0.0 0.0 0.0\n"
        "1 3 1 2 1\n"
        "0 = Number of open boundaries\n"
        "0 = Total number of open boundary nodes\n"
        "0 = Number of normal flow boundaries\n"
        "0 = Total number of land boundary nodes\n"
    )
    with pytest.raises(ValueError):
        read_fort14(bad)


def test_a_computed_depth_round_trips(tmp_path: Path) -> None:
    """Eleven significant figures is not a double.

    Depths read from a fort.14 already fit the old ``.10e`` and round-tripped
    by luck; an interpolated one does not. 5.12345678912345 m came back
    2.3e-11 m different, which is physically nothing and is still a broken
    promise -- and the local-refinement contract is checked on the written
    file (gpt-6-astra review, 2026-09-22).
    """
    import numpy as np

    from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14

    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    depths = np.array([5.12345678912345, 3.14159265358979,
                       1.0 / 3.0, 2.0 ** -20 + 7.0])
    mesh = Fort14Mesh(title="precision", nodes=nodes, depths=depths,
                      elements=np.array([[0, 1, 3], [0, 3, 2]]),
                      open_boundaries=[], land_boundaries=[])
    path = tmp_path / "precision.14"
    write_fort14(mesh, path)
    assert np.array_equal(read_fort14(path).depths, depths)


def _fort14_text(nodes, elems, ids=None, eids=None, etype=3):
    ids = ids or list(range(1, len(nodes) + 1))
    eids = eids or list(range(1, len(elems) + 1))
    lines = ["t", f"{len(elems)} {len(nodes)}"]
    lines += [f"{i} {x} {y} 5.0" for i, (x, y) in zip(ids, nodes)]
    lines += [f"{k} {etype} {a} {b} {c}" for k, (a, b, c) in zip(eids, elems)]
    lines += ["0", "0", "0", "0"]
    return "\n".join(lines) + "\n"


def test_a_single_triangle_mesh_reads(tmp_path):
    """Review of the extend tools, round 9 F11."""
    from fvcom_mesh_tools.io.fort14 import read_fort14

    p = tmp_path / "one.14"
    p.write_text(_fort14_text([(0, 0), (1, 0), (0, 1)], [(1, 2, 3)]))
    m = read_fort14(p)
    assert m.n_nodes == 3 and m.elements.tolist() == [[0, 1, 2]]


@pytest.mark.parametrize("kw, match", [
    ({"ids": [1, 1, 3, 4]}, "node ids"),
    ({"eids": [1, 1]}, "element ids"),
    ({"etype": 4}, "type 3"),
])
def test_record_ids_and_element_type_are_checked(tmp_path, kw, match):
    """Review round 9 F12."""
    from fvcom_mesh_tools.io.fort14 import read_fort14

    p = tmp_path / "bad.14"
    p.write_text(_fort14_text([(0, 0), (1, 0), (0, 1), (1, 1)], [(1, 2, 3), (2, 4, 3)], **kw))
    with pytest.raises(ValueError, match=match):
        read_fort14(p)


def test_boundary_ids_are_checked_on_read_and_write(tmp_path):
    """Review of the extend tools, round 10 F9."""
    import numpy as np

    from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14

    text = _fort14_text([(0, 0), (1, 0), (0, 1)], [(1, 2, 3)]).replace(
        "0\n0\n0\n0\n", "1\n2\n2\n0\n4\n0\n0\n")
    p = tmp_path / "b.14"
    p.write_text(text)
    with pytest.raises(ValueError, match="outside 1..3"):
        read_fort14(p)
    p.write_text(text.replace("1\n2\n2\n0\n4\n", "1\n3\n2\n1\n2\n"))
    with pytest.raises(ValueError, match="NETA"):
        read_fort14(p)
    m = Fort14Mesh("t", np.array([[0.0, 0], [1, 0], [0, 1]]), np.ones(3),
                   np.array([[0, 1, 2]]), [np.array([1.9, 2.9])], [])
    with pytest.raises(ValueError, match="whole number"):
        write_fort14(m, tmp_path / "w.14")
    assert not (tmp_path / "w.14").exists()


def test_write_fort14_leaves_an_existing_file_on_a_bad_mesh(tmp_path):
    """Review of the extend tools, round 11 F5."""
    import os

    import numpy as np

    from fvcom_mesh_tools.io.fort14 import Fort14Mesh, write_fort14

    p = tmp_path / "case.14"
    p.write_text("OLD CASE\n")
    bad = Fort14Mesh("t", np.array([[0.0, 0], [1, 0], [0, 1]]), np.ones(3),
                     np.array([[0, 1, 2, 0]]), [], [])
    with pytest.raises(ValueError, match=r"\(NE, 3\)"):
        write_fort14(bad, p)
    assert p.read_text() == "OLD CASE\n" and list(tmp_path.iterdir()) == [p]
    good = Fort14Mesh("t", bad.nodes, bad.depths, np.array([[0, 1, 2]]), [], [])
    write_fort14(good, p)
    umask = os.umask(0)
    os.umask(umask)
    assert p.stat().st_mode & 0o777 == 0o666 & ~umask
