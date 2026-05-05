from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

FIXTURE = Path(__file__).parent / "fixtures" / "tiny.fort14"


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
