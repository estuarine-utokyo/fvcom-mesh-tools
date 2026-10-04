"""Guards added in review round 29 (extend tools)."""
import numpy as np
import pytest

from fvcom_mesh_tools.extend import rfactor_smooth_free, round_depths_inside
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, write_fort14


def test_round_depths_inside_keeps_a_bound_that_is_on_the_grid():
    # F8: 0.07 * 100 is 7.000000000000001 in binary, which ceil() pushed past 0.07
    assert round_depths_inside([0.07], 0.07, 0.07, decimals=2).tolist() == [0.07]
    assert round_depths_inside([1.000001], 1.0000008, 1.000001, decimals=6).tolist() == [1.000001]
    with pytest.raises(ValueError, match="no 2-decimal depth"):
        round_depths_inside([0.071], 0.071, 0.0719, decimals=2)


@pytest.mark.parametrize("decimals", [-309, -400, -1, 16, 0.5, True])
def test_round_depths_inside_refuses_unsupported_decimals(decimals):
    # F4: a negative or huge precision gave NaN from finite depths
    with pytest.raises(ValueError, match="decimals"):
        round_depths_inside([4.0], 1.0, decimals=decimals)


def test_rfactor_smooth_free_refuses_floor_and_accumulation_overflow():
    # F3: the floor, and a sum of corrections, overflowed after the depths passed
    with pytest.raises(ValueError, match="overflow"):
        rfactor_smooth_free([1e307, 1.0], [0], [1], [False, True], rmax=0.2, hmin=1.7e308,
                            max_iter=10)
    n = 12
    h0 = np.r_[1.0, np.full(n, 4e307)]
    with pytest.raises(ValueError, match="overflow"):
        rfactor_smooth_free(h0, np.zeros(n, int), np.arange(1, n + 1),
                            np.r_[True, np.zeros(n, bool)], rmax=0.2, hmin=1.0, max_iter=10)


@pytest.mark.parametrize("which", ["nodes", "depths", "elements"])
def test_write_fort14_refuses_masked_values(tmp_path, which):
    # F7: a masked value was written as nan or "--"
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    depths = np.array([5.0, 5.0, 5.0])
    elements = np.array([[0, 1, 2]])
    arrays = {"nodes": nodes, "depths": depths, "elements": elements}
    a = arrays[which]
    mask = np.zeros(a.shape, bool)
    mask.flat[0] = True
    arrays[which] = np.ma.array(a, mask=mask)
    mesh = Fort14Mesh(title="t", nodes=arrays["nodes"], elements=arrays["elements"],
                      depths=arrays["depths"], open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError, match="masked"):
        write_fort14(mesh, tmp_path / "x.14")
    assert not (tmp_path / "x.14").exists()


def test_write_fort14_refuses_nonfinite_values(tmp_path):
    mesh = Fort14Mesh(title="t", nodes=np.array([[0.0, np.nan], [1.0, 0.0], [0.0, 1.0]]),
                      elements=np.array([[0, 1, 2]]), depths=np.array([5.0, 5.0, 5.0]),
                      open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError, match="finite"):
        write_fort14(mesh, tmp_path / "x.14")


def test_m7001_grid_turns_masked_cells_into_nan(tmp_path):
    # F1: a -9999 fill under the mask became a 9999 m depth
    netCDF4 = pytest.importorskip("netCDF4")
    from fvcom_mesh_tools.dem.m7001 import _grid

    f = tmp_path / "g.nc"
    with netCDF4.Dataset(f, "w") as ds:
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 2)
        ds.createVariable("lat", "f8", ("lat",))[:] = [0.0, 1.0]
        ds.createVariable("lon", "f8", ("lon",))[:] = [0.0, 1.0]
        z = ds.createVariable("z", "f8", ("lat", "lon"), fill_value=-9999.0)
        z[:] = np.array([[-5.0, -9999.0], [-5.0, -5.0]])
    _, _, z = _grid(f, "z")
    assert np.isnan(z[0, 1]) and z[0, 0] == -5.0


def _acceptance(tmp_path, merge, generate):
    """Run 445's acceptance check on report files made in tmp_path (round 29 F2)."""
    import json
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "notebooks" / "445_extend_mesh.py").read_text()
    start, end = src.index("def _count("), src.index("if not failed and (why := _acceptance")
    gen = tmp_path / "gen"
    gen.mkdir(exist_ok=True)
    ns = {"json": json, "OUT": tmp_path, "gen": gen, "recipe": {"case": "C"}}
    exec(src[start:end], ns)
    for kind in ("grd", "dep", "obc"):
        (tmp_path / f"C_{kind}.dat").write_text("x")
    (tmp_path / "merge.json").write_text(json.dumps(merge))
    (gen / "generate.json").write_text(json.dumps(generate))
    return ns["_acceptance_problem"]()


def test_445_acceptance_check(tmp_path):
    good = {"n_nodes": 5, "n_elements": 3, "qa": {"n_gate_total": 23, "n_gate_failed": 0},
            "problems": []}
    gen = {"inputs": {}, "settings": {}, "n_nodes": 4, "n_elements": 2}
    assert _acceptance(tmp_path, good, gen) is None
    for bad_merge, bad_gen in (
            (good, None), (good, {}), ({}, gen),
            ({**good, "qa": {"n_gate_total": 23, "n_gate_failed": 0.5}}, gen),
            ({**good, "n_nodes": -1}, gen), ({**good, "n_elements": 0}, gen),
            ({**good, "qa": {"n_gate_total": 23, "n_gate_failed": 1}}, gen),
            ({**good, "problems": ["x"]}, gen)):
        assert _acceptance(tmp_path, bad_merge, bad_gen) is not None
