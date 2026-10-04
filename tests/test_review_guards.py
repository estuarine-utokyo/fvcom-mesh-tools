"""Guards added by review rounds 29-30 (extend tools)."""
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


def _acceptance(tmp_path, merge, generate, nodes=3, products=True):
    """Run 445's acceptance check on report files made in tmp_path (rounds 29 F2, 30 F1)."""
    import json
    from pathlib import Path

    from fvcom_mesh_tools.io.fvcom_native import write_dep, write_grd, write_obc

    src = (Path(__file__).resolve().parents[1] / "notebooks" / "445_extend_mesh.py").read_text()
    start, end = src.index("def _count("), src.index("if not failed and (why := _acceptance")
    gen = tmp_path / "gen"
    gen.mkdir(exist_ok=True)
    recipe = {"case": "C", "recipe_sha256": "r" * 64, "open_boundary_sha256": "o" * 64}
    ns = {"json": json, "OUT": tmp_path, "gen": gen, "recipe": recipe}
    exec(src[start:end], ns)
    if products:
        mesh = Fort14Mesh(title="t", nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                          elements=np.array([[0, 1, 2]]), depths=np.array([5.0, 5.0, 5.0]),
                          open_boundaries=[np.array([0, 1])], land_boundaries=[])
        write_grd(mesh, tmp_path / "C_grd.dat")
        write_dep(mesh, tmp_path / "C_dep.dat")
        write_obc(mesh, tmp_path / "C_obc.dat")
    (tmp_path / "merge.json").write_text(json.dumps(merge))
    (gen / "generate.json").write_text(json.dumps(generate))
    return ns["_acceptance_problem"]()


def test_445_acceptance_check(tmp_path):
    good = {"n_nodes": 3, "n_elements": 1, "qa": {"n_gate_total": 23, "n_gate_failed": 0},
            "problems": []}
    ident = {k: "x" for k in ("base_sha256", "land_sha256", "outer_utm14_sha256")}
    ident |= {"recipe_sha256": "r" * 64, "open_boundary_sha256": "o" * 64}
    gen = {"inputs": ident, "settings": {"a": 1}, "n_nodes": 4, "n_elements": 2}
    assert _acceptance(tmp_path, good, gen) is None
    for bad_merge, bad_gen in (
            (good, None), (good, {}), ({}, gen),
            (good, {**gen, "inputs": {}}),
            (good, {**gen, "inputs": {**ident, "recipe_sha256": "other"}}),
            (good, {**gen, "settings": {}}),
            ({**good, "qa": {"n_gate_total": 23, "n_gate_failed": 0.5}}, gen),
            ({**good, "n_nodes": -1}, gen), ({**good, "n_elements": 0}, gen),
            ({**good, "n_nodes": 9}, gen),          # the products hold 3 nodes
            ({**good, "qa": {"n_gate_total": 23, "n_gate_failed": 1}}, gen),
            ({**good, "problems": ["x"]}, gen)):
        assert _acceptance(tmp_path, bad_merge, bad_gen) is not None
    (tmp_path / "C_grd.dat").write_text("x")           # an unreadable product
    assert "unreadable product" in _acceptance(tmp_path, good, gen, products=False)


def test_round_depths_inside_ignores_the_callers_decimal_context():
    # round 30 F2
    import decimal

    with decimal.localcontext() as ctx:
        ctx.prec = 6
        assert round_depths_inside([3.0], 3.0).tolist() == [3.0]
    assert round_depths_inside([1e25], 1e25).tolist() == [1e25]


def test_m7001_limiter_refuses_what_it_cannot_limit():
    # round 30 F7-F9
    from fvcom_mesh_tools.dem.m7001 import rfactor_smooth

    e = np.array([0]), np.array([1])
    assert rfactor_smooth([1.0, 1.0], *e, hmin=3.0)[0].tolist() == [3.0, 3.0]
    for bad in ([np.nan, 5.0], [1.7e308, 1e307], [-1.0, 5.0]):
        with pytest.raises(ValueError):
            rfactor_smooth(bad, *e)
    n = 62
    h = np.where(np.arange(n) < 15, 3.0, 300.0)
    ei = np.arange(n - 1)
    with pytest.raises(ValueError, match="not reached"):
        rfactor_smooth(h, ei, ei + 1, rmax=0.01, max_iter=50)


def _m7001_file(path, z, dims, var="z"):
    netCDF4 = pytest.importorskip("netCDF4")
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 2)
        ds.createVariable("lat", "f8", ("lat",))[:] = [0.0, 1.0]
        ds.createVariable("lon", "f8", ("lon",))[:] = [0.0, 1.0]
        ds.createVariable(var, "f8", dims)[:] = z


def test_m7001_grid_follows_the_declared_dimension_order(tmp_path):
    # round 30 F5: a square (lon, lat) grid passed a shape test as (lat, lon)
    from fvcom_mesh_tools.dem.m7001 import _grid

    z = np.array([[-10.0, -20.0], [-30.0, -40.0]])
    _m7001_file(tmp_path / "a.nc", z, ("lat", "lon"))
    assert _grid(tmp_path / "a.nc", "z")[2].tolist() == z.tolist()
    _m7001_file(tmp_path / "b.nc", z, ("lon", "lat"))
    assert _grid(tmp_path / "b.nc", "z")[2].tolist() == z.T.tolist()


def test_m7001_interpolation_keeps_a_valid_node_beside_a_missing_one(tmp_path, monkeypatch):
    # round 30 F6: a zero-weight NaN corner discarded the fine grid's own value
    from fvcom_mesh_tools.dem import m7001

    fine = tmp_path / "fine.nc"
    wide = tmp_path / "wide.nc"
    _m7001_file(fine, np.array([[-10.0, -20.0], [-30.0, np.nan]]), ("lat", "lon"),
                 var="elevation")
    _m7001_file(wide, np.full((2, 2), -100.0), ("lat", "lon"))
    monkeypatch.setattr(m7001, "_data_dir", lambda: tmp_path)
    monkeypatch.setattr(m7001, "_FINE", "fine.nc")
    monkeypatch.setattr(m7001, "_WIDE", "wide.nc")
    depth, source = m7001.interpolate_m7001_tp(np.array([0.0]), np.array([0.0]))
    assert depth.tolist() == [10.0] and source.tolist() == [0]
