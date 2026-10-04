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


def _acceptance(tmp_path, merge, generate, nodes=3, products=True):
    """Run 445's acceptance check on report files made in tmp_path (rounds 29 F2, 30 F1)."""
    import json
    from pathlib import Path

    from fvcom_mesh_tools.io.fvcom_native import write_dep, write_grd, write_obc

    nb = Path(__file__).resolve().parents[1] / "notebooks" / "445_extend_mesh.py"
    src = nb.read_text()
    start, end = src.index("def _count("), src.index("if not failed and (why := _acceptance")
    gen = tmp_path / "gen"
    gen.mkdir(exist_ok=True)
    recipe = {"case": "C", "recipe_sha256": "r" * 64, "open_boundary_sha256": "o" * 64}
    ns = {"json": json, "np": np, "OUT": tmp_path, "gen": gen, "recipe": recipe}
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
    good = {"n_nodes": 3, "n_elements": 1, "n_open_boundary_nodes": 2,
            "qa": {"n_gate_total": 23, "n_gate_failed": 0}, "problems": []}
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
            ({**good, "n_open_boundary_nodes": 5}, gen),
            ({**good, "qa": {"n_gate_total": 23, "n_gate_failed": 1}}, gen),
            ({**good, "problems": ["x"]}, gen)):
        assert _acceptance(tmp_path, bad_merge, bad_gen) is not None
    (tmp_path / "C_grd.dat").write_text("x")           # an unreadable product
    assert "unreadable product" in _acceptance(tmp_path, good, gen, products=False)
    # a depth file of another mesh beside a good grid (round 31 F1)
    _acceptance(tmp_path, good, gen)
    (tmp_path / "C_dep.dat").write_text("Node Number = 1\n77 88 5\n")
    assert "unreadable product" in _acceptance(tmp_path, good, gen, products=False)


def test_round_depths_inside_ignores_the_callers_decimal_context():
    # round 30 F2
    import decimal

    with decimal.localcontext() as ctx:
        ctx.prec = 6
        assert round_depths_inside([3.0], 3.0).tolist() == [3.0]
    assert round_depths_inside([1e25], 1e25).tolist() == [1e25]
    # round 31 F2: nor its traps or exponent limits
    with decimal.localcontext() as ctx:
        ctx.traps[decimal.Inexact] = True
        ctx.Emax = 3
        assert round_depths_inside([3.0000004], 3.0000004).tolist() == [3.000001]


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


def _m7001_file(path, z, dims, var="z", axes=("lat", "lon"), fill=None, axvals=None):
    netCDF4 = pytest.importorskip("netCDF4")
    with netCDF4.Dataset(path, "w") as ds:
        for name, n in zip(dims, np.shape(z)):
            ds.createDimension(name, n)
        vals = axvals or {"lat": [0.0, 1.0], "lon": [0.0, 1.0]}
        for ax, dim in zip(("lat", "lon"), axes):
            if dim not in ds.dimensions:
                ds.createDimension(dim, 2)
            ds.createVariable(ax, "f8", (dim,))[:] = vals[ax]
        v = ds.createVariable(var, "f8", dims, fill_value=fill)
        v[:] = z


@pytest.fixture
def m7001(tmp_path, monkeypatch):
    """interpolate_m7001_tp on two tiny grids: fine.nc (elevation) and wide.nc (z)."""
    from fvcom_mesh_tools.dem import m7001 as mod

    def make(fine, wide=None, **kw):
        _m7001_file(tmp_path / "fine.nc", fine, **{"dims": ("lat", "lon"), "var": "elevation",
                                                    **kw})
        _m7001_file(tmp_path / "wide.nc", np.full((2, 2), -100.0) if wide is None else wide,
                    dims=("lat", "lon"))
        monkeypatch.setattr(mod, "_data_dir", lambda: tmp_path)
        monkeypatch.setattr(mod, "_FINE", "fine.nc")
        monkeypatch.setattr(mod, "_WIDE", "wide.nc")
        return mod.interpolate_m7001_tp

    return make


def test_m7001_follows_the_declared_dimension_order(m7001):
    # round 30 F5: a square (lon, lat) grid passed a shape test as (lat, lon)
    z = np.array([[-10.0, -20.0], [-30.0, -40.0]])
    f = m7001(z.T, dims=("lon", "lat"))
    assert f(np.array([1.0]), np.array([0.0]))[0].tolist() == [20.0]      # z[lat 0, lon 1]
    f = m7001(z, dims=("lat", "lon"))
    assert f(np.array([1.0]), np.array([0.0]))[0].tolist() == [20.0]


def test_m7001_accepts_coordinate_dimension_aliases(m7001):
    # round 31 F4: lat(y), lon(x), elevation(y, x)
    f = m7001(np.array([[-10.0, -20.0], [-30.0, -40.0]]), dims=("y", "x"), axes=("y", "x"))
    assert f(np.array([0.0]), np.array([0.0]))[0].tolist() == [10.0]


def test_m7001_masked_cells_are_missing_not_depths(m7001):
    # round 29 F1: a -9999 fill under the mask became a 9999 m depth
    f = m7001(np.array([[-5.0, -9999.0], [-5.0, -5.0]]), fill=-9999.0)
    depth, source = f(np.array([1.0]), np.array([0.0]))
    assert depth.tolist() == [100.0] and source.tolist() == [1]        # the wide grid's value


def test_m7001_keeps_a_valid_node_beside_a_missing_one(m7001):
    # round 30 F6: a zero-weight NaN corner discarded the fine grid's own value
    f = m7001(np.array([[-10.0, -20.0], [-30.0, np.nan]]))
    depth, source = f(np.array([0.0]), np.array([0.0]))
    assert depth.tolist() == [10.0] and source.tolist() == [0]


def test_m7001_refuses_bad_axes_and_queries(m7001):
    # round 31 F6, F7
    z = np.array([[-10.0, -20.0], [-30.0, -40.0]])
    f = m7001(z, axvals={"lat": [0.0, 1.0], "lon": [0.0, np.inf]})
    with pytest.raises(ValueError):
        f(np.array([1.0]), np.array([0.0]))
    f = m7001(z)
    with pytest.raises(ValueError, match="masked"):
        f(np.ma.array([0.0], mask=[True]), np.array([0.0]))
    with pytest.raises(ValueError, match="matching"):
        f(np.zeros((2, 2)), np.zeros((2, 2)))


def test_m7001_node_edges_and_production_bounds_are_checked():
    # round 31 F3, F8, F9
    from fvcom_mesh_tools.dem.m7001 import node_edges, production_depths, rfactor_smooth

    with pytest.raises(ValueError):
        node_edges(np.array([[0.0, 1.0, 2.9]]))
    with pytest.raises(ValueError, match="masked"):
        node_edges(np.ma.array([[0, 1, 2]], mask=[[False, False, True]]))
    with pytest.raises(ValueError, match="masked"):
        rfactor_smooth([10.0, 100.0], np.ma.array([0], mask=[True]), np.array([1]), rmax=0.2)
    for hmax in (np.nan, -1.0, 2.0):
        with pytest.raises(ValueError, match="bounds"):
            production_depths(np.zeros(3), np.zeros(3), np.array([[0, 1, 2]]), hmin=3.0,
                              hmax=hmax)


def test_band_field_refuses_an_invalid_half_width():
    # round 31 F10
    from fvcom_mesh_tools.extend import band_field

    x, y = np.meshgrid([0.0, 1.0], [0.0, 1.0])
    line = np.array([[0.0, 0.0], [1.0, 0.0]])
    for w in (-1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="half_width_m"):
            band_field(x, y, line, np.array([10.0, 10.0]), half_width_m=w)
