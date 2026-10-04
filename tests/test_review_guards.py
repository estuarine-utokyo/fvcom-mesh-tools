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


def _acceptance(tmp_path, merge, generate, products=True, boundary=None):
    """Run 445's acceptance check on report files made in tmp_path (rounds 29 F2, 30 F1)."""
    import json
    from pathlib import Path

    from fvcom_mesh_tools.io.fvcom_native import write_dep, write_grd, write_obc

    nb = Path(__file__).resolve().parents[1] / "notebooks" / "445_extend_mesh.py"
    src = nb.read_text()
    start, end = src.index("def _count("), src.index("if not failed and (why := _acceptance")
    gen = tmp_path / "gen"
    gen.mkdir(exist_ok=True)
    from pyproj import Transformer

    back = Transformer.from_crs(32654, 4326, always_xy=True)
    lonlat = np.column_stack(back.transform(np.array([500000.0, 501000.0]),
                                            np.array([3900000.0, 3900000.0])))
    recipe = {"case": "C", "recipe_sha256": "r" * 64, "open_boundary_sha256": "o" * 64,
              "open_boundary_lonlat": lonlat if boundary is None else boundary(lonlat)}
    ns = {"json": json, "np": np, "OUT": tmp_path, "gen": gen, "recipe": recipe}
    exec(src[start:end], ns)
    if products:
        mesh = Fort14Mesh(title="t", nodes=np.array([[500000.0, 3900000.0], [501000.0, 3900000.0],
                                          [500000.0, 3901000.0]]),
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
    # the recipe's boundary, either way round, and nothing else (round 32 F1)
    assert _acceptance(tmp_path, good, gen, boundary=lambda ll: ll[::-1]) is None
    assert "recipe's boundary" in _acceptance(tmp_path, good, gen, boundary=lambda ll: ll + 1e-3)
    assert "recipe's boundary" in _acceptance(tmp_path, good, gen, boundary=lambda ll: ll[:1])
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


def test_node_edges_wants_three_columns_even_when_empty():
    from fvcom_mesh_tools.dem.m7001 import node_edges

    assert node_edges(np.empty((0, 3), int))[0].size == 0       # (0, 3) is legitimate
    with pytest.raises(ValueError, match="M, 3"):
        node_edges(np.empty((0, 4), int))


def test_production_depths_names_a_missing_node_before_sampling():
    from fvcom_mesh_tools.dem.m7001 import production_depths

    with pytest.raises(ValueError, match="name a node"):
        production_depths(np.zeros(3), np.zeros(3), np.array([[0, 1, 3]]))


def test_band_field_wants_a_real_width_and_a_line_with_length():
    from fvcom_mesh_tools.extend import band_field

    x, y = np.meshgrid([0.0, 1.0], [0.0, 1.0])
    line = np.array([[0.0, 0.0], [1.0, 0.0]])
    for w in (np.complex128(1 + 7j), np.complex128(7j), "2", None):
        with pytest.raises(ValueError, match="half_width_m"):
            band_field(x, y, line, np.array([10.0, 10.0]), half_width_m=w)
    with pytest.raises(ValueError, match="zero-length"):
        band_field(x, y, np.zeros((2, 2)), np.array([10.0, 20.0]), 2.0)


def test_obc_band_refuses_unknown_values_and_bad_controls():
    # round 32 F5, F6, F9
    from fvcom_mesh_tools.obc_band import build_obc_band

    arc = np.c_[139.0 + np.arange(6) * 0.01, np.full(6, 35.0)]
    h = np.full(6, 500.0)
    assert build_obc_band(arc, h, taper="local")["inner_ll"].shape == (4, 2)
    with pytest.raises(ValueError, match="masked"):
        build_obc_band(arc, np.ma.array(h, mask=[0, 0, 0, 0, 0, 1]), taper="local")
    with pytest.raises(ValueError, match="finite"):
        build_obc_band(arc, np.array([10, 10, np.nan, 10, 10, 10.0]), taper="local")
    for kw in ({"k_offset": -1.0}, {"k_offset": np.nan}, {"skip_ends": -1},
               {"smooth_passes": -1}, {"skip_ends": 1.5}):
        with pytest.raises(ValueError):
            build_obc_band(arc, h, **kw)
    with pytest.raises(ValueError, match="fewer than two"):
        build_obc_band(arc[:3], h[:3])
    with pytest.raises(ValueError, match="repeated"):
        build_obc_band(np.repeat(arc[:3], 2, axis=0), np.full(6, 500.0))


def test_corridor_targets_reach_the_last_endpoint_with_its_own_target():
    # round 32 F7
    from fvcom_mesh_tools.obc_band import corridor_targets

    arc = np.c_[139 + np.arange(3) * 0.0005, np.full(3, 35.0)]
    pts, tgt = corridor_targets(arc, [100.0, 100.0, 1000.0], step_m=100.0)
    assert tgt[-1] == pytest.approx(1000.0 / 1.2)
    closure = np.c_[139.001 + np.arange(3) * 0.0003, np.full(3, 35.001)]
    pts, tgt = corridor_targets(arc, [100.0, 100.0, 1000.0], closure_ll=closure,
                                h_closure_end_m=300.0, step_m=100.0)
    assert tgt[-1] == pytest.approx(300.0 / 1.2)


def test_round33_guards():
    from fvcom_mesh_tools.dem.m7001 import production_depths, rfactor_smooth
    from fvcom_mesh_tools.extend import rfactor_smooth_free
    from fvcom_mesh_tools.obc_band import apply_corridor, build_obc_band, corridor_targets

    # F2: a repeated pair of arc nodes
    arc = np.c_[139 + 0.01 * np.array([0, 1, 2, 2, 3, 4]), np.full(6, 35.0)]
    with pytest.raises(ValueError, match="repeated"):
        build_obc_band(arc, np.full(6, 500.0), taper="local")
    # F3: a closure shorter than a metre still ends on its exact target
    arc3 = np.c_[139 + np.arange(3) * 0.0005, np.full(3, 35.0)]
    short = np.array([[139.001, 35.0], [139.001, 35.0 + 0.5 / 111e3]])
    _, tgt = corridor_targets(arc3, [100.0, 100.0, 1000.0], closure_ll=short,
                              h_closure_end_m=300.0)
    assert tgt[-1] == pytest.approx(250.0)
    with pytest.raises(ValueError, match="no length"):
        corridor_targets(arc3, [100.0, 100.0, 1000.0], closure_ll=short[:1].repeat(2, 0),
                         h_closure_end_m=300.0)
    # F4: step, factor, masks
    for kw in ({"step_m": -100.0}, {"step_m": 0.0}, {"mesh_factor": np.nan}):
        with pytest.raises(ValueError):
            corridor_targets(arc3, [100.0, 100.0, 1000.0], **kw)
    with pytest.raises(ValueError, match="masked"):
        corridor_targets(arc3, np.ma.array([100.0, 100.0, 1000.0], mask=[0, 0, 1]))
    # F5: scalar controls of both limiters
    e = np.array([0]), np.array([1])
    for kw in ({"max_iter": 1.9}, {"hmin": True}, {"rmax": np.complex128(0.2 + 0.1j)}):
        with pytest.raises(ValueError):
            rfactor_smooth_free(np.array([100.0, 10.0]), *e, np.array([True, True]),
                                **{"rmax": 0.2, "hmin": 1.0, **kw})
        with pytest.raises(ValueError):
            rfactor_smooth(np.array([100.0, 10.0]), *e, **kw)
    # F6: nothing to limit
    with pytest.raises(ValueError, match="non-empty"):
        production_depths(np.zeros(3), np.zeros(3), np.empty((0, 3), int))
    # F7: gradation
    for g in (-0.2, np.nan):
        with pytest.raises(ValueError, match="gradation"):
            apply_corridor(np.zeros((1, 2)), np.zeros((1, 2)), np.zeros((1, 2)),
                           np.array([[0.0, 0.0]]), np.array([100.0]), grade=g,
                           arc_mean_lat=35.0)


def test_round34_guards():
    from fvcom_mesh_tools.obc_band import apply_corridor, build_obc_band, corridor_targets
    from fvcom_mesh_tools.obc_design import fillet, resample

    # F1: distinct arc nodes whose offsets fold the guide onto one point
    arc = 0.001 * np.array([[-1, -1], [0, -1], [1, 0], [0, 1], [-1, 0], [-2, 1]], float)
    with pytest.raises(ValueError, match="coincident"):
        build_obc_band(arc, np.full(6, 88.8), skip_ends=2, taper="local")
    # F2: masks, a NaN in the field, shapes
    z = np.zeros((1, 2))
    kw = {"grade": 0.2, "arc_mean_lat": 0.0}
    pts, tgt = np.array([[0.0, 0.0]]), np.array([100.0])
    with pytest.raises(ValueError, match="masked"):
        apply_corridor(z, z, np.ma.array([[0.1, 0.0]], mask=[[1, 0]]), pts, tgt, **kw)
    with pytest.raises(ValueError, match="finite"):
        apply_corridor(z, z, np.array([[np.nan, 0.0]]), pts, tgt, **kw)
    with pytest.raises(ValueError, match="shape"):
        apply_corridor(z, z, np.zeros((2, 2)), pts, tgt, **kw)
    # F4: a scaling that overflows
    with pytest.raises(ValueError, match="overflowed"):
        corridor_targets(np.array([[0.0, 0.0], [0.001, 0.0]]), [100.0, 100.0],
                         mesh_factor=1e-320)
    # F5: planar helpers want two columns
    line = np.array([[0.0, 0, 0], [100.0, 0, 50]])
    with pytest.raises(ValueError, match="N >= 2, 2"):
        fillet(line, [])
    with pytest.raises(ValueError, match="N >= 2, 2"):
        resample(line, 10)


def _lattice_mesh(n=6, shear=0.37):
    """A triangulated n x n lattice; row 0 is the extension, rows 1.. the base."""
    ij = [(i, j) for j in range(n) for i in range(n)]
    nodes = np.array([[1000.0 * (i + shear * j), 1000.0 * j] for i, j in ij])
    tris = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b, c, d = j * n + i, j * n + i + 1, (j + 1) * n + i, (j + 1) * n + i + 1
            tris += [[a, b, d], [a, d, c]]
    return nodes, np.array(tris), n


def test_perpendicularity_repair_never_moves_a_frozen_node():
    # round 35 F1
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local

    nodes, tris, n = _lattice_mesh()
    mesh = Fort14Mesh(title="t", nodes=nodes, elements=tris, depths=np.full(len(nodes), 10.0),
                      open_boundaries=[np.arange(n)], land_boundaries=[])
    frozen = np.arange(len(nodes)) >= n                      # rows 1.. are the base
    out, _ = align_open_boundary_local(mesh, movable=~frozen)
    assert np.array_equal(out.nodes[frozen], nodes[frozen])
    moved, _ = align_open_boundary_local(mesh)               # without the mask they move
    assert not np.array_equal(moved.nodes[frozen], nodes[frozen])
    with pytest.raises(ValueError, match="movable"):
        align_open_boundary_local(mesh, movable=np.ones(len(nodes), int))


def test_round35_guards():
    from fvcom_mesh_tools.extend import check_land_cover
    from fvcom_mesh_tools.obc_band import apply_corridor

    with pytest.raises(ValueError, match="hmin|real"):
        round_depths_inside([4.0], 3.0, np.complex128(5 + 8j))
    with pytest.raises(ValueError, match="real"):
        round_depths_inside([4.0], True)
    z = np.zeros((1, 2))
    for t in (-100.0, 0.0):
        with pytest.raises(ValueError, match="positive"):
            apply_corridor(z, z, np.full((1, 2), 1e-4), np.array([[0.0, 0.0]]),
                           np.array([t]), grade=0.2, arc_mean_lat=0.0)
    import shapely

    mesh = Fort14Mesh(title="t", nodes=np.array([[0.0, 0], [100, 0], [0, 100]]),
                      elements=np.array([[0, 1, 2]]), depths=np.ones(3),
                      open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError):
        check_land_cover(mesh, shapely.Polygon([(0, 0), (100, 0), (0, 100)]), True, erode=0.1)


def test_round36_guards():
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.dem.m7001 import rfactor_smooth
    from fvcom_mesh_tools.extend import check_no_overlap, compose_sizing
    from fvcom_mesh_tools.obc_design import fillet

    # F1: a masked permission is refused, nothing moves
    nodes, tris, n = _lattice_mesh()
    mesh = Fort14Mesh(title="t", nodes=nodes, elements=tris, depths=np.full(len(nodes), 10.0),
                      open_boundaries=[np.arange(n)], land_boundaries=[])
    movable = np.ma.array(np.ones(len(nodes), bool), mask=np.arange(len(nodes)) >= n)
    with pytest.raises(ValueError, match="masked"):
        align_open_boundary_local(mesh, movable=movable)
    # F2: complex arrays
    with pytest.raises(ValueError, match="complex"):
        round_depths_inside(np.array([4 + 9j]), 3.0)
    with pytest.raises(ValueError, match="complex"):
        rfactor_smooth_free(np.array([100 + 9j, 10 + 2j]), [0], [1], [True, True], rmax=0.2,
                            hmin=1.0)
    with pytest.raises(ValueError, match="complex"):
        rfactor_smooth(np.array([100 + 9j, 10 + 2j]), np.array([0]), np.array([1]))
    x, y = np.meshgrid([0.0, 1000.0], [0.0])
    with pytest.raises(ValueError, match="complex"):
        compose_sizing(np.array([[100 + 9j, 100 + 2j]]), x, y, grade=0.2)
    # F3: finite inputs that overflow
    with pytest.raises(ValueError, match="finite"):
        fillet(np.array([[0.0, 0.0], [1e200, 0.0], [1e200, 1e200]]), [1e199])
    # F4: a bool, or a tolerance that exempts complete overlap
    import shapely

    del shapely
    mesh2 = Fort14Mesh(title="t", nodes=np.array([[0.0, 0], [100, 0], [0, 100], [10, 10],
                                                  [30, 10], [10, 30]]),
                       elements=np.array([[0, 1, 2], [3, 4, 5]]), depths=np.ones(6),
                       open_boundaries=[], land_boundaries=[])
    for tol in (True, 1, 1.5):
        with pytest.raises(ValueError, match="rel_tol"):
            check_no_overlap(mesh2, 1, rel_tol=tol)
    with pytest.raises(ValueError, match="overlap the base"):
        check_no_overlap(mesh2, 1)


def test_round37_guards(tmp_path):
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.dem.sources import sample
    from fvcom_mesh_tools.extend import check_land_cover
    from fvcom_mesh_tools.io.fvcom_native import write_cor, write_dep, write_grd
    from fvcom_mesh_tools.obc_band import build_obc_band
    from fvcom_mesh_tools.obc_design import fillet, ray_intersection, resample

    # F1: complex inputs of the geometry and sizing interfaces
    c = np.array([[0, 0], [10, 0]], complex) + 9j
    with pytest.raises(ValueError, match="complex"):
        resample(c, 3)
    with pytest.raises(ValueError, match="complex"):
        ray_intersection(np.array([0 + 9j, 0]), [1, 0], [1, -1], [0, 1])
    with pytest.raises(ValueError, match="complex"):
        fillet([[0, 0], [10, 0], [10, 10]], [np.complex128(2 + 9j)])
    arc = np.c_[139.0 + np.arange(6) * 0.01, np.full(6, 35.0)]
    with pytest.raises(ValueError, match="complex"):
        build_obc_band(arc, np.full(6, 100 + 9j))
    with pytest.raises(ValueError, match="complex"):
        sample(["srtm15plus"], np.array([0.5 + 9j]), np.array([0.5]))
    # F2: native writers
    mesh = Fort14Mesh(title="t", nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                      elements=np.array([[0, 1, 2]]), depths=np.array([5.0, 6.0, 7.0]),
                      open_boundaries=[], land_boundaries=[])
    bad = Fort14Mesh(title="t", nodes=mesh.nodes, elements=mesh.elements,
                     depths=mesh.depths + 9j, open_boundaries=[], land_boundaries=[])
    badxy = Fort14Mesh(title="t", nodes=mesh.nodes + 9j, elements=mesh.elements,
                       depths=mesh.depths, open_boundaries=[], land_boundaries=[])
    for fn, m in ((write_dep, bad), (write_grd, badxy)):
        with pytest.raises(ValueError, match="complex"):
            fn(m, tmp_path / "x.dat")
    with pytest.raises(ValueError, match="complex"):
        write_cor(mesh, tmp_path / "x.dat", np.full(3, 35 + 9j))
    with pytest.raises(ValueError, match="complex"):
        write_fort14(bad, tmp_path / "x.14")
    assert not list(tmp_path.glob("x.*"))
    # F3: a bool erosion
    import shapely

    m1 = Fort14Mesh(title="t", nodes=np.array([[0.0, 0], [100, 0], [0, 100]]),
                    elements=np.array([[0, 1, 2]]), depths=np.ones(3),
                    open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError, match="erode"):
        check_land_cover(m1, shapely.Polygon([(0, 0), (100, 0), (0, 100)]), 0, erode=True)
    # F4: controls, and what remains read from the mesh returned
    nodes, tris, n = _lattice_mesh()
    m6 = Fort14Mesh(title="t", nodes=nodes, elements=tris, depths=np.full(len(nodes), 10.0),
                    open_boundaries=[np.arange(n)], land_boundaries=[])
    for kw in ({"min_angle": np.nan}, {"max_area_change": np.nan}, {"dev_max": np.nan},
               {"max_outer": 0}, {"max_angle": 20.0}, {"w_steps": ()}, {"seed": -1}):
        with pytest.raises(ValueError):
            align_open_boundary_local(m6, **kw)
    _, info = align_open_boundary_local(m6, movable=np.zeros(len(nodes), bool))
    assert info["remaining"]                       # nothing could move: still violating


def test_round38_guards(tmp_path):
    from fvcom_mesh_tools._checks import no_complex
    from fvcom_mesh_tools.io.fvcom_native import _check_exportable, write_spg
    from fvcom_mesh_tools.obc_band import apply_corridor
    from fvcom_mesh_tools.obc_design import fillet, resample

    # F1: a radius iterator is read once
    vertices = [[0, 0], [10, 0], [10, 10]]
    assert fillet(vertices, iter([2.0])).shape == fillet(vertices, [2.0]).shape
    # F2: complex numbers inside an object array
    obj = np.array([np.complex128(4 + 9j)], dtype=object)
    with pytest.raises(ValueError, match="complex"):
        no_complex(a=obj)
    with pytest.raises(ValueError, match="complex"):
        round_depths_inside(obj, 3)
    no_complex(a=np.array([1.0, 2.0], dtype=object))          # real objects pass
    # F3: spacing and corridor inputs
    with pytest.raises(ValueError, match="complex"):
        resample([[0, 0], [10, 0]], np.complex128(3 + 9j))
    with pytest.raises(ValueError, match="complex"):
        resample([[0, 0], [10, 0]], lambda xy: np.full(len(xy), 3 + 9j))
    z = np.zeros((1, 2))
    with pytest.raises(ValueError, match="complex"):
        apply_corridor(z + 9j, z, np.full((1, 2), 1e-4), np.array([[0.0, 0.0]]),
                       np.array([100.0]), grade=0.2, arc_mean_lat=0.0)
    # F4: sponge records
    mesh = Fort14Mesh(title="t", nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                      elements=np.array([[0, 1, 2]]), depths=np.array([5.0, 6.0, 7.0]),
                      open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError, match="complex"):
        write_spg(mesh, tmp_path / "x.spg", np.array([[0, 10, 0.001]], complex) + 9j)
    assert not list(tmp_path.glob("x.*"))
    # F5: the validator keeps its docstring
    assert _check_exportable.__doc__


def test_round39_guards(tmp_path):
    import shapely

    from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast
    from fvcom_mesh_tools.io.fvcom_native import (
        apply_obc_depth_control,
        fvcom_next_obc,
        write_grd,
    )
    from fvcom_mesh_tools.patch import improve_patch

    # F1: masked or mistyped permissions
    nodes, tris, n = _lattice_mesh()
    nodes = nodes.copy()
    every = np.ones(len(nodes), bool)
    faces = np.ones(len(tris), bool)
    masked = np.ma.array(every, mask=every)
    with pytest.raises(ValueError, match="masked"):
        improve_patch(nodes, tris, masked, faces, rounds=1)
    with pytest.raises(ValueError, match="masked"):
        improve_patch(nodes, tris, every, np.ma.array(faces, mask=faces), rounds=1)
    with pytest.raises(ValueError, match="boolean"):
        improve_patch(nodes, tris, every.astype(int), faces, rounds=1)
    # F2: masked and complex coordinates, fractional connectivity
    tri1 = np.array([[0, 1, 2]])
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    land = shapely.Polygon([(-1, -1), (2, -1), (2, 2), (-1, 2)])
    for bad in (np.ma.array(xy, mask=True), xy + 9j):
        with pytest.raises(ValueError, match="masked|complex"):
            improve_patch(bad, tri1, np.zeros(3, bool), np.zeros(1, bool), rounds=1)
        with pytest.raises(ValueError, match="masked|complex"):
            fit_boundary_to_coast(bad, tri1, land, fixed=np.array([0, 1, 2]))
    with pytest.raises(ValueError, match="whole"):
        improve_patch(xy, np.array([[0.0, 1.0, 2.9]]), np.zeros(3, bool), np.zeros(1, bool),
                      rounds=1)
    # F3: depth control and next-OBC
    with pytest.raises(ValueError, match="masked|complex"):
        fvcom_next_obc(np.ma.array(xy, mask=True), tri1, [0, 1])
    with pytest.raises(ValueError, match="whole"):
        apply_obc_depth_control(Fort14Mesh(
            title="t", nodes=xy, elements=tri1, depths=np.ones(3),
            open_boundaries=[np.array([0.5, 1.0])], land_boundaries=[]))
    # F4: an overflowed area is not an orientation verdict
    huge = Fort14Mesh(title="t", nodes=np.array([[0.0, 0.0], [2e200, 3e200], [1e200, 1e200]]),
                      elements=tri1, depths=np.ones(3), open_boundaries=[], land_boundaries=[])
    with pytest.raises(ValueError, match="not finite"):
        write_grd(huge, tmp_path / "x.dat")
    assert not list(tmp_path.glob("x.*"))


def test_round40_guards():
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast
    from fvcom_mesh_tools.io.fvcom_native import apply_obc_depth_control, fvcom_next_obc
    from tests.test_coast_fit import _land_above, _strip_mesh

    nodes, tri = _strip_mesh(nx=3, ny=2)
    land = _land_above(240, nx=3)
    kw = {"sweeps": 1, "min_water_width_frac": None, "freeze_fixed_neighbours": False}
    # F1: fixed indices are whole and in range
    for fixed in ([5.9], [999], np.ma.array([5], mask=True), np.array([5 + 9j])):
        with pytest.raises(ValueError):
            fit_boundary_to_coast(nodes, tri, land, fixed=fixed, **kw)
    assert fit_boundary_to_coast(nodes, tri, land, fixed=[5], **kw).nodes[5].tolist() \
        == nodes[5].tolist()
    # F2: movement budget and relaxation
    for bad in ({"max_move_frac": -0.25}, {"max_move_frac": np.nan}, {"relax": 0.0},
                {"relax": np.nan}, {"sweeps": 0}):
        with pytest.raises(ValueError):
            fit_boundary_to_coast(nodes, tri, land, **{**kw, **bad})
    # F3: depths
    for dep in (np.ma.array(np.full(len(nodes), 10.0), mask=True), np.full(len(nodes), 10 + 9j),
                np.full(len(nodes), np.nan), np.full(3, 10.0)):
        with pytest.raises(ValueError):
            fit_boundary_to_coast(nodes, tri, land, depths=dep, **kw)
    tri1 = np.array([[0, 1, 2]])
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    m = Fort14Mesh(title="t", nodes=xy, elements=tri1, depths=np.array([5.0, 6.0, np.nan]),
                   open_boundaries=[np.array([0, 1])], land_boundaries=[])
    with pytest.raises(ValueError, match="finite"):
        apply_obc_depth_control(m)
    # F4: coincident OBC nodes have no normal
    with pytest.raises(ValueError, match="coincide|normal"):
        fvcom_next_obc([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0]], tri1, [0, 1])
    # F5: boundary chains are checked before the repair
    lat, tris, n = _lattice_mesh()
    for chain in (np.arange(6) + 0.9, np.ma.array(np.arange(6), mask=True)):
        mesh = Fort14Mesh(title="t", nodes=lat, elements=tris, depths=np.full(len(lat), 10.0),
                          open_boundaries=[chain], land_boundaries=[])
        with pytest.raises(ValueError):
            align_open_boundary_local(mesh, max_outer=1, n_jitter=1)
    # F6: no open boundary, nothing to control
    m0 = Fort14Mesh(title="t", nodes=xy, elements=tri1, depths=np.array([5.0, 6.0, 7.0]),
                    open_boundaries=[], land_boundaries=[])
    out, change = apply_obc_depth_control(m0)
    assert out.depths.tolist() == [5.0, 6.0, 7.0] and change.size == 0


def test_round41_guards():
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast
    from fvcom_mesh_tools.io.fvcom_native import apply_obc_depth_control, fvcom_next_obc
    from fvcom_mesh_tools.patch import improve_patch
    from tests.test_coast_fit import _land_above, _strip_mesh

    tri1 = np.array([[0, 1, 2]])
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    # F2: an undefined tangent
    m = Fort14Mesh(title="t", nodes=np.array([[0.0, 0], [0, 0], [1, 1]]), elements=tri1,
                   depths=np.ones(3), open_boundaries=[np.array([0, 1])], land_boundaries=[])
    with pytest.raises(ValueError, match="coincident"):
        align_open_boundary_local(m, max_outer=1, n_jitter=1)
    # F3: coast-fit quality controls
    nodes, tri = _strip_mesh(nx=3, ny=2)
    land = _land_above(1000, nx=3)
    kw = {"sweeps": 1, "min_water_width_frac": None}
    for bad in ({"max_area_change": np.nan}, {"min_angle_deg": np.nan},
                {"max_angle_deg": 10.0}, {"min_water_width_frac": 2.0}):
        with pytest.raises(ValueError):
            fit_boundary_to_coast(nodes, tri, land, **{**kw, **bad})
    # F4: patch scoring controls
    for bad in ({"max_area_change": np.nan}, {"max_valence": np.nan}, {"only_below": np.nan}):
        with pytest.raises(ValueError):
            improve_patch(xy, tri1, np.zeros(3, bool), np.zeros(1, bool), rounds=1, **bad)
    # F5: unsigned depths, no wrap in the reported change
    mu = Fort14Mesh(title="t", nodes=xy, elements=tri1,
                    depths=np.array([250, 250, 5], dtype=np.uint16),
                    open_boundaries=[np.array([0, 1])], land_boundaries=[])
    out, change = apply_obc_depth_control(mu)
    assert out.depths.tolist() == [5, 5, 5] and change.tolist() == [-245.0, -245.0]
    # F6: a third coordinate column is not part of the normal
    nxt, _ = fvcom_next_obc(np.c_[xy, np.zeros(3)], tri1, [0, 1])
    assert nxt.tolist() == [2, 2]
    lat, tris, n = _lattice_mesh()
    mf = Fort14Mesh(title="t", nodes=lat, elements=tris.astype(float),
                    depths=np.full(len(lat), 10.0), open_boundaries=[np.arange(n)],
                    land_boundaries=[])
    out, _ = align_open_boundary_local(mf, max_outer=1, n_jitter=1)
    assert out.elements.dtype.kind == "i"


def test_round42_guards():
    import shapely

    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast
    from fvcom_mesh_tools.io.fvcom_native import apply_obc_depth_control
    from fvcom_mesh_tools.patch import improve_patch

    tri1 = np.array([[0, 1, 2]])
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    # F2: a chain that doubles back
    m = Fort14Mesh(title="t", nodes=xy, elements=tri1, depths=np.ones(3),
                   open_boundaries=[np.array([0, 1, 0])], land_boundaries=[])
    with pytest.raises(ValueError, match="doubles back"):
        align_open_boundary_local(m, max_outer=1, n_jitter=1)
    # F3: a third coordinate column
    lat, tris, n = _lattice_mesh()
    m3 = Fort14Mesh(title="t", nodes=np.c_[lat, np.zeros(len(lat))], elements=tris,
                    depths=np.full(len(lat), 10.0), open_boundaries=[np.arange(n)],
                    land_boundaries=[])
    with pytest.raises(ValueError, match="N, 2"):
        align_open_boundary_local(m3, max_outer=1, n_jitter=1)
    # F4: extended-precision depths keep their change
    d0 = np.nextafter(np.longdouble(10), np.longdouble(np.inf))
    ml = Fort14Mesh(title="t", nodes=xy, elements=tri1,
                    depths=np.array([d0, 10, 10], dtype=np.longdouble),
                    open_boundaries=[np.array([0, 1])], land_boundaries=[])
    out, change = apply_obc_depth_control(ml)
    assert change[0] == out.depths[0] - ml.depths[0] != 0
    # F5: the angle endpoints divide the score
    for lo, hi in ((0, 130), (30, 180)):
        with pytest.raises(ValueError):
            improve_patch(xy, tri1, np.ones(3, bool), np.zeros(1, bool), rounds=1,
                          min_angle_deg=lo, max_angle_deg=hi)
    # F6: an unchanged fit reports its after-statistics
    res = fit_boundary_to_coast(xy, tri1, shapely.box(-1, 2, 2, 3), fixed=[0, 1, 2],
                                depths=np.ones(3))
    assert np.isfinite([res.min_angle_after_deg, res.max_angle_after_deg, res.dt_before_s,
                        res.dt_after_s]).all()


def test_round43_guards(tmp_path):
    from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
    from fvcom_mesh_tools.io.fvcom_native import apply_obc_depth_control, write_grd

    tri1 = np.array([[0, 1, 2]])
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

    def mesh_with(depths, nodes=xy):
        return Fort14Mesh(title="t", nodes=nodes, elements=tri1, depths=depths,
                          open_boundaries=[np.array([0, 1])], land_boundaries=[])

    # F1: object arrays are refused, large unsigned integers differ exactly
    obj = np.array([np.uint16(6), np.uint16(5), np.uint16(5)], dtype=object)
    with pytest.raises(ValueError, match="float or integer"):
        apply_obc_depth_control(mesh_with(obj))
    big = np.array([2**53 + 1, 2**53, 2**53], dtype=np.uint64)
    out, change = apply_obc_depth_control(mesh_with(big))
    assert change.tolist() == [-1.0, 0.0] and out.depths[0] == 2**53
    # F2: the violation count does not depend on the coordinate scale or dtype
    counts = []
    for scale in (1.0, 1e200):
        _, info = align_open_boundary_local(mesh_with(np.ones(3), xy * scale), max_outer=1,
                                            n_jitter=1)
        counts.append(info["remaining"])
    assert counts[0] == counts[1] and counts[0]
    lat, tris, n = _lattice_mesh()
    stored = lat.astype(np.float16)
    f16 = Fort14Mesh(title="t", nodes=stored, elements=tris, depths=np.full(len(lat), 10.0),
                     open_boundaries=[np.arange(n)], land_boundaries=[])
    f64 = Fort14Mesh(title="t", nodes=stored.astype(np.float64), elements=tris,
                     depths=np.full(len(lat), 10.0), open_boundaries=[np.arange(n)],
                     land_boundaries=[])
    kw = {"max_outer": 1, "n_jitter": 1, "movable": np.zeros(len(lat), bool)}
    rem16 = align_open_boundary_local(f16, **kw)[1]["remaining"]
    assert rem16 and rem16 == align_open_boundary_local(f64, **kw)[1]["remaining"]
    # F3: integer coordinates do not wrap in the area
    ints = Fort14Mesh(title="t", nodes=np.array([[0, 0], [0, 200], [200, 0]], dtype=np.int16),
                      elements=tri1, depths=np.ones(3), open_boundaries=[],
                      land_boundaries=[])
    with pytest.raises(ValueError, match="CCW"):
        write_grd(ints, tmp_path / "x.dat")
    assert not list(tmp_path.glob("x.*"))
