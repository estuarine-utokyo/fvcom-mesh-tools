"""The review findings of docs/hires_walls_tools_review.md, as regression tests.

Written as failing probes by the reviewer (gpt-6-astra); kept here, adapted
where the fix changed the structure they probed, and passing since the fixes.

Run: PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider \
    docs/hires_walls_tools_review_probes.py
Small synthetic inputs only: no meshing, scheduler submission or FVCOM run.
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

import netCDF4
import numpy as np
import pytest
import shapely

from fvcom_mesh_tools.cli.finish_depths import main as finish_main
from fvcom_mesh_tools.cli.refine_run import main as refine_main
from fvcom_mesh_tools.dem import tokyo_bay as tb
from fvcom_mesh_tools.io.fort14 import Fort14Mesh
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.patch import _corner_walk, blunt_acute_corners, hole_polygon

ROOT = Path(__file__).resolve().parents[1]


def test_f1_blunting_adjacent_corners_has_no_deleted_endpoint():
    """F1: two adjacent 45-degree corners leave -1 edge endpoints."""
    xy = np.array([[0.0, 0.0], [100.0, 0.0], [50.0, 50.0]])
    edges = np.array([[0, 1], [1, 2], [2, 0]])
    _, new_edges, _, report = blunt_acute_corners(
        xy, edges, [-1, -1, -1], shapely.Polygon(xy),
        lambda q: np.full(len(q), 10.0),
    )
    assert report["n_corners_blunted"] == 2
    assert (new_edges >= 0).all(), new_edges


def test_f2_island_corner_walk_is_translation_invariant():
    """F2: default allclose silently drops an unclosed ring's last UTM vertex."""
    xy = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 30.0], [0.0, 30.0]])
    origin = np.array([392000.0, 3909000.0])
    local = _corner_walk(xy, 10.0, closed=True)
    utm = _corner_walk(xy + origin, 10.0, closed=True) - origin
    assert shapely.Polygon(utm).area == pytest.approx(shapely.Polygon(local).area)


def test_f3_netcdf_fill_value_falls_through_to_next_rung(tmp_path, monkeypatch):
    """F3: a finite _FillValue must stay missing when the NetCDF mask is decoded."""
    nc = tmp_path / "masked.nc"
    with netCDF4.Dataset(nc, "w") as ds:
        ds.createDimension("lat", 2)
        ds.createDimension("lon", 2)
        ds.createVariable("lat", "f8", ("lat",))[:] = [35.0, 36.0]
        ds.createVariable("lon", "f8", ("lon",))[:] = [139.0, 140.0]
        z = ds.createVariable("elevation", "f8", ("lat", "lon"), fill_value=-9999.0)
        z[:] = np.full((2, 2), -9999.0)
    monkeypatch.setenv("DATA_DIR", str(tmp_path))
    monkeypatch.setitem(tb._PRODUCTS, "m7001", (nc.name, "elevation"))
    tb._grid.cache_clear()
    try:
        _, _, z = tb._grid("m7001")
        assert np.isnan(z).all(), z
    finally:
        tb._grid.cache_clear()


def test_f4_hires_preserve_does_not_blunt_the_coast():
    """F4: execute the actual driver's blunting block with preserve selected."""
    source = (ROOT / "notebooks/420_local_refine.py").read_text()
    block = next(
        n for n in ast.parse(source).body
        if isinstance(n, ast.If)
        and any(isinstance(c, ast.Call) and isinstance(c.func, ast.Name)
                and c.func.id == "blunt_acute_corners" for c in ast.walk(n))
    )
    xy = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 30.0]])
    water = shapely.Polygon(xy)
    env = {
        "HIRES": {"coastline": "preserve"}, "cfg": {"coastline": "preserve"},
        "rc": {"pfix": xy, "egfix": np.array([[0, 1], [1, 2], [2, 0]]),
               "pfix_base": np.full(3, -1)},
        "hole": water, "h_achieved": lambda q: np.full(len(q), 10.0),
        "blunt_acute_corners": blunt_acute_corners, "hole_polygon": hole_polygon,
        "reports": {}, "say": lambda *args: None,
    }
    exec(compile(ast.Module(body=[block], type_ignores=[]), "<driver block>", "exec"), env)
    assert env["hole"].symmetric_difference(water).area < 1e-8


def _case(tmp_path, depths=(3.0, 100.0, 100.0)):
    mesh = Fort14Mesh(
        title="probe", nodes=np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]]),
        elements=np.array([[0, 1, 2]]), depths=np.asarray(depths),
        open_boundaries=[], land_boundaries=[],
    )
    export_fvcom_case(mesh, tmp_path, "probe", twodm=False, obc_depth_control=False)
    return tmp_path / "probe_grd.dat"


@pytest.mark.parametrize("patch", [False, True])
def test_f7_finisher_rejects_unconverged_product(tmp_path, patch):
    """F7: both exhausted sweeps and infeasible frozen seams return success."""
    if patch:
        _case(tmp_path / "fvcom", depths=(1.0, 3.0, 3.0))
        np.save(tmp_path / "node_map.npy", np.array([0]))
        args = [str(tmp_path), "--hmin", "3"]
    else:
        grd = _case(tmp_path)
        args = [str(grd), "--hmin", "3", "--rounds", "0"]
    assert finish_main(args) != 0


@pytest.mark.parametrize("value", ["nan", "-0.2", "1.2", "inf"])
def test_f8_finisher_rejects_invalid_rfactor(tmp_path, value):
    """F8: only zero (disabled) or a finite r in (0,1) is meaningful."""
    grd = _case(tmp_path)
    assert finish_main([str(grd), "--hmin", "3", "--rfactor", value, "--dry-run"]) == 2


def test_f9_nearest_finite_search_is_independent_of_query_batch(monkeypatch):
    """F9: a point outside the first search box can be the nearest cell."""
    lat = np.array([0.0, 0.019])
    lon = np.array([0.019, 0.021])
    z = np.array([[np.nan, -20.0], [-10.0, np.nan]])
    monkeypatch.setattr(tb, "_grid", lambda rung: (lat, lon, z))
    one, distance = tb._nearest_finite("m7001", np.array([0.0]), np.array([0.0]), 0.0)
    batch, _ = tb._nearest_finite("m7001", np.array([0.0, 0.01]), np.array([0.0, 0.0]), 0.0)
    assert one[0] == batch[0] == -20.0, (one, batch, distance)


def test_f10_refine_refuses_existing_output_without_report(tmp_path, monkeypatch):
    """F10: a failed pre-report run's directory is still an existing output."""
    recipe = tmp_path / "p.yaml"
    recipe.write_text("hires: {}\n")
    land = tmp_path / "land.shp"
    land.touch()
    out = tmp_path / "existing"
    out.mkdir()
    (out / "fill_constraints.npz").write_bytes(b"prior result")
    calls = []
    monkeypatch.setattr(subprocess, "call", lambda *a, **kw: calls.append(a) or 0)
    rc = refine_main([str(recipe), "--out", str(out), "--land", str(land)])
    assert rc == 2 and not calls


def test_f11_the_land_filter_does_not_wait_for_a_base_coastline():
    """F11: the filter (islands, walls) runs on the hires/resolve branch even
    when the hole has no free base coastline -- it may not be conditioned on
    `shore`, which only exists where there is a base coastline to re-cut."""
    source = (ROOT / "notebooks/420_local_refine.py").read_text()
    tree = ast.parse(source)
    block = next(
        n for n in tree.body
        if isinstance(n, ast.If)
        and any(isinstance(c, ast.Call) and isinstance(c.func, ast.Name)
                and c.func.id == "filter_shoreline_local" for c in ast.walk(n)))
    names = {x.id for x in ast.walk(block.test) if isinstance(x, ast.Name)}
    assert "shore" not in names and "_keep" in names
    isl = next(
        n for n in tree.body
        if isinstance(n, ast.If)
        and any(isinstance(c, ast.Call) and isinstance(c.func, ast.Name)
                and c.func.id == "island_rings" for c in ast.walk(n)))
    assert {x.id for x in ast.walk(isl.test) if isinstance(x, ast.Name)} >= {"_land_filtered"}
    from fvcom_mesh_tools.patch import island_rings

    water = shapely.box(0, 0, 1000, 1000)
    land = shapely.box(400, 400, 600, 600)
    env = {"HIRES": {"coastline": "resolve"}, "_land_filtered": True, "_filtered": land,
           "hole": water, "reports": {}, "island_rings": island_rings,
           "h_achieved": lambda q: np.full(len(q), 30.0), "hole_polygon": hole_polygon,
           "rc": {"pfix": np.zeros((0, 2)), "egfix": np.zeros((0, 2), dtype=np.int64),
                  "pfix_base": np.zeros(0, dtype=np.int64), "curves": []},
           "np": np, "say": lambda *a: None}
    try:
        exec(compile(ast.Module(body=[isl], type_ignores=[]), "<driver block>", "exec"), env)
    except ValueError:
        pass            # hole_polygon of a lone island ring; the count is what matters
    assert env["reports"]["islands_added"]["n_islands_added"] == 1


def test_f14_the_workflow_refuses_a_comma_before_submitting(tmp_path):
    recipe = tmp_path / "a,b.yaml"
    recipe.write_text("refine: []\n")
    r = subprocess.run(["bash", str(ROOT / "jobs/octopus/refine_workflow.sh"), str(recipe)],
                       text=True, capture_output=True, check=False)
    assert r.returncode == 2 and "comma" in r.stdout


@pytest.mark.parametrize("script,marker", [
    ("421_finish_and_run.sh", "ACCEPTED"), ("423_m2_smoke.sh", "STAGED"),
    ("412_m2_run.sh", "SMOKE_OK"), ("413_m2_analysis.sh", "RUN_OK")])
def test_f5_every_stage_requires_the_previous_stage_s_marker(script, marker):
    text = (ROOT / "jobs/octopus" / script).read_text()
    assert marker in text and "exit 2" in text
