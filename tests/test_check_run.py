"""fmesh-check-run: a zero exit code is not an FVCOM run that finished."""

from __future__ import annotations

import json

import netCDF4
import numpy as np
import pytest

from fvcom_mesh_tools.cli.check_run import check_run, main


def _run(tmp_path, *, tada=True, times=("2020-01-01T00:00:00.000000",
                                         "2020-01-02T00:00:00.000000"),
         zeta=0.1, end="2020-01-02 00:00:00", extra_log="", n=3, grid=None):
    run = tmp_path / "run"
    (run / "output").mkdir(parents=True)
    (run / "m2_run.nml").write_text(" START_DATE = '2020-01-01 00:00:00',\n"
                                    f" END_DATE = '{end}',\n"
                                    " NC_OUT_INTERVAL = 'seconds = 86400.0',\n"
                                    + (f" INPUT_DIR = '{run}/input/',\n"
                                       " GRID_FILE = 'm2_grd.dat',\n" if grid else ""))
    if grid:
        # a whole grid: n_nodes points, every cell on nodes 1, 2, 3
        (run / "input").mkdir()
        (run / "input" / "m2_grd.dat").write_text(
            f"Node Number = {grid[0]}\nCell Number = {grid[1]}\n"
            + "".join(f"{k} 1 2 3\n" for k in range(1, grid[1] + 1))
            + "".join(f"{k} {1000.0 * k} {500.0 * (k % 2)}\n" for k in range(1, grid[0] + 1)))
    (run / "fvcom.log").write_text("step ...\n" + extra_log + ("TADA!\n" if tada else ""))
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", n)
        ds.createDimension("nele", n)
        t = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        z = ds.createVariable("zeta", "f4", ("time", "node"))
        u = ds.createVariable("ua", "f4", ("time", "nele"))
        v = ds.createVariable("va", "f4", ("time", "nele"))
        # the grid the run was on (FVCOM's x, y, nv), matching the fixture's
        ds.createDimension("three", 3)
        ds.createVariable("x", "f8", ("node",))[:] = [1000.0 * (k + 1) for k in range(n)]
        ds.createVariable("y", "f8", ("node",))[:] = [500.0 * ((k + 1) % 2) for k in range(n)]
        ds.createVariable("nv", "i4", ("three", "nele"))[:] = np.tile([[3], [2], [1]], (1, n))
        for k, s in enumerate(times):
            t[k] = np.array(list(s.ljust(26)), dtype="S1")
            z[k] = np.full(n, zeta)
            u[k] = np.zeros(n)
            v[k] = np.zeros(n)
    return run


def test_a_finished_run_passes_and_writes_its_marker(tmp_path):
    run = _run(tmp_path)
    assert check_run(run)["ok"]
    marker = tmp_path / "RUN_OK"
    assert main([str(run), "--marker", str(marker)]) == 0
    assert json.loads(marker.read_text())["ok"]


def test_no_tada_fails_even_with_a_clean_exit(tmp_path):
    run = _run(tmp_path, tada=False, extra_log="STOP: integration ended early\n")
    info = check_run(run)
    assert not info["ok"] and any("TADA" in r for r in info["reasons"])
    marker = tmp_path / "RUN_OK"
    assert main([str(run), "--marker", str(marker)]) == 1 and not marker.exists()


def test_output_that_stops_short_of_end_date_fails(tmp_path):
    run = _run(tmp_path, end="2020-01-05 00:00:00")
    assert any("before END_DATE" in r for r in check_run(run)["reasons"])


def test_non_finite_output_fails(tmp_path):
    run = _run(tmp_path, zeta=np.nan)
    assert any("zeta" in r for r in check_run(run)["reasons"])


def test_an_unreadable_placeholder_fails(tmp_path):
    run = _run(tmp_path)
    (run / "output" / "m2_0002.nc").write_bytes(b"not a completed NetCDF run")
    assert any("cannot be read" in r for r in check_run(run)["reasons"])


def test_a_fatal_word_in_the_log_fails(tmp_path):
    run = _run(tmp_path, extra_log="forrtl: severe (174): SIGSEGV, segmentation fault\n")
    assert any("segmentation" in r for r in check_run(run)["reasons"])


def test_a_gap_in_the_output_fails_against_the_declared_interval(tmp_path):
    run = _run(tmp_path, times=("2020-01-01T00:00:00.000000", "2020-01-03T00:00:00.000000"),
               end="2020-01-03 00:00:00")
    assert any("gap" in r for r in check_run(run)["reasons"])


def test_the_history_output_must_carry_ua_and_va(tmp_path):
    run = _run(tmp_path)
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds.renameVariable("ua", "u_other")
    assert any("lacks ua" in r for r in check_run(run)["reasons"])


def test_a_failed_recheck_removes_the_old_marker(tmp_path):
    run = _run(tmp_path)
    marker = tmp_path / "RUN_OK"
    assert main([str(run), "--marker", str(marker)]) == 0 and marker.exists()
    (run / "fvcom.log").write_text("STOP: integration ended early\n")
    assert main([str(run), "--marker", str(marker)]) == 1 and not marker.exists()


def test_fortran_reals_and_cycles_are_read_as_fvcom_does():
    from fvcom_mesh_tools.cli.check_run import _fortran_float, _interval

    assert _fortran_float("1.") == 1.0 and _fortran_float("1.5d0") == 1.5
    nml = " EXTSTEP_SECONDS = 2.0d0, ISPLIT = 10,\n"
    assert _interval("cycles = 180", nml).total_seconds() == 3600.0
    assert _interval("days = 0.5", nml).total_seconds() == 43200.0


def test_comments_are_stripped_outside_quotes_only():
    from fvcom_mesh_tools.cli.check_run import _nml_value, _strip_comments

    text = (" ! NC_OUT_INTERVAL = 'days = 1',\n"
            " NC_OUT_INTERVAL = 'seconds = 3600.', ! hourly\n"
            " CASE_TITLE = 'a ! inside quotes',\n")
    clean = _strip_comments(text)
    assert _nml_value(clean, "NC_OUT_INTERVAL") == "seconds = 3600."
    assert _nml_value(clean, "CASE_TITLE") == "a ! inside quotes"


def test_history_without_spatial_data_fails(tmp_path):
    """Review of the extend tools, round 4 F14: (time, 0) arrays passed."""
    info = check_run(_run(tmp_path, n=0))
    assert not info["ok"] and any("shape" in r for r in info["reasons"])


def test_history_must_be_on_the_staged_mesh(tmp_path):
    assert check_run(_run(tmp_path / "a", grid=(3, 3)))["ok"]
    info = check_run(_run(tmp_path / "b", grid=(4, 3)))
    assert not info["ok"] and any("staged mesh" in r for r in info["reasons"])


def test_a_named_grid_that_cannot_be_read_fails(tmp_path):
    """Review round 5 F6."""
    run = _run(tmp_path, grid=(3, 3))
    (run / "input" / "m2_grd.dat").unlink()
    info = check_run(run)
    assert not info["ok"] and any("cannot be read" in r for r in info["reasons"])


def test_the_history_is_looked_for_in_the_namelists_output_dir(tmp_path):
    """Review round 6 F3: an old history in run/output passed although the
    namelist sends the output elsewhere."""
    run = _run(tmp_path)
    nml = run / "m2_run.nml"
    nml.write_text(nml.read_text() + " OUTPUT_DIR = 'current_output/',\n")
    (run / "current_output").mkdir()
    info = check_run(run)
    assert not info["ok"] and info["output_dir"].endswith("current_output")
    for f in (run / "output").iterdir():
        f.rename(run / "current_output" / f.name)
    assert check_run(run)["ok"]


def test_quoted_values_round_trip_with_the_namelist_writer():
    """Review round 7 F3: OUTPUT_DIR 'a''b/' was read as 'a'."""
    from fvcom_mesh_tools.cli.check_run import _nml_value
    from fvcom_mesh_tools.io.fvcom_namelist import fortran_string

    for d in ("/x/current'case/", "/x/plain/"):
        assert _nml_value(f" OUTPUT_DIR = {fortran_string(d)},\n", "OUTPUT_DIR") == d
    assert _nml_value(' OUTPUT_DIR = "/x/a""b/",\n', "OUTPUT_DIR") == '/x/a"b/'


def test_an_assignment_inside_another_value_is_not_read():
    """Review round 8 F8."""
    import pytest

    from fvcom_mesh_tools.cli.check_run import _nml_value

    text = ("&NML_CASE\n CASE_TITLE = \"OUTPUT_DIR = '/x/stale/'\",\n/\n"
            "&NML_IO\n OUTPUT_DIR = '/x/current/',\n/\n")
    assert _nml_value(text, "OUTPUT_DIR") == "/x/current/"
    assert _nml_value(" A = 'it''s OUTPUT_DIR = 1',\n OUTPUT_DIR = 'b/',\n", "OUTPUT_DIR") == "b/"
    with pytest.raises(ValueError, match="more than once"):
        _nml_value(" OUTPUT_DIR = 'a/',\n OUTPUT_DIR = 'b/',\n", "OUTPUT_DIR")


def test_a_non_finite_three_dimensional_field_fails(tmp_path):
    """Review round 10 F4: NaN u and infinite w passed beside finite zeta."""
    run = _run(tmp_path)
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds.createDimension("siglay", 2)
        ds.createVariable("u", "f4", ("time", "siglay", "node"))[:] = np.nan
        ds.createVariable("w", "f4", ("time", "siglay", "node"))[:] = np.inf
        ds.createVariable("iint", "i4", ("time",))[:] = [1, 2]
    info = check_run(run)
    assert not info["ok"]
    assert any(" u is not finite" in r for r in info["reasons"])
    assert any(" w is not finite" in r for r in info["reasons"])


def test_an_empty_three_dimensional_field_fails(tmp_path):
    """Review round 11 F4."""
    run = _run(tmp_path)
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds.createDimension("siglay", 0)
        ds.createVariable("u", "f4", ("time", "siglay", "node"))
    info = check_run(run)
    assert not info["ok"] and any("empty" in r for r in info["reasons"])


def test_time_not_first_and_a_zero_length_run_fail(tmp_path):
    """Review round 12 F4 and F5."""
    run = _run(tmp_path / "a")
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds.createDimension("siglay", 2)
        ds.createVariable("u", "f4", ("siglay", "time", "node"))[:] = np.nan
    info = check_run(run)
    assert not info["ok"] and any("time at position" in r for r in info["reasons"])
    run = _run(tmp_path / "b", times=("2020-01-01T00:00:00.000000",), end="2020-01-01 00:00:00")
    info = check_run(run)
    assert not info["ok"] and any("not after the start" in r for r in info["reasons"])


def test_transposed_barotropic_fields_fail(tmp_path):
    """Review round 13 F4: zeta(node, time) of a square shape passed."""
    run = _run(tmp_path, n=2, times=("2020-01-01T00:00:00.000000",
                                     "2020-01-02T00:00:00.000000"))
    f = run / "output" / "m2_0001.nc"
    f.unlink()
    with netCDF4.Dataset(f, "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", 2)
        ds.createDimension("nele", 2)
        t = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        for k, s in enumerate(("2020-01-01T00:00:00.000000", "2020-01-02T00:00:00.000000")):
            t[k] = np.array(list(s.ljust(26)), dtype="S1")
        ds.createVariable("zeta", "f4", ("node", "time"))[:] = 0.1
        ds.createVariable("ua", "f4", ("time", "nele"))[:] = 0.0
        ds.createVariable("va", "f4", ("time", "nele"))[:] = 0.0
    info = check_run(run)
    assert not info["ok"] and any("zeta has dimensions" in r for r in info["reasons"])


def test_an_initial_record_alone_is_not_a_run(tmp_path):
    """Review round 13 F5: a 10-minute run with only the 00:00 record."""
    run = _run(tmp_path, times=("2020-01-01T00:00:00.000000",), end="2020-01-01 00:10:00")
    info = check_run(run)
    assert not info["ok"] and any("no output after the start" in r for r in info["reasons"])


def test_sub_second_records_are_distinct(tmp_path):
    """Review round 14 F7."""
    times = tuple(f"2020-01-01T00:00:00.{us:06d}" for us in (0, 250000, 500000, 750000)) + (
        "2020-01-01T00:00:01.000000",)
    run = _run(tmp_path, times=times, end="2020-01-01 00:00:01")
    nml = run / "m2_run.nml"
    nml.write_text(nml.read_text().replace("seconds = 86400.0", "seconds = 0.25"))
    info = check_run(run)
    assert not any("do not increase" in r for r in info["reasons"]), info["reasons"]


def test_history_from_another_mesh_of_the_same_size_fails(tmp_path):
    """Review round 19 F7: counts matched, coordinates did not."""
    run = _run(tmp_path, grid=(3, 3))
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds["x"][:] = ds["x"][:] + 100000.0
    info = check_run(run)
    assert not info["ok"] and any("x is not the staged" in r for r in info["reasons"])


def _stage_depths(run, n, h=5.0):
    nml = run / "m2_run.nml"
    nml.write_text(nml.read_text() + " DEPTH_FILE = 'm2_dep.dat',\n")
    (run / "input" / "m2_dep.dat").write_text(
        f"Node Number = {n}\n" + "".join(f"{1000.0 * k} {500.0 * (k % 2)} {h}\n"
                                         for k in range(1, n + 1)))


def test_identity_is_checked_on_every_stack_with_depths_masks_and_float32(tmp_path):
    """Review round 20 F2-F5."""
    import shutil

    # F3: h from another bathymetry
    run = _run(tmp_path / "a", grid=(3, 3))
    _stage_depths(run, 3)
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds.createVariable("h", "f8", ("node",))[:] = 500.0
    info = check_run(run)
    assert not info["ok"] and any("another bathymetry" in r for r in info["reasons"])
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        ds["h"][:] = 5.0
    assert check_run(run)["ok"]
    # F2: a second stack on another mesh
    second = run / "output" / "m2_0002.nc"
    shutil.copy(run / "output" / "m2_0001.nc", second)
    with netCDF4.Dataset(second, "a") as ds:
        ds["x"][:] = ds["x"][:] + 100000.0
    info = check_run(run)
    assert not info["ok"] and any(r.startswith("m2_0002.nc: x") for r in info["reasons"])
    second.unlink()
    # F4: a masked coordinate
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "a") as ds:
        x = ds["x"][:]
        ds["x"][:] = np.ma.masked_array(x, mask=[True, False, False])
    info = check_run(run)
    assert not info["ok"] and any("masked" in r for r in info["reasons"])


def test_float32_history_coordinates_are_accepted(tmp_path):
    """Review round 20 F5: UTM coordinates stored as float32."""
    run = _run(tmp_path, grid=(3, 3))
    g = run / "input" / "m2_grd.dat"
    xs = [380000.123456789 + 1000 * k for k in range(1, 4)]
    ys = [3900000.123456789 + 500 * (k % 2) for k in range(1, 4)]
    lines = g.read_text().splitlines()
    lines[-3:] = [f"{k} {x!r} {y!r}" for k, (x, y) in enumerate(zip(xs, ys), 1)]
    g.write_text("\n".join(lines) + "\n")
    f = run / "output" / "m2_0001.nc"
    with netCDF4.Dataset(f, "a") as ds:
        ds.renameVariable("x", "x64")
        ds.renameVariable("y", "y64")
        ds.createVariable("x", "f4", ("node",))[:] = np.float32(xs)
        ds.createVariable("y", "f4", ("node",))[:] = np.float32(ys)
    assert check_run(run)["ok"], check_run(run)["reasons"]


def test_history_past_end_date_and_bad_timing_controls_fail(tmp_path):
    """Review round 21 F1 and F9."""
    from fvcom_mesh_tools.cli.check_run import _interval

    days = tuple(f"2020-01-0{d}T00:00:00.000000" for d in range(1, 6))
    info = check_run(_run(tmp_path / "a", times=days, end="2020-01-02 00:00:00"))
    assert not info["ok"] and any("past END_DATE" in r for r in info["reasons"])
    info = check_run(_run(tmp_path / "b", end="garbage"))
    assert not info["ok"] and any("END_DATE" in r for r in info["reasons"])
    for v in ("seconds = 0", "seconds = 1e999"):
        with pytest.raises(ValueError):
            _interval(v, "")


def test_start_date_is_checked_on_its_own_and_cycle_controls_too(tmp_path):
    """Review round 22 F2 and F4."""
    from fvcom_mesh_tools.cli.check_run import _interval

    run = _run(tmp_path)
    nml = run / "m2_run.nml"
    nml.write_text(nml.read_text().replace("2020-01-01 00:00:00", "garbage")
                   + " NC_FIRST_OUT = '2020-01-01 00:00:00',\n")
    info = check_run(run)
    assert not info["ok"] and any("START_DATE" in r for r in info["reasons"])
    for text in (" EXTSTEP_SECONDS = 1, ISPLIT = 1e999,", " EXTSTEP_SECONDS = -1, ISPLIT = -10,",
                 " EXTSTEP_SECONDS = 1, ISPLIT = 1.9,"):
        with pytest.raises(ValueError):
            _interval("cycles = 1", text)
    with pytest.raises(ValueError):
        _interval("seconds = 1e-9", "")
