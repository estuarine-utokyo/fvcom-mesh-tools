"""fmesh-check-run: a zero exit code is not an FVCOM run that finished."""

from __future__ import annotations

import json

import netCDF4
import numpy as np

from fvcom_mesh_tools.cli.check_run import check_run, main


def _run(tmp_path, *, tada=True, times=("2020-01-01T00:00:00.000000",
                                         "2020-01-02T00:00:00.000000"),
         zeta=0.1, end="2020-01-02 00:00:00", extra_log="", n=3, grid=None):
    run = tmp_path / "run"
    (run / "output").mkdir(parents=True)
    (run / "m2_run.nml").write_text(f" END_DATE = '{end}',\n"
                                    " NC_OUT_INTERVAL = 'seconds = 86400.0',\n"
                                    + (f" INPUT_DIR = '{run}/input/',\n"
                                       " GRID_FILE = 'm2_grd.dat',\n" if grid else ""))
    if grid:
        (run / "input").mkdir()
        (run / "input" / "m2_grd.dat").write_text(
            f"Node Number = {grid[0]}\nCell Number = {grid[1]}\n")
    (run / "fvcom.log").write_text("step ...\n" + extra_log + ("TADA!\n" if tada else ""))
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", n)
        t = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        z = ds.createVariable("zeta", "f4", ("time", "node"))
        u = ds.createVariable("ua", "f4", ("time", "node"))
        v = ds.createVariable("va", "f4", ("time", "node"))
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
