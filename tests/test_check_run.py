"""fmesh-check-run: a zero exit code is not an FVCOM run that finished."""

from __future__ import annotations

import json

import netCDF4
import numpy as np

from fvcom_mesh_tools.cli.check_run import check_run, main


def _run(tmp_path, *, tada=True, times=("2020-01-01T00:00:00.000000",
                                         "2020-01-02T00:00:00.000000"),
         zeta=0.1, end="2020-01-02 00:00:00", extra_log=""):
    run = tmp_path / "run"
    (run / "output").mkdir(parents=True)
    (run / "m2_run.nml").write_text(f" END_DATE = '{end}',\n")
    (run / "fvcom.log").write_text("step ...\n" + extra_log + ("TADA!\n" if tada else ""))
    with netCDF4.Dataset(run / "output" / "m2_0001.nc", "w") as ds:
        ds.createDimension("time", None)
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", 3)
        t = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        z = ds.createVariable("zeta", "f4", ("time", "node"))
        for k, s in enumerate(times):
            t[k] = np.array(list(s.ljust(26)), dtype="S1")
            z[k] = np.full(3, zeta)
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
