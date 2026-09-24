#!/usr/bin/env python3
"""Second-review regression expectations; intentionally fail on 18d7c85.

Run with pytest -q -p no:cacheprovider docs/hires_walls_tools_review_2_probes.py.
Only tiny temporary fixtures and extracted shell stages are executed.
No scheduler, mesher, MPI executable, or FVCOM binary is invoked.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import pytest

from fvcom_mesh_tools.cli.check_run import check_run, main
from fvcom_mesh_tools.dem import tokyo_bay as tb

ROOT = Path(__file__).resolve().parents[1]


def _run(tmp_path, fields=("zeta", "ua", "va")):
    run = tmp_path / "base"
    (run / "output").mkdir(parents=True)
    (run / "m2_run.nml").write_text("END_DATE = '2020-01-02 00:00:00',\n")
    (run / "fvcom.log").write_text("TADA!\n")
    with netCDF4.Dataset(run / "output/m2_0001.nc", "w") as ds:
        ds.createDimension("time", 2)
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", 3)
        times = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        for k, stamp in enumerate(("2020-01-01T00:00:00.000000",
                                   "2020-01-02T00:00:00.000000")):
            times[k] = np.asarray(list(stamp), dtype="S1")
        for field in fields:
            ds.createVariable(field, "f4", ("time", "node"))[:] = 0.1
    return run


@pytest.mark.parametrize("script,start,stop,marker", [
    ("421_finish_and_run.sh", '[ -f "$OUTDIR/ACCEPTED" ]',
     '# exits 3', "STAGED"),
    ("423_m2_smoke.sh", '[ -f "$RUN_ROOT/STAGED" ]',
     'python - "$RUN_ROOT"', "SMOKE_OK"),
    ("412_m2_run.sh", '[ -f "$RUN_ROOT/$CASE/m2_run.nml" ]',
     'set +u', "base/RUN_OK"),
])
def test_r1_failed_prerequisite_invalidates_old_marker(tmp_path, script, start, stop, marker):
    """Exercise real prerequisite/cleanup statements, before any heavy work."""
    old = tmp_path / marker
    old.parent.mkdir(parents=True, exist_ok=True)
    old.write_text("previous successful attempt\n")
    source = (ROOT / "jobs/octopus" / script).read_text()
    block = source[source.index(start):source.index(stop, source.index(start))]
    env = {**os.environ, "OUTDIR": str(tmp_path / "unaccepted"),
           "RUN_ROOT": str(tmp_path), "CASE": "base", "FMESH_REQUIRE_SMOKE": "1"}
    result = subprocess.run(["bash", "-c", "set -euo pipefail\n" + block],
                            env=env, capture_output=True, text=True, check=False)
    assert result.returncode == 2, result.stdout + result.stderr
    assert not old.exists(), f"{script} exited 2 but left {marker}"


def test_r1_checker_failure_invalidates_old_marker(tmp_path):
    run = _run(tmp_path)
    marker = run / "RUN_OK"
    assert main([str(run), "--marker", str(marker)]) == 0
    (run / "fvcom.log").write_text("STOP: integration ended early\n")
    assert main([str(run), "--marker", str(marker)]) == 1
    assert not marker.exists(), "Failed recheck left a success marker"


def test_r2_nonzero_solver_exit_cannot_write_run_ok(tmp_path):
    """Actual 412 tail, stub only MPI and environment modules; real checker."""
    run = _run(tmp_path)
    source = (ROOT / "jobs/octopus/412_m2_run.sh").read_text()
    block = source[source.index("status=0\n"):]
    stubs = """set -euo pipefail
module() { :; }
conda() { :; }
mpiexec() { echo 'TADA!'; return 7; }
"""
    result = subprocess.run(
        ["bash", "-c", stubs + block],
        env={**os.environ, "RUN_ROOT": str(tmp_path), "CASE": "base",
             "RANKS": "1", "FVCOM": "unused-by-stub"},
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 1, result.stdout + result.stderr
    assert not (run / "RUN_OK").exists(), result.stdout


@pytest.mark.parametrize("fields", [(), ("zeta",), ("zeta", "ua")])
def test_r3_required_output_fields_must_exist(tmp_path, fields):
    run = _run(tmp_path, fields)
    verdict = check_run(run)
    assert not verdict["ok"], verdict


def test_r4_large_gap_cannot_replace_configured_output_interval(tmp_path):
    run = _run(tmp_path)
    (run / "m2_run.nml").write_text(
        "END_DATE = '2020-01-03 00:00:00',\n"
        "NC_OUT_INTERVAL = 'seconds=3600.',\n"
    )
    verdict = check_run(run)
    assert not verdict["ok"], verdict  # last output is 24 hourly records short


def test_r5_depth_choice_is_independent_of_batch_latitude(monkeypatch):
    """Near a distance tie within Tokyo Bay, batching changes the chosen depth."""
    monkeypatch.setattr(tb, "_interp", lambda rung, lon, lat: np.full(len(lon), np.nan))
    monkeypatch.setattr(tb, "_grid", lambda rung: (
        np.array([35.0, 35.00818]), np.array([139.0, 139.01]),
        np.array([[np.nan, -10.0], [-20.0, np.nan]]),
    ))
    one = tb.sample([139.0], [35.0])[0][0]
    batched = tb.sample([139.0, 139.0], [35.0, 35.6])[0][0]
    assert one == batched, (one, batched)


def test_r5_sounding_distance_is_independent_of_batch_latitude(tmp_path, monkeypatch):
    monkeypatch.setattr(tb, "_data_dir", lambda: tmp_path)
    monkeypatch.setattr(tb, "_SOUNDINGS", "soundings.parquet")
    (tmp_path / "soundings.parquet").touch()
    monkeypatch.setattr(pd, "read_parquet", lambda *a, **kw: pd.DataFrame({
        "lon": [139.01], "lat": [35.0], "z_tp": [-10.0],
    }))
    one = tb.sounding_distance([139.0], [35.0])[0]
    batched = tb.sounding_distance([139.0, 139.0], [35.0, 35.6])[0]
    assert one == pytest.approx(batched, abs=1e-6), (one, batched)
