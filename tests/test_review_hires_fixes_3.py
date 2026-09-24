"""The third review's findings (docs/hires_walls_tools_review_3.md) as regression
tests: written as failing probes by the reviewer (gpt-6-astra), passing since
the fixes. No mesher, solver, scheduler or laboratory data."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path
from threading import Barrier, BrokenBarrierError, Lock

import netCDF4
import numpy as np
import pytest

from fvcom_mesh_tools.cli import refine_run
from fvcom_mesh_tools.cli.check_run import check_run
from fvcom_mesh_tools.cli.refine_depths import finish_depths


def _run(tmp_path, interval="seconds = 3600."):
    run = tmp_path / "run"
    (run / "output").mkdir(parents=True)
    (run / "fvcom.log").write_text("TADA!\n")
    (run / "m2_run.nml").write_text(
        "&NML_CASE\n START_DATE = '2020-01-01 00:00:00',\n"
        " END_DATE = '2020-01-02 00:00:00',\n/\n"
        "&NML_INTEGRATION\n EXTSTEP_SECONDS = 1., ISPLIT = 10,\n/\n"
        "&NML_NETCDF\n NC_ON = T, NC_OUTPUT_STACK = 12,\n"
        " NC_FIRST_OUT = '2020-01-01 00:00:00',\n"
        f" NC_OUT_INTERVAL = '{interval}',\n/\n"
    )
    return run


def _history(run, name, hours):
    with netCDF4.Dataset(run / "output" / name, "w") as ds:
        ds.createDimension("time", len(hours))
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", 3)
        ds.createDimension("nele", 1)
        t = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        for i, hour in enumerate(hours):
            stamp = datetime(2020, 1, 1) + timedelta(hours=hour)
            t[i] = np.array(list(stamp.strftime("%Y-%m-%dT%H:%M:%S.%f")), dtype="S1")
        for name, dim in (("zeta", "node"), ("ua", "nele"), ("va", "nele")):
            ds.createVariable(name, "f4", ("time", dim))[:] = 0.1


def test_t1_restart_is_not_history(tmp_path):
    run = _run(tmp_path)
    _history(run, "m2_restart_0001.nc", [24])
    verdict = check_run(run)
    assert not verdict["ok"], verdict


def test_t1_restart_must_not_fill_history_gap(tmp_path):
    run = _run(tmp_path)
    _history(run, "m2_0001.nc", range(12))
    _history(run, "m2_0003.nc", [24])
    for hour in range(12, 24):
        _history(run, f"m2_restart_{hour:04}.nc", [hour])
    verdict = check_run(run)
    assert not verdict["ok"], verdict


@pytest.mark.parametrize("defect", ["missing_first_stack", "reversed_stack"])
def test_t2_history_requires_start_coverage_and_order(tmp_path, defect):
    run = _run(tmp_path)
    if defect == "missing_first_stack":
        _history(run, "m2_0002.nc", range(12, 24))
        _history(run, "m2_0003.nc", [24])
    else:
        _history(run, "m2_0001.nc", list(range(12))[::-1])
        _history(run, "m2_0002.nc", range(12, 24))
        _history(run, "m2_0003.nc", [24])
    verdict = check_run(run)
    assert not verdict["ok"], verdict


@pytest.mark.parametrize("syntax", ["cycles", "double_quotes"])
def test_t3_unparsed_cadence_must_not_disable_gap_check(tmp_path, syntax):
    interval = "cycles = 360" if syntax == "cycles" else "seconds = 3600."
    run = _run(tmp_path, interval)
    if syntax == "double_quotes":
        nml = run / "m2_run.nml"
        nml.write_text(nml.read_text().replace(f"'{interval}'", f'"{interval}"'))
    # 360 internal cycles of 10 s, or 3600 seconds, both mean one hour.
    _history(run, "m2_0001.nc", [0, 24])
    verdict = check_run(run)
    assert not verdict["ok"], verdict


def test_t4_concurrent_refinements_reserve_output(tmp_path, monkeypatch):
    recipe = tmp_path / "recipe.yaml"
    land = tmp_path / "land.shp"
    recipe.touch()
    land.touch()
    out = tmp_path / "shared"
    barrier = Barrier(2)
    lock = Lock()
    writes = []

    def generator(cmd, *, env, cwd):
        # Hold the first child before its first write, while the second caller
        # executes the real wrapper's output guard. No meshing is performed.
        try:
            barrier.wait(timeout=2)
        except BrokenBarrierError:
            pass  # A fixed wrapper admits only one caller to this stub.
        with lock:
            target = Path(env["LR_OUT"])
            target.mkdir(parents=True, exist_ok=True)
            report = target / "report.json"
            writes.append(report.exists())
            report.write_text("stub product\n")
        return 0

    monkeypatch.setattr(refine_run.subprocess, "call", generator)
    args = [str(recipe), "--land", str(land), "--out", str(out)]
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: refine_run.main(args), range(2)))
    assert sorted(results) == [0, 2], (results, writes)


def test_t5_new_frozen_diagonal_is_not_an_inherited_violation():
    # A unit-scale square, original diagonal 1--3; all original edges have
    # r <= 0.2. Flip to diagonal 0--2 without moving any retained node/depth.
    depths = np.array([3., 4.5, 6.75, 4.5])
    base = np.array([[0, 1, 3], [1, 2, 3]])
    changed = np.array([[0, 1, 2], [0, 2, 3]])
    frozen = np.zeros(4, dtype=bool)
    _, before = finish_depths(depths, base, frozen, hmin=3, hmax=300, rfactor=0.2)
    assert before["max_r_frozen_pair"] <= 0.2
    _, after = finish_depths(depths, changed, frozen, hmin=3, hmax=300, rfactor=0.2)
    assert not after["converged_at_write_precision"], after
