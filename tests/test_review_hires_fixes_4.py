"""The fourth review's findings (docs/hires_walls_tools_review_4.md) as
regression tests: written as failing probes by the reviewer (gpt-6-astra),
passing since the fixes."""

import json
import os
import time
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from fvcom_mesh_tools.cli.check_run import check_run
from fvcom_mesh_tools.cli.refine_depths import finish_depths, inherited_edges_of
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14

ROOT = Path(__file__).resolve().parents[1]


def make_run(root, hours, *, first=0, end=24, prefix="", case="m2"):
    (root / "output").mkdir(parents=True)
    (root / "fvcom.log").write_text("TADA!\n")
    origin = datetime(2021, 1, 1)
    def stamp(h):
        return (origin + timedelta(hours=h)).strftime("%Y-%m-%d %H:%M:%S")

    (root / f"{case}_run.nml").write_text(
        prefix + "&NML_CASE\n START_DATE = '2021-01-01 00:00:00',\n"
        f" END_DATE = '{stamp(end)}',\n/\n&NML_NETCDF\n"
        f" NC_FIRST_OUT = '{stamp(first)}',\n"
        " NC_OUT_INTERVAL = 'seconds = 3600.',\n/\n"
    )
    with netCDF4.Dataset(root / "output" / f"{case}_0001.nc", "w") as ds:
        ds.createDimension("time", len(hours))
        ds.createDimension("DateStrLen", 26)
        ds.createDimension("node", 4)
        ds.createDimension("nele", 2)
        times = ds.createVariable("Times", "S1", ("time", "DateStrLen"))
        for i, hour in enumerate(hours):
            value = (origin + timedelta(hours=hour)).strftime("%Y-%m-%dT%H:%M:%S.%f")
            times[i] = np.asarray(list(value), dtype="S1")
        for name, dim in (("zeta", "node"), ("ua", "nele"), ("va", "nele")):
            ds.createVariable(name, "f4", ("time", dim))[:] = 0.1
    return root


def test_u1_commented_interval_cannot_hide_missing_records(tmp_path):
    run = make_run(tmp_path, [0, 24], prefix="! NC_OUT_INTERVAL = 'days = 1',\n")
    verdict = check_run(run)
    assert not verdict["ok"], verdict


def test_u2_valid_single_record_history_returns_a_verdict(tmp_path):
    # Output starts at END_DATE: a single scheduled record is complete.
    run = make_run(tmp_path, [24], first=24)
    try:
        verdict = check_run(run)
    except ValueError as exc:
        pytest.fail(f"one valid record crashes the checker: {exc}")
    assert verdict["ok"], verdict


def square_mesh(elements):
    return Fort14Mesh(
        title="tiny provenance witness",
        nodes=np.array([[390000., 3900000.], [390100., 3900000.],
                        [390100., 3900100.], [390000., 3900100.]]),
        depths=np.array([3., 4.5, 6.75, 4.5]),
        elements=np.asarray(elements), open_boundaries=[], land_boundaries=[],
    )


def test_u3_replaced_base_cannot_excuse_a_new_frozen_edge(tmp_path):
    base = tmp_path / "base.14"
    original = square_mesh([[0, 1, 3], [1, 2, 3]])
    changed = square_mesh([[0, 1, 2], [0, 2, 3]])
    write_fort14(original, base)
    (tmp_path / "report.json").write_text(json.dumps({"base_mesh": str(base)}))
    node_map = np.arange(4)
    inherited = inherited_edges_of(tmp_path, node_map)
    _, before = finish_depths(changed.depths, changed.elements, np.zeros(4, bool),
                             hmin=3, hmax=300, rfactor=0.2, inherited_edges=inherited)
    assert not before["converged_at_write_precision"]
    # Only the external base file changes. The refinement, mapping and report
    # still refer to the earlier base; a path alone is not its identity.
    write_fort14(changed, base)
    inherited = inherited_edges_of(tmp_path, node_map)
    _, after = finish_depths(changed.depths, changed.elements, np.zeros(4, bool),
                            hmin=3, hmax=300, rfactor=0.2, inherited_edges=inherited)
    assert not after["converged_at_write_precision"], after


def reservation_block(out):
    """Execute only the actual driver's admission block, before any meshing."""
    source = (ROOT / "notebooks/420_local_refine.py").read_text()
    block = source[source.index('_res = OUT / ".reserved"'):source.index('LAND = Path(')]
    exec(compile(block, "420_local_refine.py:reservation", "exec"),
         {"OUT": out, "os": os, "time": time})


def test_u4_stale_delegated_token_must_not_reclaim_nonempty_output(tmp_path, monkeypatch):
    (tmp_path / "report.json").write_text('{"previous_attempt": true}\n')
    monkeypatch.setenv("LR_RESERVATION", "missing-previous-reservation")
    with pytest.raises(SystemExit):
        reservation_block(tmp_path)


def test_u3_the_recorded_base_is_honoured_when_it_is_unchanged(tmp_path):
    import hashlib

    base = tmp_path / "base.14"
    original = square_mesh([[0, 1, 3], [1, 2, 3]])
    write_fort14(original, base)
    (tmp_path / "report.json").write_text(json.dumps({
        "base_mesh": str(base), "base_mesh_sha256": hashlib.sha256(base.read_bytes()).hexdigest()}))
    inherited = inherited_edges_of(tmp_path, np.arange(4))
    assert inherited is not None and len(inherited) == 5


def test_u4_a_fresh_direct_run_reserves_an_empty_directory(tmp_path, monkeypatch):
    monkeypatch.delenv("LR_RESERVATION", raising=False)
    out = tmp_path / "fresh"
    reservation_block(out)
    assert (out / ".reserved").exists()
    with pytest.raises(SystemExit):
        reservation_block(out)                   # taken now
