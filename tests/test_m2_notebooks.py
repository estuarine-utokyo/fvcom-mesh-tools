"""Small analytical checks for the numbered M2 experiment scripts."""

import importlib.util
from pathlib import Path

import netCDF4
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "notebooks" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prep = load("383_m2_case_prep")
analysis = load("384_m2_analysis")


def test_harmonics_overtides_and_wrapped_phase():
    period = 44714.16
    t = np.arange(3 * 86400, 11 * 86400 + 1, 1800.0)
    w = 2 * np.pi * t / period
    z = 0.23 + 0.8 * np.cos(w - np.deg2rad(359)) + 0.1 * np.cos(2 * w + 0.5)
    z += 0.03 * np.sin(3 * w)
    amp, phase, coef = analysis.harmonic_fit(t, z, period)
    assert amp == pytest.approx(0.8, abs=1e-12)
    assert phase == pytest.approx(359, abs=1e-10)
    assert coef[0] == pytest.approx(0.23, abs=1e-12)
    assert analysis.phase_difference(1, 359) == pytest.approx(2)
    with pytest.raises(ValueError, match="Nonfinite"):
        analysis.harmonic_fit(t, np.full_like(t, np.nan), period)
    with pytest.raises(ValueError, match="Insufficient"):
        analysis.harmonic_fit(t[:3], z[:3], period)


def test_detided_trend():
    t = np.arange(0, 8 * 86400, 1800.0)
    z = 1e8 + 2e5 * np.cos(2 * np.pi * t / 44714) + 2500 * t / 86400
    _, _, coef = analysis.harmonic_fit(t, z, 44714, trend=True)
    assert coef[-1] == pytest.approx(2500, abs=1e-6)


def test_arc_position_and_ascii_spectral_ordinals():
    arc = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0]])
    position, offset = prep.arc_position(np.array([[5.0, 2.0], [12.0, 5.0]]), arc)
    np.testing.assert_allclose(position, [5, 15])
    np.testing.assert_allclose(offset, [2, 2])
    text = prep.spectral_text(44714, [0.3, 0.4], [359, 1])
    lines = text.splitlines()
    assert lines[:2] == ["Tidal Component Number = 1", "1 = M2 44714.0000000000"]
    for key in ("Amplitude", "Phase", "Eref"):
        start = lines.index(key)
        assert [line.split()[0] for line in lines[start + 1 : start + 3]] == ["1", "2"]
        assert lines[start + 3] == key
    with pytest.raises(ValueError, match="Repeated"):
        prep.arc_position(np.array([[1.0, 1.0]]), np.zeros((2, 2)))


def fixture_run(tmp_path, corrupt=False, incomplete=False):
    xy = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
    tri = np.array([[0, 1, 2]])
    np.savez(tmp_path / "mesh.npz", xy=xy, tri=tri, depth=np.full(3, 10.0), obc=[0, 1])
    (tmp_path / "output").mkdir()
    t = np.arange(0, 11 * 86400 + 1, 1800.0)
    if incomplete:
        t = t[:-1]
    with netCDF4.Dataset(tmp_path / "output/m2_0001.nc", "w") as ds:
        for key, n in [("time", len(t)), ("node", 3), ("nele", 1), ("three", 3), ("siglay", 1)]:
            ds.createDimension(key, n)
        v = ds.createVariable("time", "f8", ("time",))
        v.units = "seconds since 2021-01-01 00:00:00"
        v[:] = t
        ds.createVariable("nv", "i4", ("three", "nele"))[:] = [[1], [3], [2]]
        ds.createVariable("h", "f8", ("node",))[:] = 10
        z = ds.createVariable("zeta", "f8", ("time", "node"))
        z[:] = np.tile(0.5 * np.cos(2 * np.pi * t / 44714)[:, None], (1, 3))
        if corrupt:
            z[4, 0] = np.nan
        for key in ("u", "v"):
            ds.createVariable(key, "f8", ("time", "siglay", "nele"))[:] = 0.1
    return dict(
        epoch_utc="2021-01-01 00:00:00",
        spinup_seconds=3 * 86400,
        duration_seconds=11 * 86400,
        period_seconds=44714,
    )


def test_output_health_and_volume(tmp_path):
    meta = fixture_run(tmp_path)
    *_, health = analysis.read_run(tmp_path, meta)
    assert health["complete_pass"] and health["finite_pass"]
    assert health["initial_volume_m3"] == pytest.approx(21)
    assert health["detided_analysis_drift_m3_per_day"] == pytest.approx(0, abs=1e-10)
    assert health["max_horizontal_speed_m_s"] == pytest.approx(np.sqrt(0.02))


def test_output_rejects_truncation(tmp_path):
    meta = fixture_run(tmp_path, incomplete=True)
    with pytest.raises(ValueError, match="incomplete run"):
        analysis.read_run(tmp_path, meta)


def test_output_records_nan(tmp_path):
    meta = fixture_run(tmp_path, corrupt=True)
    *_, health = analysis.read_run(tmp_path, meta)
    assert not health["finite_pass"]
    assert health["nonfinite_or_masked"]["zeta"] == 1
