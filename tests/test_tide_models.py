"""Ocean tide model constants (tide_models.py)."""

from datetime import datetime

import numpy as np
import pytest

from fvcom_mesh_tools.tide_models import (
    fvcom_spectral,
    load_nao,
    read_nao,
    sample_constants,
    spectral_text,
)


def _write_nao(path, amp_cm, phase_deg, xmin=139.0, ymax=35.0, d="1/   12"):
    """amp_cm, phase_deg: (ny, nx) with row 0 = NORTH; NaN = land."""
    ny, nx = amp_cm.shape
    head = [
        "Model name : TEST", "Content    : M2  amplitude and phase",
        "Unit       : 0.01  cm for amplitude, 0.01 deg for phase", "Date       : x",
        f"xmin = {xmin:7.2f}  xmax = {xmin + (nx - 1) / 12:7.2f}  ymin = "
        f"{ymax - (ny - 1) / 12:7.2f}  ymax = {ymax:7.2f}",
        f"dx   = {d}  dy   = {d}  mend = {nx:7d}  nend = {ny:7d}",
        "Default value = 999999  Format = (10i6)",
    ]
    body = []
    for j in range(ny):
        for field, scale in ((amp_cm, 100), (phase_deg, 100)):
            vals = [999999 if np.isnan(v) else int(round(v * scale)) for v in field[j]]
            body += ["".join(f"{v:6d}" for v in vals[k:k + 10]) for k in range(0, nx, 10)]
    path.write_text("\n".join(head + body) + "\n")


def test_read_nao_puts_row_zero_in_the_north_and_scales_units(tmp_path):
    amp = np.array([[10.0, 20.0, np.nan], [30.0, 40.0, 50.0]])      # north row first
    pha = np.array([[1.0, 2.0, np.nan], [3.0, 4.0, 5.0]])
    _write_nao(tmp_path / "m2_j.nao", amp, pha)
    g = read_nao(tmp_path / "m2_j.nao")
    assert g["constituent"] == "M2"
    assert np.allclose(g["lat"], [35.0 - 1 / 12, 35.0])             # ascending
    assert g["amp"][1, 0] == pytest.approx(0.10) and g["amp"][0, 0] == pytest.approx(0.30)
    assert np.isnan(g["amp"][1, 2]) and g["phase"][0, 2] == pytest.approx(5.0)


def test_read_nao_refuses_a_file_of_the_wrong_length(tmp_path):
    _write_nao(tmp_path / "x.nao", np.ones((2, 3)), np.ones((2, 3)))
    p = tmp_path / "x.nao"
    p.write_text("\n".join(p.read_text().splitlines()[:-1]) + "\n")
    with pytest.raises(ValueError, match="expected"):
        read_nao(p)


def test_sample_interpolates_through_the_phase_wrap_and_fills_at_the_coast(tmp_path):
    amp = np.full((3, 3), 50.0)
    pha = np.array([[350.0, 10.0, 10.0], [350.0, 10.0, 10.0], [350.0, 10.0, 10.0]])
    amp[0, 2] = np.nan
    _write_nao(tmp_path / "m2_j.nao", amp, pha)
    g = read_nao(tmp_path / "m2_j.nao")
    a, p, filled = sample_constants(g, np.array([139.0 + 0.5 / 12]), np.array([34.9]))
    # a vector mean of 350 and 10 deg: phase 0, amplitude cos(10 deg) of the input
    assert a[0] == pytest.approx(0.5 * np.cos(np.radians(10)), rel=1e-3)
    assert min(p[0], 360 - p[0]) < 1.0   # not 180
    # a point in the cell touching the land node takes the nearest ocean node
    a, p, filled = sample_constants(g, np.array([139.0 + 1.5 / 12]), np.array([35.0 - 0.5 / 12]))
    assert filled[0] and a[0] == pytest.approx(0.5)
    a, p, filled = sample_constants(g, np.array([150.0]), np.array([10.0]))
    assert np.isnan(a[0]) and not filled[0]


def test_load_nao_names_its_constituents(tmp_path):
    _write_nao(tmp_path / "m2_j.nao", np.ones((2, 2)), np.ones((2, 2)))
    assert set(load_nao(tmp_path, ["m2"])) == {"M2"}
    with pytest.raises(ValueError, match="no constituent"):
        load_nao(tmp_path, ["X9"])
    with pytest.raises(FileNotFoundError):
        load_nao(tmp_path, ["S2"])


def test_fvcom_spectral_round_trips_through_utide():
    # An FVCOM-style series built from the converted constants must analyse
    # back to the Greenwich constants it came from.
    import utide

    names = ["M2", "K1", "O1"]
    amp = np.array([[0.5], [0.25], [0.2]])
    phase = np.array([[150.0], [200.0], [300.0]])
    start = datetime(2021, 1, 1)
    period, a, ph = fvcom_spectral(names, amp, phase, start, datetime(2021, 1, 31), 35.0)
    t = np.arange(0, 60 * 86400, 1800.0)
    eta = sum(a[k, 0] * np.cos(2 * np.pi * t / period[k] - np.radians(ph[k, 0]))
              for k in range(3))
    days = start.toordinal() + t / 86400   # Python ordinal days, as utide's epoch="python"
    coef = utide.solve(days, eta, lat=35.0, constit=names, method="ols", conf_int="none",
                       verbose=False, epoch="python")
    got = dict(zip(coef.name, zip(coef.A, coef.g)))
    for k, c in enumerate(names):
        assert abs(got[c][0] - amp[k, 0]) < 0.003
        assert abs((got[c][1] - phase[k, 0] + 180) % 360 - 180) < 0.3


def test_spectral_text_rejects_bad_input():
    with pytest.raises(ValueError):
        spectral_text(["M2"], [44714.0], [[np.nan]], [[0.0]], "2021-01-01 00:00:00")
    with pytest.raises(ValueError):
        spectral_text(["M2", "K1"], [1.0, 2.0], [[1.0]], [[0.0]], "2021-01-01 00:00:00")
