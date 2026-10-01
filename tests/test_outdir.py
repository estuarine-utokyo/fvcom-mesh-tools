"""Atomic output-directory reservation (outdir.reserve)."""

import pytest

from fvcom_mesh_tools.outdir import MARKER, reserve


def test_reserve_takes_a_new_or_empty_directory_once(tmp_path):
    d = reserve(tmp_path / "run")
    assert (d / MARKER).exists()
    with pytest.raises(SystemExit, match="reserved by another run"):
        reserve(d)                                   # a second taker loses
    (tmp_path / "empty").mkdir()
    assert reserve(tmp_path / "empty") == (tmp_path / "empty").resolve()


def test_reserve_refuses_a_directory_with_content(tmp_path):
    (tmp_path / "used").mkdir()
    (tmp_path / "used" / "old.nc").write_text("x")
    with pytest.raises(SystemExit, match="not empty"):
        reserve(tmp_path / "used")
    (tmp_path / "file").write_text("x")
    with pytest.raises(SystemExit, match="not a directory"):
        reserve(tmp_path / "file")
