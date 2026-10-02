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


def test_a_stage_writes_only_into_its_own_runs_reservation(tmp_path, monkeypatch):
    """Review of the extend tools, round 19 F2."""
    from fvcom_mesh_tools.outdir import TOKEN_ENV, claim, new_token, reserve

    token = new_token()
    out = reserve(tmp_path / "build", token=token)
    monkeypatch.setenv(TOKEN_ENV, token)
    assert claim(out) == out
    assert claim(out / "generate", owner=out) == out / "generate"
    monkeypatch.setenv(TOKEN_ENV, new_token())
    with pytest.raises(SystemExit, match="another run"):
        claim(out)
    monkeypatch.delenv(TOKEN_ENV)
    with pytest.raises(SystemExit):          # alone: an earlier output is refused
        claim(out)
    assert claim(tmp_path / "fresh") == (tmp_path / "fresh").resolve()
