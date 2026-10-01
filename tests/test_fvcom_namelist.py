"""FVCOM namelist values and moved cases (review of the extend tools, round 6)."""

from __future__ import annotations

import pytest

from fvcom_mesh_tools.io import fvcom_namelist
from fvcom_mesh_tools.io.fvcom_namelist import (
    check_fvcom_dirs,
    fortran_string,
    relocate_case,
    set_value,
)


def test_strings_are_quoted_for_fortran():
    assert fortran_string("/a/run'x/") == "'/a/run''x/'"
    with pytest.raises(ValueError):
        fortran_string("/a/\tb")
    with pytest.raises(ValueError):
        fortran_string("/a/漢")


def test_a_backslash_in_a_value_stays_literal():
    out = set_value(" INPUT_DIR = 'x',\n", "INPUT_DIR", fortran_string("/a\\tb/"))
    assert out == " INPUT_DIR = '/a\\tb/',\n"


def test_directories_longer_than_fvcom_keeps_are_refused(tmp_path, monkeypatch):
    # the test's own directory is the yardstick: one byte more is refused
    monkeypatch.setattr(fvcom_namelist, "FVCOM_DIR_MAX", len(f"{tmp_path}/"))
    check_fvcom_dirs(tmp_path)
    with pytest.raises(ValueError, match="bytes of a run directory"):
        check_fvcom_dirs(tmp_path / "x")


def _case(root):
    case = root / "refined"
    (case / "input").mkdir(parents=True)
    (case / "output").mkdir()
    (case / "output" / "m2_0001.nc").write_text("old")
    (case / "m2_run.nml").write_text(
        f" START_DATE = '2021-01-01 00:00:00',\n END_DATE = '2021-01-21 00:00:00',\n"
        f" INPUT_DIR = '{case}/input/',\n OUTPUT_DIR = '{case}/output/',\n")
    return case


def test_relocate_points_the_namelist_at_the_copy(tmp_path, monkeypatch):
    monkeypatch.setattr(fvcom_namelist, "FVCOM_DIR_MAX", 4096)
    case = _case(tmp_path)
    nml = relocate_case(case, tmp_path / "smoke" / "refined", end_date="2021-01-03 00:00:00")
    text = nml.read_text()
    dst = tmp_path / "smoke" / "refined"
    assert f"INPUT_DIR = '{dst}/input/'," in text and f"OUTPUT_DIR = '{dst}/output/'," in text
    assert "END_DATE = '2021-01-03 00:00:00'," in text
    assert list((dst / "output").iterdir()) == []


def test_relocate_refuses_a_moved_path_fvcom_would_cut(tmp_path, monkeypatch):
    """Round 6 F11: the original passed the guard, the moved one did not."""
    case = _case(tmp_path)
    monkeypatch.setattr(fvcom_namelist, "FVCOM_DIR_MAX", len(f"{case}/output/"))
    check_fvcom_dirs(case / "input", case / "output")
    with pytest.raises(ValueError, match="bytes of a run directory"):
        relocate_case(case, tmp_path / "smoke" / "refined")
    assert not (tmp_path / "smoke").exists()


def test_relocate_refuses_to_overlap_its_source(tmp_path, monkeypatch):
    """Review round 7 F6: relocating a case onto itself deleted it."""
    monkeypatch.setattr(fvcom_namelist, "FVCOM_DIR_MAX", 4096)
    case = _case(tmp_path)
    link = tmp_path / "link"
    link.symlink_to(case)
    for dst in (case, link, tmp_path, case / "smoke"):
        with pytest.raises(ValueError, match="overlaps"):
            relocate_case(case, dst)
    assert (case / "m2_run.nml").exists()


def test_relocate_keeps_the_previous_copy_when_the_move_fails(tmp_path, monkeypatch):
    """Review round 8 F3."""
    from pathlib import Path

    monkeypatch.setattr(fvcom_namelist, "FVCOM_DIR_MAX", 4096)
    case = _case(tmp_path)
    dst = tmp_path / "smoke" / "refined"
    relocate_case(case, dst)
    (dst / "marker").write_text("previous")
    real = Path.rename

    def flaky(self, target):
        if self.name == "case":
            raise OSError("injected")
        return real(self, target)

    monkeypatch.setattr(Path, "rename", flaky)
    with pytest.raises(OSError, match="injected"):
        relocate_case(case, dst)
    assert (dst / "marker").read_text() == "previous"
