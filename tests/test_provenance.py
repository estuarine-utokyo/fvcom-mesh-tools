from __future__ import annotations

import hashlib
import subprocess

import pytest

from fvcom_mesh_tools.provenance import collect, file_sha256, git_state


def test_file_sha256_matches_hashlib_and_tolerates_a_missing_file(tmp_path):
    f = tmp_path / "a.txt"
    f.write_bytes(b"coast")
    assert file_sha256(f) == hashlib.sha256(b"coast").hexdigest()
    assert file_sha256(tmp_path / "missing") is None


def _git(root, *args):
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    try:
        _git(tmp_path, "init", "-q")
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("git is not available")
    _git(tmp_path, "config", "user.email", "t@example.com")
    _git(tmp_path, "config", "user.name", "t")
    (tmp_path / "code.py").write_text("x = 1\n")
    _git(tmp_path, "add", "code.py")
    _git(tmp_path, "commit", "-q", "-m", "c")
    return tmp_path


def test_git_state_reports_the_commit_and_uncommitted_changes(repo):
    clean = git_state(repo / "code.py")
    assert len(clean["commit"]) == 40 and clean["dirty"] == []
    (repo / "code.py").write_text("x = 2\n")
    assert git_state(repo)["dirty"] == ["code.py"]


def test_git_state_is_none_outside_a_work_tree(tmp_path):
    assert git_state(tmp_path) is None


def test_collect_hashes_a_shapefile_set_together_and_records_libraries(repo, tmp_path):
    parts = []
    for ext in (".shp", ".dbf"):
        p = tmp_path / f"land{ext}"
        p.write_bytes(ext.encode())
        parts.append(p)
    rec = collect(code={"tool": repo}, files={"land": parts, "gone": tmp_path / "x"},
                  libraries=("numpy", "no_such_library_xyz"))
    want = hashlib.sha256("".join(file_sha256(p) for p in parts).encode()).hexdigest()
    assert rec["files"]["land"]["sha256"] == want
    assert rec["files"]["gone"]["sha256"] is None
    assert rec["libraries"]["numpy"] and rec["libraries"]["no_such_library_xyz"] is None
    assert rec["code"]["tool"]["commit"]
    # the set's hash changes when any member does
    parts[1].write_bytes(b"other")
    assert collect(files={"land": parts})["files"]["land"]["sha256"] != want


def test_collect_imports_no_library(monkeypatch):
    """review round 1: importing oceanmesh (GPL) from this package is barred."""
    import importlib

    def refuse(name, *a, **k):
        raise AssertionError(f"imported {name}")

    monkeypatch.setattr(importlib, "import_module", refuse)
    rec = collect(libraries=("numpy", "no_such_library_xyz"))
    assert rec["libraries"]["numpy"] and rec["libraries"]["no_such_library_xyz"] is None


def test_git_state_counts_untracked_files_and_whether_the_path_is_tracked(repo):
    new = repo / "driver.py"
    new.write_text("print(1)\n")
    st = git_state(new)
    assert "driver.py" in st["dirty"] and st["path_tracked"] is False
    assert git_state(repo / "code.py")["path_tracked"] is True


def test_dataset_files_takes_a_shapefiles_sidecars_and_a_geojson_alone(tmp_path):
    from fvcom_mesh_tools.provenance import dataset_files

    for ext in (".shp", ".shx", ".dbf", ".prj"):
        (tmp_path / f"land{ext}").write_text(ext)
    got = {p.suffix for p in dataset_files(tmp_path / "land.shp")}
    assert got == {".shp", ".shx", ".dbf", ".prj"}
    gj = tmp_path / "area.geojson"
    gj.write_text("{}")
    assert dataset_files(gj) == [gj]
