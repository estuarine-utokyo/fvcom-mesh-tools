"""Coastline rules review, round 21."""

import subprocess

from shapely.affinity import scale
from shapely.geometry import Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings

A = box(100, 100, 900, 900).difference(box(100, 100, 110, 500))
B = box(-100, 200, 94, 450)


def test_a_lake_kept_by_two_simultaneous_flips_survives_either_order():
    c = scale(A, xfact=-1, yfact=1, origin=(-3, 0))
    shell = box(-1100, -100, 1100, 1100)
    water = box(-2000, -1000, 2000, 2000)
    for lakes in ([A, c, B], [B, A, c]):
        land = Polygon(shell.exterior.coords, [lk.exterior.coords for lk in lakes])
        rings, rep = island_rings(land, water, 30.0)
        assert rep["n_lakes_added"] == 3, rep["skipped"]
        assert _with_rings(water, rings).intersection(B).area > 48000.0


def test_an_invalid_first_group_does_not_end_the_search():
    shell = box(-300, -100, 1100, 1100)
    water = box(-1000, -1000, 2000, 2000).difference(box(-50, 280, 0, 330))
    land = Polygon(shell.exterior.coords, [A.exterior.coords, B.exterior.coords])
    rings, rep = island_rings(land, water, 30.0)
    assert rep["n_islands_added"] == 1 and rep["n_lakes_added"] == 2


def test_an_untracked_path_in_a_checkout_is_hashed(tmp_path):
    def git(*a):
        subprocess.run(["git", *a], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "t@example.com")
    git("config", "user.name", "t")
    (tmp_path / ".gitignore").write_text("build/\n")
    git("add", ".gitignore")
    git("commit", "-qm", "ignore")
    pkg = tmp_path / "build" / "pkg"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("x = 1\n")
    first = provenance.code_state("nope", pkg / "__init__.py")
    (pkg / "__init__.py").write_text("x = 2\n")
    second = provenance.code_state("nope", pkg / "__init__.py")
    assert first["git"] and not first["commit_identifies_code"]
    assert first["source_sha256"] and first["source_sha256"] != second["source_sha256"]
