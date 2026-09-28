"""Coastline rules review, round 16."""

import subprocess

import numpy as np
from shapely.geometry import LineString, MultiPolygon, Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import _clear_of, island_rings
from fvcom_mesh_tools.walls import close_wall_pockets


def _const(v):
    f = lambda q: np.full(len(np.atleast_2d(q)), float(v))  # noqa: E731
    f.size_min = f.size_max = float(v)
    return f


def test_clearance_is_judged_on_the_ring_delivered():
    land = Polygon([(100, 100), (900, 100), (900, 900), (600, 900),
                    (600, 550), (400, 550), (400, 900), (100, 900)])
    water = box(0, 0, 1000, 1000).difference(box(450, 640, 550, 1100))

    def h(q):
        return np.clip(30 + 0.6 * (np.atleast_2d(q)[:, 1] - 100), 30, 300)

    rings, rep = island_rings(land, water, h, fine_h=60.000001)
    for r in rings:
        assert _clear_of(r, water.boundary, h, 0.5)[3]       # possible


def test_islands_are_judged_against_each_other():
    two = MultiPolygon([box(100, 100, 300, 300), box(100, 301, 300, 500)])
    rings, rep = island_rings(two, box(0, 0, 1000, 1000), _const(30.0))
    assert rep["n_islands_added"] == 1
    assert any("another ring" in s["why"] for s in rep["skipped"])


def test_pockets_a_metre_apart_do_not_both_become_islands():
    h = _const(200.0)
    walls = [LineString(box(0, 500, 400, 650).exterior.coords),
             LineString(box(0, 651, 400, 801).exterior.coords)]
    added, _, _ = close_wall_pockets(walls, box(-1000, -1000, 2000, 0), h,
                                     min_h=60, size_floor=200)
    rings, rep = island_rings(added, box(-500, 100, 1000, 1500), h, fine_h=60)
    assert rep["n_islands_added"] <= 1


def test_untracked_files_are_seen_whatever_git_is_told(tmp_path, monkeypatch):
    def git(*a):
        subprocess.run(["git", *a], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "t@example.com")
    git("config", "user.name", "t")
    (tmp_path / "entry.py").write_text("import helper\n")
    git("add", "entry.py")
    git("commit", "-qm", "entry")
    (tmp_path / "helper.py").write_text("x = 1\n")
    git("config", "status.showUntrackedFiles", "no")
    st = provenance.git_state(tmp_path / "entry.py")
    assert st["path_tracked"] and "helper.py" in st["dirty"]
