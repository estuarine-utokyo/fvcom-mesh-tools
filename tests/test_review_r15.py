"""Coastline rules review, round 15."""

import os

import numpy as np
from shapely.geometry import Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import island_rings


def test_a_resampled_island_stays_inside_the_water():
    water = box(0, 0, 1000, 1000).difference(box(450, 600, 550, 1100))
    land = Polygon([(100, 100), (900, 100), (900, 900), (600, 900),
                    (600, 550), (400, 550), (400, 900), (100, 900)])

    def h(q):
        return np.clip(30 + 0.6 * (np.atleast_2d(q)[:, 1] - 100), 30, 300)

    rings, rep = island_rings(land, water, h, fine_h=60.000001)
    for r in rings:
        assert Polygon(r).within(water)


def test_clearance_is_bounded_between_samples():
    def h(q):
        y = np.atleast_2d(q)[:, 1]
        return 30 + np.maximum(0, 1 - np.abs(y - 502))

    rings, rep = island_rings(box(15.2, 100, 900, 900), box(0, 0, 1000, 1000), h)
    assert rep["n_islands_added"] == 1 and rep["tight"]


def test_a_link_back_to_an_ancestor_is_part_of_the_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)
    digests = []
    for k, link in enumerate((False, True)):
        pkg = tmp_path / str(k) / "pkg"
        (pkg / "sub").mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        (pkg / "mod.py").write_text("x = 1\n")
        (pkg / "sub" / "__init__.py").write_text("")
        if link:
            os.symlink(pkg, pkg / "sub" / "up")
        digests.append(provenance.code_state("nope", pkg)["source_sha256"])
    assert None not in digests and digests[0] != digests[1]
