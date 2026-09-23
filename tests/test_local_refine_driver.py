"""Tests for the seed state machine in ``notebooks/420_local_refine.py``.

The driver is a notebook, not an importable module: it loads a base mesh, a
land shapefile and GPL DistMesh at import time. So the pieces worth testing
are lifted out of its source and run against stubs. That is not elegant, and
the alternative was no coverage at all for the part of the generator that
decides which mesh is delivered -- where an adversarial review found three
defects that a production run would only show by luck (gpt-6-astra, second
review, 2026-09-22).
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14

DRIVER = Path(__file__).resolve().parents[1] / "notebooks" / "420_local_refine.py"


def driver_function(name: str, env: dict):
    """Compile one top-level function out of the driver into ``env``."""
    tree = ast.parse(DRIVER.read_text())
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[fn], type_ignores=[]), "<driver>", "exec"), env)
    return env[name]


def driver_block(start: str, stop: str, env: dict):
    """Run one top-level block of the driver, delimited by two source lines."""
    src = DRIVER.read_text().split("\n")
    i = next(k for k, ln in enumerate(src) if ln.startswith(start))
    j = next(k for k, ln in enumerate(src) if ln.startswith(stop))
    exec(compile("\n".join(src[i:j + 1]), "<driver>", "exec"), env)
    return env


def square_mesh():
    nodes = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]])
    elements = np.array([[0, 1, 2], [0, 2, 3]])
    return Fort14Mesh(title="base", nodes=nodes, depths=np.full(4, 8.0),
                      elements=elements, open_boundaries=[np.array([0, 1])],
                      land_boundaries=[(20, np.array([1, 2, 3, 0]))])


def serialise_env(tmp_path, seen):
    base = square_mesh()
    return {
        "np": np, "json": json, "Path": Path, "base": base, "sel": None,
        "rc": None, "recipe": Path("recipe.yaml"), "Fort14Mesh": Fort14Mesh,
        "om": SimpleNamespace(boundary_loops=lambda _t: [np.arange(4)]),
        "say": lambda *a: None, "write_fort14": write_fort14,
        "read_fort14": lambda p: (seen.append(p), read_fort14(p))[1],
        "verify_patch": lambda *a, **kw: {"ok": True},
        "boundary_after_patch": lambda *a, **kw: set(),
        "shapely": __import__("shapely"),
        "out14": tmp_path / "stale.14",
    }


def test_serialise_reads_the_path_it_was_given():
    """It wrote its argument and read a module-level global.

    Latent while every caller passed the same path, and exactly the kind of
    thing that surfaces the first time a per-attempt filename is used.
    """
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    seen: list[Path] = []
    env = serialise_env(tmp, seen)
    serialise = driver_function("serialise", env)
    base = env["base"]
    candidate = (base.nodes, base.elements, base.depths, np.arange(4))
    serialise(candidate, {"want_boundary": set()}, tmp / "candidate.14")
    assert seen == [tmp / "candidate.14"]


def seed_search_env(tmp_path, attempts, qa_failures, misses=None):
    """A stub search: ``attempts`` maps a seed to a candidate or an exception.

    ``misses`` maps a seed to the regions that seed left coarser than their
    target, which is a rejection of its own: a mesh can pass every gate and
    still not be the mesh that was asked for.
    """
    called: list[int] = []
    written: list[int] = []
    misses = misses or {}

    def attempt(seed):
        called.append(seed)
        outcome = attempts[seed]
        if isinstance(outcome, Exception):
            raise outcome
        return outcome, {"seed": seed}

    def serialise(candidate, out, path):
        written.append(int(candidate[0][0]))
        path.write_text(str(written[-1]))
        return SimpleNamespace(nodes=np.zeros((1, 2)), elements=np.zeros((1, 3)),
                               depths=np.zeros(1), n_nodes=1, n_elements=1), None

    def run_qa(written_mesh, **kw):
        n = qa_failures[written[-1]]
        return SimpleNamespace(
            n_gate_total=21, n_gate_failed=n,
            checks=[SimpleNamespace(check_id="c1", requirement=">= 30",
                                    observed="bad", status="fail")] * n)

    def achieved_per_region(written_mesh):
        names = list(misses.get(written[-1], []))
        return ({n: {"miss": "median 90.0 m is coarser than the 30 m target",
                     "target_h_m": 30.0, "n_edges": 4, "median_m": 90.0,
                     "p90_m": 95.0, "max_m": 99.0} for n in names}, names)

    return {
        "np": np, "json": json, "os": SimpleNamespace(environ={}),
        "achieved_per_region": achieved_per_region,
        "Path": Path, "OUT": tmp_path, "cfg": {"base_mesh": Path("b.14")},
        "recipe": SimpleNamespace(stem="r"), "reports": {},
        "say": lambda *a: None, "attempt": attempt, "serialise": serialise,
        "run_qa": run_qa, "read_fort14": lambda p: None,
        "sel": SimpleNamespace(retained=np.zeros((1, 3), dtype=int)),
        "introduced_violations": lambda checks, n, elements=None: [
            {"check": c.check_id, "kind": "element", "id": 0}
            for c in checks if c.status == "fail"],
        "_called": called, "_written": written,
        # The hires branch is one `if` at the top of the driver, so the search
        # block reads its flag.  False is the default path: everything the
        # recipe did before the branch existed.
        "_LADDER": False,
    }


def with_patch_violations(env: dict) -> dict:
    """Compile the driver's own attribution helper into a search environment.

    Not a stub: the rule it applies -- that the minimum-depth check is a
    report and not a gate on the hires branch -- is the thing under test, and
    a stub of it would test the test.
    """
    driver_function("patch_violations", env)
    driver_function("wall_pairs", env)
    return env


def run_search(env):
    with_patch_violations(env)
    return driver_block("# --------------------------------------------------------- the seed",
                        'say(f"accepted seed', env)


def test_a_seed_that_cannot_be_built_does_not_end_the_search():
    """A fill that loses a constrained point is one seed's problem.

    ``stitch_patch`` raises for a lost or merged fixed point, and with no
    boundary around ``attempt`` a bad first seed ended the whole run -- which
    is the opposite of what a search is for.
    """
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    env = seed_search_env(
        tmp,
        attempts={0: ValueError("constrained points are absent"),
                  1: (np.array([1]), None, None, np.arange(1))},
        qa_failures={1: 0})
    env["os"].environ = {"LR_SEEDS": "0,1"}
    run_search(env)
    assert env["_called"] == [0, 1]
    assert env["reports"]["seed"] == 1
    assert json.loads((tmp / "report.json").read_text())["attempts"][0]["error"]


def test_the_report_names_the_mesh_that_is_on_disk():
    """When no seed passes, the kept candidate must still be the written one.

    The old code re-wrote only when the accepted seed was not the first tried,
    so a failed run left seed 1's file beside seed 0's node map and QA.
    """
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    env = seed_search_env(
        tmp,
        attempts={0: (np.array([0]), None, None, np.arange(1)),
                  1: (np.array([1]), None, None, np.arange(1))},
        qa_failures={0: 1, 1: 2})
    env["os"].environ = {"LR_SEEDS": "0,1"}
    run_search(env)
    assert env["reports"]["seed"] == 0
    assert env["_written"][-1] == 0, "the last file written is not the kept one"


def test_a_seed_that_misses_the_target_resolution_is_not_accepted():
    """QA does not know what was asked for.

    Every gate can pass on a mesh that left a region at its base size -- a
    core too thin to hold an edge midpoint reported nothing at all and the
    run still succeeded (fourth review). A seed that misses is worse than one
    that passes, and if every seed misses the run must fail rather than write.
    """
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    env = seed_search_env(
        tmp,
        attempts={0: (np.array([0]), None, None, np.arange(1)),
                  1: (np.array([1]), None, None, np.arange(1))},
        qa_failures={0: 0, 1: 1},
        misses={0: ["futtsu_nori"]})
    env["os"].environ = {"LR_SEEDS": "0,1"}
    run_search(env)
    assert env["_called"] == [0, 1], "a resolution miss must not end the search"
    assert env["reports"]["seed"] == 1, (
        "seed 0 passed every gate but did not deliver the target; seed 1 did")

    tmp2 = Path(tempfile.mkdtemp())
    env = seed_search_env(
        tmp2,
        attempts={0: (np.array([0]), None, None, np.arange(1))},
        qa_failures={0: 0}, misses={0: ["futtsu_nori"]})
    env["os"].environ = {"LR_SEEDS": "0"}
    with pytest.raises(SystemExit, match="resolution"):
        run_search(env)


def test_no_candidate_at_all_still_leaves_a_report():
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    env = seed_search_env(tmp, attempts={0: ValueError("no"), 1: ValueError("no")},
                          qa_failures={})
    env["os"].environ = {"LR_SEEDS": "0,1"}
    with pytest.raises(SystemExit, match="no seed produced"):
        run_search(env)
    assert len(json.loads((tmp / "report.json").read_text())["attempts"]) == 2


def test_the_m2_step_is_one_the_output_interval_divides():
    """FVCOM aborts before its first step otherwise.

    `NC_OUT_INTERVAL must be an integer number of internal time steps` --
    every rank, on both cases. The CFL allowance was 1.7276 s, rounding to
    three significant figures gave 1.72, an internal step of 17.2 s, and
    1800 / 17.2 = 104.65. The earlier experiments used 1.5 s and passed by
    luck rather than by rule.
    """
    import runpy

    root = Path(__file__).resolve().parents[1]
    env = runpy.run_path(str(root / "notebooks/414_refine_m2_prep.py"))
    dividing_step = env["dividing_step"]
    isplit, interval = 10, 1800.0
    for allowance in (1.7276, 8.221, 5.0, 0.83, 3.1, 0.257):
        dte = dividing_step(interval, isplit, allowance)
        assert dte <= allowance, "a step may not exceed what the mesh allows"
        steps = interval / (dte * isplit)
        assert abs(steps - round(steps)) < 1e-9, (
            f"{dte} s leaves {steps} steps per output, which FVCOM refuses")
    assert dividing_step(interval, isplit, 1.7276) == 1.5


def test_the_minimum_depth_is_reported_not_gated_on_the_hires_branch(tmp_path):
    """The branch writes the depths as the source gives them, so a tidal flat
    fails run_qa's 2 m floor.  That is a tidal flat under wetting and drying,
    not a defect the seed can be blamed for; gating it would refuse the very
    field the option exists to deliver (owner, 2026-09-23).
    """
    env = seed_search_env(tmp_path, attempts={0: (np.array([0]), None, None, np.arange(1))},
                     qa_failures={0: 1})

    def run_qa(written_mesh, **kw):
        return SimpleNamespace(
            n_gate_total=21, n_gate_failed=1,
            checks=[SimpleNamespace(check_id="min_depth_clip",
                                    requirement=">= 2 m",
                                    observed="0.07 m", status="fail")])

    env["run_qa"] = run_qa
    env["os"].environ = {"LR_SEEDS": "0"}
    env["_LADDER"] = True
    run_search(env)
    attempt = env["reports"]["attempts"][0]
    assert attempt["qa"]["n_introduced"] == 0, (
        "a shallow node must not be counted against the seed on this branch")
    assert attempt["min_depth_reported_not_gated"] == 1, (
        "and the absence of the gate must not be the absence of the report")

    env2 = seed_search_env(tmp_path, attempts={0: (np.array([0]), None, None, np.arange(1))},
                      qa_failures={0: 1})
    env2["run_qa"] = run_qa
    env2["os"].environ = {"LR_SEEDS": "0"}
    run_search(env2)                      # _LADDER False: the default branch
    assert env2["reports"]["attempts"][0]["qa"]["n_introduced"] == 1, (
        "off the branch the floor is still a gate")


def test_wall_pairs_are_every_pair_a_split_made_coincident():
    env = {"np": np}
    wall_pairs = driver_function("wall_pairs", env)
    assert wall_pairs({}) is None
    co = np.array([0, 1, 2, 3, 1, 3, 3])   # node 1 has one copy, node 3 two
    assert sorted(wall_pairs({"copy_of": co})) == [(1, 4), (3, 5), (3, 6), (5, 6)]


def test_a_wall_that_closes_water_off_is_found_by_its_piece():
    """Two walls that cross into a small triangle strand the water inside it."""
    import importlib

    walls = importlib.import_module("fvcom_mesh_tools.walls")
    env = {"np": np}
    walls_cut_off = driver_function("walls_cut_off", env)
    n = 9
    g = np.arange(n) * 100.0
    gx, gy = np.meshgrid(g, g)
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            tri += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    tri = np.asarray(tri)
    node = lambda i, j: j * n + i                            # noqa: E731
    # a closed square of walls around the cell (3..5, 3..5), plus a pier
    ring = [node(3, 3), node(4, 3), node(5, 3), node(5, 4), node(5, 5),
            node(4, 5), node(3, 5), node(3, 4), node(3, 3)]
    pier = [node(7, 0), node(7, 1), node(7, 2)]
    we = np.vstack([walls.wall_edges_from_path(ring), walls.wall_edges_from_path(pier)])
    out, t2, copy_of, _ = walls.split_along_walls(xy, tri, we)
    obc = {node(0, j) for j in range(n)}
    bad = walls_cut_off(t2, len(out), copy_of, we, obc)
    assert sorted(bad) == list(range(8)), "the ring's eight edges, not the pier's"
    one = walls_cut_off(t2, len(out), copy_of, we, obc, xy=xy)
    assert len(one) == 1 and one[0] < 8, (
        "given coordinates, only the shortest touching edge is withdrawn")
