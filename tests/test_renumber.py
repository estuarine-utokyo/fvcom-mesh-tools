"""Renumbering changes the order and nothing else."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.renumber import (
    METHODS,
    element_permutation,
    hilbert_permutation,
    locality_stats,
    morton_permutation,
    rcm_permutation,
    renumber_mesh,
)


def grid_mesh(nx: int = 12, ny: int = 12, h: float = 100.0):
    gx, gy = np.meshgrid(np.arange(nx) * h, np.arange(ny) * h)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            tri += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    return nodes, np.asarray(tri, dtype=np.int64)


def shuffled_mesh(seed: int = 0):
    """A mesh whose numbering has been deliberately scrambled."""
    nodes, tri = grid_mesh()
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(nodes))
    inv = np.empty(len(perm), dtype=np.int64)
    inv[perm] = np.arange(len(perm))
    outer = np.concatenate([
        np.arange(12), np.arange(23, 144, 12), np.arange(142, 131, -1),
        np.arange(132, 11, -12)])
    return Fort14Mesh(
        title="scrambled", nodes=nodes[perm], depths=np.full(len(nodes), 8.0),
        elements=inv[tri], open_boundaries=[inv[np.array([0, 1])]],
        land_boundaries=[(20, inv[outer[1:]])])


@pytest.mark.parametrize("method", METHODS)
def test_renumbering_moves_the_order_and_nothing_else(method):
    """The mesh has to come out the same mesh, triangle for triangle."""
    mesh = shuffled_mesh()
    out, rep = renumber_mesh(mesh, method, only_if_better=False)
    perm = rep["node_perm"]
    assert sorted(perm.tolist()) == list(range(mesh.n_nodes)), (
        "a permutation has to be one")
    assert np.array_equal(out.nodes, mesh.nodes[perm])
    assert np.array_equal(out.depths, mesh.depths[perm])

    def canon(tri, back=None):
        t = tri if back is None else back[tri]
        return np.unique(np.sort(t, axis=1), axis=0)

    assert np.array_equal(canon(mesh.elements), canon(out.elements, perm)), (
        "the set of triangles must survive, as node coordinates")
    # the boundaries are node lists, so they must be relabelled in place
    assert np.array_equal(mesh.nodes[mesh.open_boundaries[0]],
                          out.nodes[out.open_boundaries[0]])
    assert np.array_equal(mesh.nodes[mesh.land_boundaries[0][1]],
                          out.nodes[out.land_boundaries[0][1]])


@pytest.mark.parametrize("method", METHODS)
def test_a_scrambled_mesh_gets_better(method):
    mesh = shuffled_mesh()
    _, rep = renumber_mesh(mesh, method, only_if_better=False)
    assert rep["after"]["mean_edge_span"] < rep["before"]["mean_edge_span"], (
        f"{method} did not improve a randomly numbered mesh")


def test_a_mesh_that_is_already_well_numbered_is_left_alone():
    """The production base has a bandwidth of 80 and RCM takes it to 128.

    A tool that always reports an improvement has not measured one. The
    guard also has to look past the mean: Morton on that mesh improves the
    mean edge span and takes the 99th percentile from 74 to 610, because the
    Z curve jumps the domain at every power-of-two boundary.
    """
    nodes, tri = grid_mesh()
    perm = rcm_permutation(len(nodes), tri)
    inv = np.empty(len(perm), dtype=np.int64)
    inv[perm] = np.arange(len(perm))
    nodes, tri = nodes[perm], inv[tri]
    mesh = Fort14Mesh(title="tidy", nodes=nodes, depths=np.full(len(nodes), 8.0),
                      elements=tri, open_boundaries=[], land_boundaries=[])
    for method in METHODS:
        out, rep = renumber_mesh(mesh, method)
        if not rep["applied"]:
            assert out is mesh
            assert np.array_equal(rep["node_perm"], np.arange(mesh.n_nodes))
            assert "already well numbered" in rep["note"]
        else:
            assert rep["mean_edge_span_ratio"] < 1.0
            assert rep["p99_edge_span_ratio"] < 1.0


def test_the_element_order_follows_the_nodes():
    mesh = shuffled_mesh()
    _, rep = renumber_mesh(mesh, "rcm", only_if_better=False)
    assert rep["after"]["mean_element_span"] < rep["before"]["mean_element_span"]
    ep = rep["element_perm"]
    assert sorted(ep.tolist()) == list(range(mesh.n_elements))


def test_the_orderings_are_permutations_of_the_right_length():
    nodes, tri = grid_mesh()
    for perm in (rcm_permutation(len(nodes), tri), morton_permutation(nodes),
                 hilbert_permutation(nodes)):
        assert sorted(perm.tolist()) == list(range(len(nodes)))
    ep = element_permutation(tri, None, nodes, "centroid")
    assert sorted(ep.tolist()) == list(range(len(tri)))


def test_rcm_beats_both_space_filling_curves_on_these_meshes():
    """Measured, against the expectation rather than for it.

    The textbook argument is that Hilbert beats Morton because the Z curve
    jumps the width of the domain at every power-of-two boundary. On these
    meshes it does not: on a 33x33 grid Morton's mean edge span is 22.4 and
    Hilbert's 25.3, and Morton's p99 is the better of the two as well. What
    is true on every grid measured is that RCM beats both on every metric --
    bandwidth 33 against 545 and 879 -- because it optimises the quantity
    being measured and a space-filling curve optimises a proxy for it.
    """
    for n in (17, 33, 65):
        nodes, tri = grid_mesh(n, n)
        stats = {}
        for name, perm in (("morton", morton_permutation(nodes)),
                           ("hilbert", hilbert_permutation(nodes)),
                           ("rcm", rcm_permutation(len(nodes), tri))):
            inv = np.empty(len(perm), dtype=np.int64)
            inv[perm] = np.arange(len(perm))
            stats[name] = locality_stats(len(nodes), inv[tri])
        assert stats["rcm"]["bandwidth"] <= n, (
            "RCM on a structured grid should reach the grid's own width")
        for other in ("morton", "hilbert"):
            assert stats["rcm"]["bandwidth"] < stats[other]["bandwidth"]
            assert stats["rcm"]["mean_edge_span"] < stats[other]["mean_edge_span"]
            assert stats["rcm"]["p99_edge_span"] < stats[other]["p99_edge_span"]


def test_the_written_case_carries_every_file_and_the_permutation(tmp_path: Path):
    """A renumbered case with no record of how is one nobody can check."""
    from fvcom_mesh_tools.cli.renumber import main

    mesh = shuffled_mesh()
    src = tmp_path / "src"
    written = export_fvcom_case(mesh, src, "c", cor=mesh.nodes[:, 1] * 0 + 35.0,
                                sponge=[(int(mesh.open_boundaries[0][0]),
                                         1000.0, 0.001)],
                                obc_depth_control=False)
    assert "spg" in written and "cor" in written
    out = tmp_path / "out"
    assert main([str(src / "c"), "--outdir", str(out), "--casename", "r"]) == 0

    back = read_fvcom_case(out / "r_grd.dat", out / "r_dep.dat", out / "r_obc.dat")
    rec = json.loads((out / "r_renumber.json").read_text())
    perm = np.array(rec["node_perm_old_ids"], dtype=np.int64)
    assert np.array_equal(back.nodes, mesh.nodes[perm])
    assert np.array_equal(back.depths, mesh.depths[perm])
    # the sponge is keyed by node id and has to have moved with it
    spg = [ln.split() for ln in (out / "r_spg.dat").read_text().splitlines()[1:]
           if ln.strip()]
    inv = np.empty(len(perm), dtype=np.int64)
    inv[perm] = np.arange(len(perm))
    assert int(spg[0][0]) - 1 == int(inv[int(mesh.open_boundaries[0][0])])
    cor = [ln.split() for ln in (out / "r_cor.dat").read_text().splitlines()[1:]
           if ln.strip()]
    assert len(cor) == mesh.n_nodes


def test_a_case_that_would_not_improve_is_refused_rather_than_written(tmp_path):
    """Renumbering an already renumbered case buys nothing, and says so.

    A lexicographic grid is NOT already well numbered -- RCM takes its mean
    edge span from 8.5 to 5.9 -- so the mesh here is one that has been
    through RCM once. The production base behaves the same way for the same
    reason: somebody numbered it well already.
    """
    from fvcom_mesh_tools.cli.renumber import main

    nodes, tri = grid_mesh()
    perm = rcm_permutation(len(nodes), tri)
    inv = np.empty(len(perm), dtype=np.int64)
    inv[perm] = np.arange(len(perm))
    nodes, tri = nodes[perm], inv[tri]
    mesh = Fort14Mesh(title="tidy", nodes=nodes, depths=np.full(len(nodes), 8.0),
                      elements=tri, open_boundaries=[np.array([0, 1])],
                      land_boundaries=[])
    src = tmp_path / "src"
    export_fvcom_case(mesh, src, "c", obc_depth_control=False)
    out = tmp_path / "out"
    assert main([str(src / "c"), "--outdir", str(out)]) == 1
    assert not (out / "c_grd.dat").exists()
    assert main([str(src / "c"), "--outdir", str(out), "--force"]) == 0
    assert (out / "c_grd.dat").exists()


def test_the_sweep_starts_where_it_is_told_and_sms_is_the_default():
    """SMS selects the open-boundary nodestring and renumbers from it.

    That is not a detail of taste: measured on the production base, seeding
    the sweep at the open boundary gives a bandwidth of 78 against the
    mesh's own 80, and letting scipy choose its own start gives 128. It also
    leaves the open boundary as a contiguous block at the end of the
    numbering, which is what the convention is FOR.
    """
    from fvcom_mesh_tools.renumber import open_boundary_seeds

    mesh = shuffled_mesh()
    obc = open_boundary_seeds(mesh)
    assert obc is not None and obc.size == 2

    out, rep = renumber_mesh(mesh, "rcm", seeds="obc", only_if_better=False)
    assert rep["seeds"] == "obc"
    new_obc = np.asarray(out.open_boundaries[0])
    assert new_obc.max() == out.n_nodes - 1, (
        "an OBC-seeded RCM sweep is reversed, so the open boundary lands last")
    assert new_obc.max() - new_obc.min() == new_obc.size - 1, (
        "and lands as one contiguous block")

    auto, _ = renumber_mesh(mesh, "rcm", seeds="auto", only_if_better=False)
    assert not np.array_equal(auto.nodes, out.nodes), (
        "the automatic start is a different ordering, not the same one")


def test_seeding_the_whole_boundary_is_offered_and_is_worse():
    """A selection marks where to start, not what to number first.

    A front as wide as the coastline produces a bandwidth as wide as the
    coastline: on the production base, 781 against the mesh's own 80.
    """
    mesh = shuffled_mesh()
    wide, _ = renumber_mesh(mesh, "rcm", seeds="boundary", only_if_better=False)
    narrow, rep = renumber_mesh(mesh, "rcm", seeds="obc", only_if_better=False)
    from fvcom_mesh_tools.renumber import locality_stats
    a = locality_stats(mesh.n_nodes, wide.elements)
    b = locality_stats(mesh.n_nodes, narrow.elements)
    assert b["bandwidth"] <= a["bandwidth"]


def test_a_seed_that_is_not_a_node_is_refused():
    mesh = shuffled_mesh()
    with pytest.raises(ValueError, match="not a node"):
        renumber_mesh(mesh, "rcm", seeds=[mesh.n_nodes + 5], only_if_better=False)
