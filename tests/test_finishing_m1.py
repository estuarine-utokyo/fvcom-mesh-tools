"""Stage-2 M1: detector + planner tests."""
import numpy as np
import pytest

from fvcom_mesh_tools.autofinish import detect_violations, plan_patches


def _lattice(nx=8, ny=8):
    xs, ys = np.meshgrid(np.linspace(0, 1, nx), np.linspace(0, 1, ny))
    P = np.column_stack([xs.ravel(), ys.ravel()])
    T = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            T.append([a, a + 1, a + nx])
            T.append([a + 1, a + nx + 1, a + nx])
    return P, np.asarray(T)


def test_clean_lattice_no_violations():
    P, T = _lattice()
    det = detect_violations(P, T)
    assert all(not det[c]["elements"]
               for c in ("c1", "c2", "c4", "pinch"))
    assert plan_patches(P, T, det) == []


def test_c1_c4_detected_and_clustered():
    P, T = _lattice(14, 14)
    interior = 7 * 14 + 7         # central node, k-ring stays interior
    P = P.copy()
    P[interior] += [0.030, 0.028]  # skews incident elements
    det = detect_violations(P, T)
    assert det["c1"]["elements"] or det["c4"]["elements"]
    patches = plan_patches(P, T, det)
    assert len(patches) == 1      # one cluster
    p = patches[0]
    assert p["class"] in ("micro", "patch")
    assert p["n_elements"] > len(p["seed_elements"])  # k-ring grew


def test_bound_class_and_obc_lock():
    P, T = _lattice()
    P = P.copy()
    P[1] += [0.04, 0.045]         # boundary node -> boundary patch
    det = detect_violations(P, T)
    if not (det["c1"]["elements"] or det["c4"]["elements"]):
        pytest.skip("no violation created")
    patches = plan_patches(P, T, det)
    assert patches[0]["class"] == "bound"
    patches2 = plan_patches(P, T, det, obc_nodes=[0, 1, 2])
    assert patches2[0]["class"] == "obc-locked"


def test_pinch_detection():
    # two squares sharing ONE node -> pinch
    P1, T1 = _lattice(3, 3)
    P2 = P1 + np.array([1.0, 1.0]) - P1[8]  # weld at node 8
    P = np.vstack([P1, P2[1:]])
    remap = np.arange(9) + 8
    remap[0] = 8
    T2 = remap[T1]
    T = np.vstack([T1, T2])
    det = detect_violations(P, T)
    assert det["pinch"]["elements"]


def test_max_patches_cap():
    P, T = _lattice(30, 30)
    P = P.copy()
    for v in (10 * 30 + 10, 20 * 30 + 20):   # two far clusters
        P[v] += [0.012, 0.011]
    det = detect_violations(P, T)
    patches = plan_patches(P, T, det, max_patches=1)
    assert "note" in patches[-1]


def test_directive_refine_and_obc_protection():
    import types

    from pyproj import Transformer

    from fvcom_mesh_tools.autofinish import apply_directives

    P, T = _lattice(12, 12)
    P = P * 20000.0  # ~20 km square in fake UTM
    tr = Transformer.from_crs("EPSG:32654", "EPSG:4326",
                              always_xy=True)

    def to_lonlat(ring):
        x, y = tr.transform(ring[:, 0], ring[:, 1])
        return np.column_stack([x, y]).tolist()

    mesh = types.SimpleNamespace(
        nodes=P.copy(), elements=T.copy(),
        depths=np.full(len(P), 10.0),
        open_boundaries=[np.array([0, 1, 2, 3])],
    )
    inner = np.array([[6000., 6000.], [14000., 6000.],
                      [14000., 14000.], [6000., 14000.]])
    mesh, led = apply_directives(
        mesh, [{"polygon": to_lonlat(inner), "target_h_m": 700.0}],
        utm_epsg=32654)
    assert led[0]["outcome"].startswith("applied")
    assert len(mesh.elements) > len(T)          # refined
    assert np.isfinite(mesh.depths).all()
    assert len(mesh.open_boundaries[0]) == 4    # obc preserved
    # seam leftovers are healed by the auto stage (design order)
    from fvcom_mesh_tools.autofinish import execute_patches, plan_patches

    det = detect_violations(mesh.nodes, mesh.elements)
    patches = plan_patches(mesh.nodes, mesh.elements, det)
    mesh.nodes, _ = execute_patches(mesh.nodes, mesh.elements,
                                    patches, log=lambda *a: None)
    det2 = detect_violations(mesh.nodes, mesh.elements)
    n_after = sum(len(det2[c]["elements"]) for c in ("c1", "c2"))
    n_before = sum(len(det[c]["elements"]) for c in ("c1", "c2"))
    assert n_after <= n_before

    # directive over the OBC -> skipped
    over = np.array([[-1000., -1000.], [8000., -1000.],
                     [8000., 3000.], [-1000., 3000.]])
    mesh, led2 = apply_directives(
        mesh, [{"polygon": to_lonlat(over), "target_h_m": 500.0}],
        utm_epsg=32654)
    assert "skipped" in led2[0]["outcome"]


def test_collapse_short_boundary_edge():
    # A 137 m boundary step between ~400 m neighbours (G8-c5,
    # run 6191386): the triangle spanning it is a sliver no node
    # move can fix. The collapse merges the lower-valence
    # endpoint into the other and the local worst angle improves.
    from fvcom_mesh_tools.algorithms.obc_finish import (
        collapse_short_boundary_edges,
    )
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh

    nodes = np.array([
        [0.0, 0.0],        # A
        [400.0, 0.0],      # B  (short edge B-C = 137 m)
        [537.0, 0.0],      # C
        [900.0, 0.0],      # D
        [300.0, 380.0],    # E  (far apex -> sliver B,C,E)
        [750.0, 380.0],    # F
    ])
    els = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 5], [2, 5, 4]])
    mesh = Fort14Mesh(
        title="t", nodes=nodes, depths=np.full(6, 5.0),
        elements=els, open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 3]))])
    m2, info = collapse_short_boundary_edges(mesh)
    assert info["collapsed"] == 1
    assert len(m2.elements) == 3
    assert len(m2.nodes) == 5

    def worst(nds, e3):
        w = 180.0
        for t in e3:
            q = nds[t]
            for k in range(3):
                u = q[(k + 1) % 3] - q[k]
                v = q[(k + 2) % 3] - q[k]
                c = np.dot(u, v) / (np.linalg.norm(u)
                                    * np.linalg.norm(v))
                w = min(w, np.degrees(np.arccos(
                    np.clip(c, -1, 1))))
        return w

    assert worst(m2.nodes, m2.elements) > worst(nodes, els)
    # land-boundary chain keeps only surviving nodes, in order
    assert all(len(ids) >= 2 for _, ids in m2.land_boundaries)
