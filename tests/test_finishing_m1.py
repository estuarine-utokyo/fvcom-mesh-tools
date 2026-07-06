"""Stage-2 M1: detector + planner tests."""
import numpy as np
import pytest

from fvcom_mesh_tools.finishing import detect_violations, plan_patches


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
