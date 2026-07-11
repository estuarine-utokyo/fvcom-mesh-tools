# Finishing chain on the UTM sample-repro mesh (59e-proven tools):
#   1. align_open_boundary_local -- per-violating-node OBC
#      perpendicularity, quality-gated (no C1/C2/C4/C5 regression)
#   2. phase_h_finish -- stochastic local fixer + targeted
#      vertex_remove for C1/C4 residuals (seeded, reproducible,
#      OBC nodes frozen)
#   3. compact_nodes -- vertex_remove can orphan nodes
#   4. diagonal flips at still-violating OBC nodes: a wedge split
#      2->3 creates a near-perpendicular interior edge (this is
#      exactly the sample's SE-corner structure: corner angles
#      48/34/41 deg with the middle edge 7.7 deg from perp); a
#      pure node move cannot reach it. Quality-gated, topology-only.
#   5. verify the 13 constrained OBC nodes did not move; STOP if
#      they did (the input arc is sacrosanct)
import numpy as np
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.algorithms.perp_local import align_open_boundary_local
from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish

SRC = "outputs/sample_repro/sample_repro_utm.14"
DST = "outputs/sample_repro/sample_repro_final.14"


def _tri_angles(p0, p1, p2):
    out = []
    for x, y, z in ((p0, p1, p2), (p1, p2, p0), (p2, p0, p1)):
        u, v = y - x, z - x
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return out


def _area(p0, p1, p2):
    return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1])
                  - (p1[1] - p0[1]) * (p2[0] - p0[0]))


def flip_for_obc_perp(mesh, dev_max=20.0, min_angle=30.0,
                      max_angle=130.0, max_area_change=0.5):
    """IN-PLACE diagonal flips at OBC nodes whose best incident
    edge is still > dev_max from perpendicular: flip the opposite
    edge of a 1-ring element so its far vertex connects to the OBC
    node. Accept only if the new edge lands within dev_max and both
    new triangles pass C1/C2/positive-area and their mutual area
    ratio passes C4. Node positions never change."""
    nodes, els = mesh.nodes, mesh.elements
    seg = np.asarray(mesh.open_boundaries[0], int)
    tang = boundary_tangents(nodes[seg])
    fixed, unfixed = [], []
    for k, v in enumerate(seg):
        v = int(v)
        t_v = tang[k]

        def _dev_to(w):
            e = nodes[w] - nodes[v]
            e = e / np.linalg.norm(e)
            return abs(90.0 - np.degrees(
                np.arccos(min(1.0, abs(float(np.dot(e, t_v)))))))

        ring = np.where((els == v).any(axis=1))[0]
        devs = [_dev_to(int(w)) for ei in ring for w in els[ei]
                if int(w) != v]
        if not devs or min(devs) <= dev_max:
            continue
        best = None
        for ei in ring:
            a, b = [int(x) for x in els[ei] if int(x) != v]
            nb = [int(ej) for ej in np.where(
                ((els == a).any(axis=1)) & ((els == b).any(axis=1))
            )[0] if int(ej) != int(ei)]
            if not nb:
                continue
            ej = nb[0]
            m = int([x for x in els[ej] if int(x) not in (a, b)][0])
            if m == v:
                continue
            dev = _dev_to(m)
            if dev > dev_max - 1.0:
                continue
            t1, t2 = [v, a, m], [v, m, b]
            if _area(*nodes[t1]) < 0:
                t1 = [v, m, a]
            if _area(*nodes[t2]) < 0:
                t2 = [v, b, m]
            A1, A2 = _area(*nodes[t1]), _area(*nodes[t2])
            if A1 <= 0 or A2 <= 0:
                continue
            ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
            if min(ang) < min_angle or max(ang) > max_angle:
                continue
            if abs(A1 - A2) / max(A1, A2) > max_area_change:
                continue
            if best is None or dev < best[0]:
                best = (dev, int(ei), int(ej), t1, t2)
        if best is None:
            unfixed.append(v)
            continue
        _, ei, ej, t1, t2 = best
        els[ei] = t1
        els[ej] = t2
        fixed.append((v, round(best[0], 1)))
    return {"fixed": fixed, "unfixed": unfixed}


mesh = read_fort14(SRC)
arc0 = mesh.nodes[np.asarray(mesh.open_boundaries[0], int)].copy()

mesh, pinfo = align_open_boundary_local(mesh)
print(f"[fin] perp-local: remaining={pinfo.get('remaining')} "
      f"passes={pinfo.get('passes')}", flush=True)

mesh, hinfo = phase_h_finish(mesh, freeze_open_boundary=True)
print(f"[fin] phase_h_finish: "
      f"{ {k: v for k, v in hinfo.items() if not isinstance(v, (list, dict))} }",
      flush=True)

mesh, cinfo = compact_nodes(mesh)
print(f"[fin] compact_nodes: {cinfo}", flush=True)

finfo = flip_for_obc_perp(mesh)
print(f"[fin] obc flips: {finfo}", flush=True)


def fix_r4(mesh, min_angle=30.0, max_angle=130.0):
    """IN-PLACE: elements carrying an OBC edge whose third node is
    ALSO a boundary node (ISONB sum 5, FVCOM-fatal R4): flip one of
    the element's internal edges so the neighbour's interior node
    replaces the boundary third-node. The OBC edge itself is kept.
    Quality-gated."""
    nodes, els = mesh.nodes, mesh.elements
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    ee.sort(axis=1)
    uniq, cnt = np.unique(ee, axis=0, return_counts=True)
    bnd = np.zeros(len(nodes), bool)
    bnd[uniq[cnt == 1].ravel()] = True
    ob = set(int(v) for v in np.asarray(mesh.open_boundaries[0]))
    fixed, unfixed = [], []
    for ei in range(len(els)):
        tri = [int(x) for x in els[ei]]
        if not all(bnd[v] for v in tri):
            continue
        obc_edge = [(a, b) for a, b in ((tri[0], tri[1]),
                                        (tri[1], tri[2]),
                                        (tri[2], tri[0]))
                    if a in ob and b in ob]
        if not obc_edge:
            continue
        o1, o2 = obc_edge[0]
        w = [v for v in tri if v not in (o1, o2)][0]
        done = False
        for oo in (o1, o2):
            other = o2 if oo == o1 else o1
            nb = [int(ej) for ej in np.where(
                ((els == oo).any(axis=1)) & ((els == w).any(axis=1))
            )[0] if ej != ei]
            if not nb:
                continue
            ej = nb[0]
            m = int([x for x in els[ej] if int(x) not in (oo, w)][0])
            if bnd[m]:
                continue
            t1, t2 = [other, oo, m], [other, m, w]
            if _area(*nodes[t1]) < 0:
                t1 = [other, m, oo]
            if _area(*nodes[t2]) < 0:
                t2 = [other, w, m]
            A1, A2 = _area(*nodes[t1]), _area(*nodes[t2])
            if A1 <= 0 or A2 <= 0:
                continue
            ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
            # R4 is FVCOM-fatal, C1 is a quality target: accept
            # the flip even with a small (<=2.5 deg) min-angle
            # regression -- the far node is interior and the
            # phase_h pass right after polishes it back above 30
            cur = _tri_angles(*nodes[tri]) + _tri_angles(
                *nodes[[int(x) for x in els[ej]]])
            if (min(ang) < min(min_angle, min(cur) - 2.5)
                    or max(ang) > max(max_angle, max(cur) + 5.0)):
                continue
            els[ei] = t1
            els[ej] = t2
            fixed.append(ei)
            done = True
            break
        if not done:
            unfixed.append(ei)
    return {"fixed": fixed, "unfixed": unfixed}


rinfo = fix_r4(mesh)
print(f"[fin] r4 flips: {rinfo}", flush=True)

# polish any C1 residual the R4 flips left (land nodes on the
# straight artificial closure may slide along it; OBC frozen)
mesh, hinfo2 = phase_h_finish(mesh, freeze_open_boundary=True)
print(f"[fin] phase_h_finish(2): delta_total="
      f"{hinfo2.get('delta_total')}", flush=True)
mesh, cinfo2 = compact_nodes(mesh)
print(f"[fin] compact_nodes(2): {cinfo2}", flush=True)

# LAST: R4 cells whose only flip candidate fails the angle gate
# (single-internal-edge end cells): split that internal edge at the
# fraction (35-65%) that maximises the min angle over the 4 new
# sub-triangles. Deterministic; runs after all node-moving stages
# so nothing can undo it.
rinfo2 = fix_r4(mesh)
if rinfo2["unfixed"]:
    import dataclasses
    nodes2, els2, dep2 = mesh.nodes, mesh.elements, mesh.depths
    ee2 = np.vstack([els2[:, [0, 1]], els2[:, [1, 2]],
                     els2[:, [2, 0]]])
    ee2.sort(axis=1)
    uq2, ct2 = np.unique(ee2, axis=0, return_counts=True)
    bnd2 = np.zeros(len(nodes2), bool)
    bnd2[uq2[ct2 == 1].ravel()] = True
    ob2 = set(int(v) for v in np.asarray(mesh.open_boundaries[0]))
    done_ins = []
    for ei in list(rinfo2["unfixed"]):
        tri = [int(x) for x in els2[ei]]
        obe = [(a3, b3) for a3, b3 in ((tri[0], tri[1]),
                                       (tri[1], tri[2]),
                                       (tri[2], tri[0]))
               if a3 in ob2 and b3 in ob2]
        if not obe:
            continue
        o1, o2 = obe[0]
        w = [v for v in tri if v not in (o1, o2)][0]
        # the internal edge (has a neighbour element)
        cand = None
        for oo in (o1, o2):
            nb = [int(ej) for ej in np.where(
                ((els2 == oo).any(axis=1))
                & ((els2 == w).any(axis=1)))[0] if ej != ei]
            if nb:
                cand = (oo, o2 if oo == o1 else o1, nb[0])
        if cand is None:
            continue
        oo, other, ej = cand
        mfar = int([x for x in els2[ej]
                    if int(x) not in (oo, w)][0])
        # external neighbour areas for the C4 part of the score
        def _ext_area(a4, b4):
            nb4 = [int(ek) for ek in np.where(
                ((els2 == a4).any(axis=1))
                & ((els2 == b4).any(axis=1)))[0]
                if ek not in (ei, ej)]
            if not nb4:
                return None
            t4 = [int(x) for x in els2[nb4[0]]]
            return abs(_area(*nodes2[t4]))

        ext = {(other, oo): _ext_area(other, oo),
               (other, w): _ext_area(other, w),
               (oo, mfar): _ext_area(oo, mfar),
               (mfar, w): _ext_area(mfar, w)}
        best = None
        for fr in np.linspace(0.30, 0.75, 19):
            sN = (1 - fr) * nodes2[oo] + fr * nodes2[w]
            tris = [[other, oo, -1], [other, -1, w],
                    [oo, mfar, -1] if _area(
                        nodes2[oo], nodes2[mfar], sN) > 0
                    else [mfar, oo, -1],
                    [mfar, w, -1] if _area(
                        nodes2[mfar], nodes2[w], sN) > 0
                    else [w, mfar, -1]]
            angs = []
            areas4 = []
            ok = True
            for t3 in tris:
                P3 = [nodes2[v] if v >= 0 else sN for v in t3]
                A3 = _area(*P3)
                if A3 <= 0:
                    ok = False
                    break
                areas4.append(A3)
                angs += _tri_angles(*P3)
            if not ok:
                continue
            # C4 across the internal seams and vs external nbrs
            c4v = []
            for (ia, ib) in ((0, 1), (0, 2), (1, 3), (2, 3)):
                c4v.append(abs(areas4[ia] - areas4[ib])
                           / max(areas4[ia], areas4[ib]))
            for k4, key in enumerate(((other, oo), (other, w),
                                      (oo, mfar), (mfar, w))):
                Ae = ext[key]
                if Ae is not None:
                    c4v.append(abs(areas4[k4] - Ae)
                               / max(areas4[k4], Ae))
            feas = (min(angs) >= 30.0 and max(angs) <= 130.0
                    and max(c4v) <= 0.5)
            score = (1 if feas else 0, min(angs) - max(c4v))
            if best is None or score > best[0]:
                best = (score, fr, sN, tris)
        if best is None:
            print(f"[fin] r4 edge-split: NO candidate for elem "
                  f"{ei}", flush=True)
            continue
        score, fr, sN, tris = best
        mn = score[1]
        si = len(nodes2)
        nodes2 = np.vstack([nodes2, sN[None, :]])
        dep2 = np.append(dep2, 0.5 * (dep2[oo] + dep2[w]))
        tt = [[si if v < 0 else v for v in t3] for t3 in tris]
        els2 = els2.copy()
        els2[ei] = tt[0]
        els2[ej] = tt[2]
        els2 = np.vstack([els2, [tt[1]], [tt[3]]])
        done_ins.append((ei, round(fr, 2), round(mn, 1)))
        mesh = dataclasses.replace(mesh, nodes=nodes2,
                                   depths=dep2, elements=els2)
        nodes2, els2, dep2 = mesh.nodes, mesh.elements, mesh.depths
    print(f"[fin] r4 edge-split fallback: {done_ins}", flush=True)
    # polish the area-ratio seams the split leaves (moves only --
    # no vertex_remove, so the inserted node cannot be deleted)
    from fvcom_mesh_tools.mesh_clean_phase_h import (
        _stochastic_local_fix_round,
    )
    st = _stochastic_local_fix_round(
        mesh, np.random.default_rng(4242),
        min_angle_target=30.0, max_angle_target=130.0,
        area_ratio_target=0.5, max_valence=8,
        max_tries_per_fail=500, perturbation_sigma=0.3,
        max_outer_passes=5, coastline_projector=None,
        freeze_open_boundary=True)
    print(f"[fin] post-split polish: {st}", flush=True)


def flip_c4_edges(mesh, max_area_change=0.5, min_angle=30.0,
                  max_angle=130.0):
    """IN-PLACE: flip internal edges whose adjacent-element area
    ratio exceeds the C4 bound -- the swapped diagonal mixes the
    small and large triangle, rebalancing areas. Quality-gated."""
    nodes, els = mesh.nodes, mesh.elements
    fixed, unfixed = [], []
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, inv, ct = np.unique(ee, axis=0, return_inverse=True,
                            return_counts=True)
    for k in np.where(ct == 2)[0]:
        eids = np.where(inv == k)[0] % len(els)
        ei, ej = int(eids[0]), int(eids[1])
        A1 = abs(_area(*nodes[[int(x) for x in els[ei]]]))
        A2 = abs(_area(*nodes[[int(x) for x in els[ej]]]))
        if abs(A1 - A2) / max(A1, A2) <= max_area_change:
            continue
        a4, b4 = [int(v) for v in uq[k]]
        m1 = int([x for x in els[ei] if int(x) not in (a4, b4)][0])
        m2 = int([x for x in els[ej] if int(x) not in (a4, b4)][0])
        t1, t2 = [m1, a4, m2], [m1, m2, b4]
        if _area(*nodes[t1]) < 0:
            t1 = [m1, m2, a4]
        if _area(*nodes[t2]) < 0:
            t2 = [m1, b4, m2]
        B1, B2 = _area(*nodes[t1]), _area(*nodes[t2])
        if B1 <= 0 or B2 <= 0:
            unfixed.append((a4, b4))
            continue
        ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
        if (min(ang) < min_angle or max(ang) > max_angle
                or abs(B1 - B2) / max(B1, B2) > max_area_change):
            unfixed.append((a4, b4))
            continue
        els[ei] = t1
        els[ej] = t2
        fixed.append((a4, b4))
    return {"fixed": fixed, "unfixed": unfixed}


c4info = flip_c4_edges(mesh)
print(f"[fin] c4 flips: {c4info}", flush=True)

arc1 = mesh.nodes[np.asarray(mesh.open_boundaries[0], int)]
mv = float(np.hypot(*(arc1 - arc0).T).max())
print(f"[fin] OBC arc max displacement = {mv:.4f} m", flush=True)
if mv > 1e-3:
    raise RuntimeError(
        f"finishing moved constrained OBC nodes (max {mv:.3f} m). "
        "The input arc must stay fixed -- rerun with the offending "
        "stage disabled or pin the OBC nodes.")

write_fort14(mesh, DST)
print(f"[fin] wrote {DST} NP={mesh.n_nodes:,} NE={len(mesh.elements):,}",
      flush=True)
