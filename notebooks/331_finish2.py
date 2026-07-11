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
