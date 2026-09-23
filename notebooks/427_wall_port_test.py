# A split wall in a case KNOWN to run: the finished Kimitsu port mesh.
#
# The synthetic channel of notebook 426 went unstable beside its open
# boundary in BOTH its split and unsplit form, so it could not tell a wall
# from a forcing problem.  This changes one thing only in a case that has
# already integrated for 20 days: three paths of existing interior edges are
# split into walls, and nothing else moves -- the OBC node ids, the sponge,
# the tide, the namelist and the element order are all the staged case's.
#
#   seal        coast -> water -> coast: cuts a piece of the domain off, and a
#               piece with no open boundary must keep its level at zero
#   pier        coast -> an interior node: one root, one free tip
#   breakwater  interior -> interior: two free tips
#
#   python notebooks/427_wall_port_test.py prep <staged refined case dir> <out dir>
#   python notebooks/427_wall_port_test.py analyze <out dir> <unsplit run dir>
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402
from fvcom_mesh_tools.walls import split_along_walls, wall_edges_from_path  # noqa: E402

CENTRE = np.array([393010.0, 3909480.0])      # the port region, UTM 54N


def edges_of(tri):
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return u, c


def path_between(xy, interior, a, b, banned):
    """Shortest path over INTERIOR edges, avoiding the banned nodes."""
    w = np.linalg.norm(xy[interior[:, 0]] - xy[interior[:, 1]], axis=1)
    ok = ~(np.isin(interior[:, 0], list(banned)) | np.isin(interior[:, 1], list(banned)))
    e, w = interior[ok], w[ok]
    n = len(xy)
    g = coo_matrix((np.r_[w, w], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                   shape=(n, n)).tocsr()
    _, pred = dijkstra(g, indices=a, return_predecessors=True)
    if pred[b] < 0:
        raise ValueError(f"no interior path from {a} to {b}")
    path = [b]
    while path[-1] != a:
        path.append(int(pred[path[-1]]))
    return path[::-1]


def prep(case: Path, out: Path):
    inp = case / "input"
    m = read_fvcom_case(inp / "m2_grd.dat", inp / "m2_dep.dat", inp / "m2_obc.dat")
    xy, tri = m.nodes[:, :2], m.elements
    u, c = edges_of(tri)
    interior, bnd = u[c == 2], u[c == 1]
    bnodes = np.unique(bnd)
    d = np.linalg.norm(xy - CENTRE, axis=1)
    near_b = bnodes[d[bnodes] < 1200.0]
    inner = np.setdiff1d(np.flatnonzero(d < 1200.0), bnodes)

    # seal: two coast nodes 700-1100 m apart; the path between them bulges
    # into the water and cuts off whatever lies between it and the coast
    best = None
    for i in near_b[::3]:
        dd = np.linalg.norm(xy[near_b] - xy[i], axis=1)
        for j in near_b[(dd > 700) & (dd < 1100)][:5]:
            try:
                p = path_between(xy, interior, int(i), int(j), set(bnodes) - {int(i), int(j)})
            except ValueError:
                continue
            if best is None or len(p) < len(best):
                best = p
    if best is None:
        raise SystemExit("no coast-to-coast interior path found")
    seal = best
    used = set(seal)
    # pier: a coast node well away from the seal, to an interior node ~250 m out
    far = near_b[np.linalg.norm(xy[near_b] - xy[seal[len(seal) // 2]], axis=1) > 600]
    pier = None
    for r in far:
        cand = inner[(np.abs(np.linalg.norm(xy[inner] - xy[r], axis=1) - 250) < 40)]
        cand = [q for q in cand if q not in used]
        if not cand:
            continue
        try:
            p = path_between(xy, interior, int(r), int(cand[0]),
                             (set(bnodes) | used) - {int(r)})
        except ValueError:
            continue
        pier = p
        break
    used |= set(pier)
    # breakwater: two interior nodes ~300 m apart, touching nothing
    brk = None
    free = [q for q in inner if q not in used]
    for a in free[::7]:
        cand = [q for q in free if abs(np.linalg.norm(xy[q] - xy[a]) - 300) < 40]
        if not cand:
            continue
        try:
            p = path_between(xy, interior, int(a), int(cand[0]), set(bnodes) | used)
        except ValueError:
            continue
        if not (set(p) & used):
            brk = p
            break
    walls = {"seal": seal, "pier": pier, "breakwater": brk}
    ew = np.vstack([wall_edges_from_path(p) for p in walls.values()])
    sxy, stri, copy_of, rep = split_along_walls(xy, tri, ew)
    print("[427] walls: " + ", ".join(f"{k} {len(v) - 1} edges" for k, v in walls.items()))
    print(f"[427] split: {json.dumps({k: v for k, v in rep.items() if k != 'pairs'})}")
    if rep["n_components"] != 2:
        raise SystemExit("the seal did not cut the domain in two")

    # which component carries the open boundary; the other is the sealed one
    obc = np.asarray(m.open_boundaries[0])
    comp = _components(stri, len(sxy))
    sealed_comp = [k for k in set(comp.values()) if k != comp[int(obc[0])]][0]
    sealed = np.array(sorted(n for n, k in comp.items() if k == sealed_comp))
    print(f"[427] sealed region: {len(sealed)} node(s)")

    if out.exists():
        raise SystemExit(f"exists: {out}")
    shutil.copytree(case, out, ignore=shutil.ignore_patterns("output", "*.log"))
    (out / "output").mkdir()
    oi = out / "input"
    # the staged file's own layout: "k a b c" per cell, "k x y" per node
    with (oi / "m2_grd.dat").open("w") as f:
        f.write(f"Node Number = {len(sxy)}\nCell Number = {len(stri)}\n")
        for k, (a, b, cc) in enumerate(stri + 1, 1):
            f.write(f"{k} {a} {b} {cc}\n")
        for k, (x, y) in enumerate(sxy, 1):
            f.write(f"{k} {x:.8f} {y:.8f}\n")
    dep = m.depths[copy_of]
    (oi / "m2_dep.dat").write_text(f"Node Number = {len(sxy)}\n" + "".join(
        f"{x:.8f} {y:.8f} {h:.6f}\n" for (x, y), h in zip(sxy, dep)))
    cor_lines = (inp / "m2_cor.dat").read_text().splitlines()
    lat = np.array([float(ln.split()[-1]) for ln in cor_lines[1:] if ln.strip()])
    (oi / "m2_cor.dat").write_text(f"Node Number = {len(sxy)}\n" + "".join(
        f"{x:.8f} {y:.8f} {la:.6f}\n" for (x, y), la in zip(sxy, lat[copy_of])))
    nml = (out / "m2_run.nml").read_text().replace(str(case), str(out))
    (out / "m2_run.nml").write_text(nml)
    (out / "walls.json").write_text(json.dumps({
        "walls": walls, "split": {k: v for k, v in rep.items() if k != "pairs"},
        "pairs": rep["pairs"], "sealed_nodes_split_ids": sealed.tolist(),
        "sealed_nodes_original_ids": sorted(set(copy_of[sealed].tolist())),
        "n_original_nodes": int(len(xy))}, indent=1))


def _components(tri, n):
    parent = list(range(n))

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for a, b, c in tri.tolist():
        for x, y in ((a, b), (b, c)):
            ra, rb = find(x), find(y)
            if ra != rb:
                parent[ra] = rb
    return {k: find(k) for k in range(n)}


def m2_amp(t, z, period):
    w = 2 * np.pi / period
    a = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(a, z, rcond=None)
    return np.hypot(coef[1], coef[2])


def analyze(out: Path, control: Path):
    import netCDF4

    meta = json.loads((out / "walls.json").read_text())
    period = 12.4206012 * 3600.0
    res = {}
    for name, run, ids in (("split", out, meta["sealed_nodes_split_ids"]),
                           ("unsplit", control, meta["sealed_nodes_original_ids"])):
        with netCDF4.Dataset(run / "output" / "m2_0001.nc") as ds:
            t = np.asarray(ds["time"][:], float) * 86400.0
            z = np.asarray(ds["zeta"][:], float)
            sp = np.hypot(np.asarray(ds["ua"][:], float), np.asarray(ds["va"][:], float))
            nv = np.asarray(ds["nv"][:], int).T - 1
        keep = t >= t[0] + 86400.0
        amp = m2_amp(t[keep], z[keep], period)
        ids = np.asarray(ids)
        others = np.setdiff1d(np.arange(z.shape[1]), ids)
        r = {"records": int(len(t)),
             "finite": bool(np.isfinite(z).all() and np.isfinite(sp).all()),
             "sealed_zeta_max_abs_m": float(np.abs(z[:, ids]).max()),
             "sealed_m2_amp_max_m": float(amp[ids].max()),
             "open_m2_amp_median_m": float(np.median(amp[others])),
             "speed_p99_ms": float(np.percentile(sp[keep], 99)),
             "speed_max_ms": float(sp[keep].max())}
        if name == "split":
            tips = meta["split"]["free_tips"]
            at_tip = np.isin(nv, tips).any(axis=1)
            r["tip_speed_max_ms"] = float(sp[keep][:, at_tip].max())
        res[name] = r
        print(f"[427] {name}: {json.dumps(r)}")
    s, u = res["split"], res["unsplit"]
    verdict = {
        "split_sealed_region_is_still": s["sealed_m2_amp_max_m"] < 1e-4,
        "unsplit_same_region_sees_the_tide":
            u["sealed_m2_amp_max_m"] > 0.5 * u["open_m2_amp_median_m"],
        "split_run_is_finite": s["finite"],
        "tips_are_not_a_hotspot": s["tip_speed_max_ms"] <= 3 * s["speed_p99_ms"],
        "outside_tide_unchanged_mm": 1000 * (s["open_m2_amp_median_m"] - u["open_m2_amp_median_m"]),
    }
    print(f"[427] VERDICT {json.dumps(verdict)}")
    (out / "result.json").write_text(json.dumps({"cases": res, "verdict": verdict}, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "prep":
        prep(Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve())
    else:
        analyze(Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve())
