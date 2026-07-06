"""Stage-2 M3: declarative local-resolution directives (the
SMS-style edit, as recipe data). Each directive re-meshes a
polygon via oceanmesh.remesh_patch at target_h_m; per-directive
atomic with local acceptance, OBC protected, depths re-interpolated
for new nodes."""

from __future__ import annotations

import numpy as np

from .detect import detect_violations


def apply_directives(mesh, directives, utm_epsg=32654, log=print):
    """Apply recipe ``finishing.directives`` to a Fort14Mesh whose
    nodes are in EPSG:``utm_epsg``. Polygons are lon/lat. Returns
    (mesh, ledger). Topology changes; open-boundary node ids are
    remapped by exact coordinates (merge keeps them verbatim)."""
    import shapely
    from oceanmesh import remesh_patch
    from pyproj import Transformer
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    from scipy.spatial import cKDTree
    from shapely.geometry import Polygon

    tr = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                              always_xy=True)
    ledger = []
    for k, d in enumerate(directives or []):
        rec = {"directive": k, "op": d.get("op", "refine"),
               "target_h_m": d.get("target_h_m")}
        ring = np.asarray(d["polygon"], float)
        x, y = tr.transform(ring[:, 0], ring[:, 1])
        poly_utm = np.column_stack([x, y])
        pg = Polygon(poly_utm)
        obc_ids = [int(v) for b in mesh.open_boundaries for v in b]
        if obc_ids and shapely.contains_xy(
                pg, mesh.nodes[obc_ids, 0],
                mesh.nodes[obc_ids, 1]).any():
            rec["outcome"] = "skipped (contains OBC nodes)"
            log(f"[finishing] directive {k}: {rec['outcome']}")
            ledger.append(rec)
            continue
        P0, T0, D0 = mesh.nodes, mesh.elements, mesh.depths
        try:
            P1, T1 = remesh_patch(
                P0, T0, poly_utm,
                target_h=float(d["target_h_m"]), seed=42 + k,
            )
        except Exception as err:  # atomic: keep the old mesh
            rec["outcome"] = f"failed ({err})"
            log(f"[finishing] directive {k}: {rec['outcome']}")
            ledger.append(rec)
            continue
        # gross-failure gate only: borderline seam angles
        # (20-30 deg) are left for the auto-repair stage that runs
        # AFTER directives (design: "auto-repair heals their
        # seams"). Revert on angles < 20 deg or mass breakage.
        det = detect_violations(P1, T1,
                                thresholds={"c1_min_deg": 20.0})
        cen = P1[T1].mean(axis=1)
        halo = pg.buffer(3.0 * float(d["target_h_m"]))
        gross = [ie for c in ("c1", "c2")
                 for ie in det[c]["elements"]
                 if shapely.contains_xy(halo, cen[ie, 0],
                                        cen[ie, 1])]
        if len(gross) > 0:
            rec["outcome"] = (f"reverted ({len(gross)} gross "
                              "local failures)")
            log(f"[finishing] directive {k}: {rec['outcome']}")
            ledger.append(rec)
            continue
        # depths for new nodes; obc id remap by exact coordinates
        lin = LinearNDInterpolator(P0, D0)
        near = NearestNDInterpolator(P0, D0)
        D1 = lin(P1)
        miss = ~np.isfinite(D1)
        if miss.any():
            D1[miss] = near(P1[miss])
        tree = cKDTree(P1)
        new_obc = []
        ok = True
        for b in mesh.open_boundaries:
            dd, jj = tree.query(P0[np.asarray(b, int)])
            if dd.max() > 1e-6:
                ok = False
                break
            new_obc.append(jj.astype(int))
        if not ok:
            rec["outcome"] = "reverted (OBC nodes not preserved)"
            log(f"[finishing] directive {k}: {rec['outcome']}")
            ledger.append(rec)
            continue
        mesh.nodes, mesh.elements, mesh.depths = P1, T1, D1
        mesh.open_boundaries = new_obc
        rec["outcome"] = (f"applied (NE {len(T0)} -> {len(T1)})")
        log(f"[finishing] directive {k}: {rec['outcome']}")
        ledger.append(rec)
    return mesh, ledger
