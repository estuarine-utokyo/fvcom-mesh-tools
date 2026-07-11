# Find duplicate NON-degenerate edges in the SDF polygon set and
# verify total-crossing parity for the broken rows.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh import Shoreline
from oceanmesh import edges as om_edges

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox_poly = np.array([
    [-71.6, 42.7], [-64, 30], [-80, 24], [-85, 38], [-71.6, 42.7]])
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox_poly, 1e3*DEG)
for label, arrs in (("outer+inner", (shore.boubox, shore.inner)),):
    poly = np.vstack([np.asarray(a) for a in arrs if len(a)])
    e = om_edges.get_poly_edges(poly)
    pn = np.nan_to_num(poly)
    a = pn[e[:, 0]]; b = pn[e[:, 1]]
    L = np.hypot(*(a - b).T)
    key = np.round(
        np.sort(np.stack([a, b], axis=1), axis=1).reshape(len(e), 4),
        12,
    )
    uniq, first_idx, counts = np.unique(
        key, axis=0, return_index=True, return_counts=True)
    dup = counts > 1
    nz = uniq[dup]
    # report only non-degenerate duplicates
    ndeg = (np.hypot(nz[:, 2] - nz[:, 0], nz[:, 3] - nz[:, 1])
            > 1e-12)
    print(f"[dup] {label}: edges={len(e)} zero-len="
          f"{int((L<1e-12).sum())} dup-classes={int(dup.sum())} "
          f"non-degenerate-dups={int(ndeg.sum())}", flush=True)
    for row in nz[ndeg][:12]:
        print(f"  DUP edge ({row[0]:.8f},{row[1]:.8f})->"
              f"({row[2]:.8f},{row[3]:.8f}) "
              f"y-span=({min(row[1],row[3]):.6f},"
              f"{max(row[1],row[3]):.6f})", flush=True)
    # total crossings parity on the four broken rows
    for qy in (25.086886886886887, 25.54334334334334,
               25.801601601601602, 25.9978):
        ca = a[:, 1] > qy
        cb = b[:, 1] > qy
        cross = ca != cb
        print(f"[dup] row y={qy:.6f}: total crossings ="
              f" {int(cross.sum())} (parity "
              f"{'ODD-BAD' if cross.sum() % 2 else 'even-ok'})",
              flush=True)
