"""``fmesh-refine-depths``: the step after a hires refinement.

The hires branch writes the depths the source gives, and says so: nothing is
floored, capped or smoothed there, so a tidal flat comes out at -4.54 m and
the case is not a finished bathymetry (owner, 2026-09-23). This is the step
that finishes it, and it is separate for a reason -- a refinement costs a
seed search over the whole fill-repair-stitch loop, and trying a different
minimum depth must not cost that again.

Three operations, in the order the production recipe uses them:

1. **floor** at ``--hmin``. On a nori ground this is the dominant one: inside
   `banzu_nori` 22 % of the M7001 soundings are shallower than 3 m.
2. **cap** at ``--hmax``.
3. **r-factor smoothing**, over the patch only, with every frozen depth held
   fixed. This is the "smooth last, including the transition" of the original
   request, and it is here rather than in the refinement because a floor is
   what makes ``|hi-hj| / (hi+hj)`` defined on every edge: before it, a node
   above the datum makes that expression a sign error dressed as a gradient.

**No frozen depth moves.** ``node_map`` says which nodes the patch created,
and only those are movable. The report says how far each moved and what the
limiter could not fix -- an edge between two frozen nodes is unfixable here,
and one the patch created between two retained nodes is the patch's.

    fmesh-refine-depths outputs/refine_kimitsu_port_hires --hmin 3 --rfactor 0.2
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.refine import limit_rfactor


def _read_cor(path: Path) -> np.ndarray:
    rows = [ln.split() for ln in path.read_text().splitlines()[1:] if ln.strip()]
    return np.array([float(r[-1]) for r in rows], dtype=float)


def finish_depths(depths, elements, movable, *, hmin, hmax, rfactor):
    """Floor, cap, then smooth -- and report each separately.

    The three are reported apart because they are three different claims
    about the seabed, and a single "depths changed" number hides which one
    did it. Stage 4 can also drive a node onto a bound that stage 1 did not,
    so the clipping counts are taken twice.
    """
    h0 = np.asarray(depths, dtype=float)
    if not np.isfinite(h0).all():
        raise ValueError("the depths to finish are not all finite")
    if not (np.isfinite(hmin) and hmin > 0):
        raise ValueError("hmin must be finite and positive")
    if not (np.isfinite(hmax) and hmax > hmin):
        raise ValueError("hmax must be finite and greater than hmin")
    free = np.asarray(movable, dtype=bool)

    clipped = h0.copy()
    lo = free & (h0 < hmin)
    hi = free & (h0 > hmax)
    clipped[lo] = hmin
    clipped[hi] = hmax
    report = {
        "hmin_m": float(hmin), "hmax_m": float(hmax), "rfactor": float(rfactor),
        "n_movable": int(free.sum()),
        "n_floored": int(lo.sum()), "n_capped": int(hi.sum()),
        "floored_fraction_of_movable": float(lo.sum() / max(free.sum(), 1)),
        "raw_min_m": float(h0.min()), "raw_max_m": float(h0.max()),
        "n_below_zero_before": int((h0 <= 0).sum()),
    }
    # A frozen depth below the floor is the BASE's, and this may not touch
    # it: the contract is that the patch changes nothing outside itself, and
    # a base that sits on its own floor is not this step's business.
    report["n_frozen_below_hmin"] = int((~free & (h0 < hmin)).sum())
    # But a frozen depth at or BELOW ZERO stops the whole step, because the
    # r-factor is not defined on it and this may not move it.  Said here, with
    # the count, rather than left to surface as the limiter's own message
    # about a mesh the caller did not choose.
    used = np.unique(np.asarray(elements, dtype=np.int64))
    dry_frozen = (~free & (clipped <= 0.0))[used] if used.size else np.zeros(0, bool)
    if dry_frozen.any():
        raise ValueError(
            f"{int(dry_frozen.sum())} FROZEN node(s) are at or above the datum. "
            "The r-factor is not defined there and this step may not move a "
            "frozen depth, so the base mesh itself would have to be floored "
            "first -- which is a change to the model, not to the patch")

    out, rinfo = limit_rfactor(elements, clipped, free, float(rfactor),
                               depth_min=float(hmin), depth_max=float(hmax))
    report["rfactor_report"] = rinfo
    # `limit_rfactor` judges at 1e-9, which is tighter than the depth file it
    # will be written to.  `TokyoBay_dep_m7001tp_rfac0p2_cap300.dat` has a
    # max r of exactly 0.2000 and 497 of its 8,858 edges sit ON that bound;
    # written to six decimals they come back as 0.2 + 3e-8, and the patch was
    # reported NOT CONVERGED for 374 edges it is forbidden to touch.  So the
    # honest count is taken again at the file's own precision, and both go in
    # the report rather than one replacing the other.
    tri = np.asarray(elements, dtype=np.int64)
    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                     tri[:, [2, 0]]]), axis=1), axis=0)
    r = np.abs(out[e[:, 0]] - out[e[:, 1]]) / (out[e[:, 0]] + out[e[:, 1]])
    tol = float(rfactor) * 1e-6 + 1e-9
    touched = free[e].any(axis=1)
    over = r > float(rfactor) + tol
    report["write_precision_tolerance"] = tol
    report["max_r_movable"] = float(r[touched].max()) if touched.any() else 0.0
    report["max_r_frozen_pair"] = float(r[~touched].max()) if (~touched).any() else 0.0
    report["n_over_movable_at_tolerance"] = int((over & touched).sum())
    report["n_over_frozen_pair_at_tolerance"] = int((over & ~touched).sum())
    report["converged_at_write_precision"] = bool(not (over & touched).any())
    report["n_floored_after_smoothing"] = int((free & (out <= hmin + 1e-9)).sum())
    report["n_capped_after_smoothing"] = int((free & (out >= hmax - 1e-9)).sum())
    moved = np.abs(out - h0)
    report["max_change_m"] = float(moved.max())
    report["max_frozen_change_m"] = float(moved[~free].max()) if (~free).any() else 0.0
    report["final_min_m"] = float(out.min())
    report["final_max_m"] = float(out.max())
    return out, report


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-refine-depths",
        description="Floor, cap and smooth a hires refinement's depths.")
    p.add_argument("refine_out", type=Path,
                   help="a refinement output directory (holds node_map.npy "
                        "and fvcom/)")
    p.add_argument("--hmin", type=float, required=True, help="minimum depth, m")
    p.add_argument("--hmax", type=float, default=300.0, help="maximum depth, m")
    p.add_argument("--rfactor", type=float, default=0.2,
                   help="Beckmann-Haidvogel limit on |dh|/(hi+hj)")
    p.add_argument("--outdir", type=Path, default=None,
                   help="where to write (default: <refine_out>/fvcom_finished)")
    p.add_argument("--dry-run", action="store_true", help="report, write nothing")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.refine_out.resolve()
    nm_path = root / "node_map.npy"
    fvdir = root / "fvcom"
    if not nm_path.exists() or not fvdir.is_dir():
        print(f"not a refinement output directory: {root}", file=sys.stderr)
        return 2
    grd = next(fvdir.glob("*_grd.dat"), None)
    if grd is None:
        print(f"no _grd.dat in {fvdir}", file=sys.stderr)
        return 2
    name = grd.name[: -len("_grd.dat")]
    dep = fvdir / f"{name}_dep.dat"
    obc = fvdir / f"{name}_obc.dat"
    mesh = read_fvcom_case(grd, dep, obc if obc.exists() else None)

    node_map = np.load(nm_path)
    movable = np.ones(mesh.n_nodes, dtype=bool)
    movable[node_map[node_map >= 0]] = False
    print(f"[finish] {name}: {mesh.n_nodes:,} nodes, "
          f"{int(movable.sum()):,} the patch created")

    depths, rep = finish_depths(mesh.depths, mesh.elements, movable,
                                hmin=args.hmin, hmax=args.hmax,
                                rfactor=args.rfactor)
    r = rep["rfactor_report"]
    print(f"[finish] floor {args.hmin:g} m: {rep['n_floored']:,} node(s) "
          f"({100 * rep['floored_fraction_of_movable']:.0f} % of the patch), "
          f"cap {args.hmax:g} m: {rep['n_capped']:,}")
    print(f"[finish] r <= {args.rfactor:g}: {r['n_depths_changed']:,} depth(s) "
          f"moved, worst {r['max_depth_change_m']:.2f} m in {r['rounds']} "
          f"round(s). At the depth file's own precision "
          f"({rep['write_precision_tolerance']:.1e}): "
          f"{'CONVERGED' if rep['converged_at_write_precision'] else 'NOT CONVERGED'}"
          f", {rep['n_over_movable_at_tolerance']} movable and "
          f"{rep['n_over_frozen_pair_at_tolerance']} frozen-pair edge(s) over; "
          f"max r movable {rep['max_r_movable']:.6f}, frozen "
          f"{rep['max_r_frozen_pair']:.6f}")
    print(f"[finish] depths {rep['raw_min_m']:.2f}..{rep['raw_max_m']:.2f} -> "
          f"{rep['final_min_m']:.2f}..{rep['final_max_m']:.2f} m; "
          f"frozen moved {rep['max_frozen_change_m']:.3g} m")
    if rep["max_frozen_change_m"] > 0:
        print("[finish] a frozen depth moved -- that is the contract broken",
              file=sys.stderr)
        return 1
    if args.dry_run:
        return 0

    outdir = args.outdir or root / "fvcom_finished"
    outdir.mkdir(parents=True, exist_ok=True)
    cor_path = fvdir / f"{name}_cor.dat"
    cor = _read_cor(cor_path) if cor_path.exists() else None
    import dataclasses

    written = export_fvcom_case(
        dataclasses.replace(mesh, depths=depths), outdir, name, cor=cor,
        obc_type=getattr(mesh, "obc_type", 1), obc_depth_control=False)
    (outdir / f"{name}_finish.json").write_text(
        json.dumps(rep, indent=1, default=float) + "\n")
    print(f"[finish] wrote {', '.join(sorted(written))} in {outdir}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
