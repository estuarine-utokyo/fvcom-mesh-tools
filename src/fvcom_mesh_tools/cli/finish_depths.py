"""``fmesh-finish-depths``: floor, cap and r-factor-smooth an FVCOM depth file.

The same three operations that make the run-ready depth files of
``TB-FVCOM/input/goto2023/grid`` -- ``TokyoBay_dep_m7001tp_rfac0p2_cap300.dat``
is "M7001, floored at 3 m, smoothed to r <= 0.2, capped at 300 m" -- with the
three numbers as options:

    --hmin     the minimum depth, m           (``min<N>m`` in the file name)
    --hmax     the maximum depth, m           (``cap<NNN>``)
    --rfactor  the Beckmann-Haidvogel limit   (``rfac<r>``; 0 skips smoothing)
               on |h_i - h_j| / (h_i + h_j) over every edge, the usual
               condition for a stable sigma-coordinate pressure gradient

The output is named the TB-FVCOM way, ``<case>_dep_min3m_rfac0p2_cap300.dat``,
next to a JSON report of what each operation moved.

Two kinds of input:

* **a refinement output directory** (``fmesh-refine``'s ``outputs/refine_*``).
  Only the nodes the patch created move; every depth of the frozen base mesh
  is kept bit for bit, which is the refinement's contract. ``--whole-mesh``
  lifts that and processes every node. The full case (grd, dep, obc, cor) is
  also written to ``<dir>/fvcom_finished`` for staging an FVCOM run.
* **an FVCOM case**, given as its grd file (``--dep`` picks the depth file
  to start from; default ``<case>_dep.dat``). Every node moves.

The smoothing is TB-FVCOM's own (``--method equal``, the default): floor,
smooth the uncapped field with equal-and-opposite moves, then cap -- from
``TokyoBay_dep_m7001tp_raw.dat`` it reproduces
``TokyoBay_dep_m7001tp_rfac0p2_cap300.dat``. ``--method limit`` is the
refinement's own limiter (floor and cap first, then pull each edge to the
bound), which moves fewer depths but gives a different field.

Examples::

    fmesh-finish-depths outputs/refine_kimitsu_port_hires --hmin 3 --hmax 300 --rfactor 0.2
    fmesh-finish-depths ~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat \\
        --dep TokyoBay_dep_m7001tp_raw.dat --hmin 3 --hmax 300 --rfactor 0.2 --outdir /tmp/dep
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.cli.refine_depths import _read_cor, finish_depths
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case


def _num(x: float) -> str:
    """3 -> '3', 2.5 -> '2p5', 0.2 -> '0p2' (the TB-FVCOM spelling)."""
    return f"{x:g}".replace(".", "p")


def variant_tag(hmin: float, hmax: float, rfactor: float | None) -> str:
    """``min3m_rfac0p2_cap300`` -- floor, smoothing, cap, in that order."""
    parts = [f"min{_num(hmin)}m"]
    if rfactor is not None and rfactor > 0:
        parts.append(f"rfac{_num(rfactor)}")
    parts.append(f"cap{int(round(hmax)):03d}" if float(hmax).is_integer()
                 else f"cap{_num(hmax)}")
    return "_".join(parts)


def write_dep(path: Path, nodes: np.ndarray, depths: np.ndarray) -> None:
    """An FVCOM ``_dep.dat``: a node-count header, then ``x y depth`` rows."""
    with Path(path).open("w") as f:
        f.write(f"Node Number = {len(depths)}\n")
        for (x, y), h in zip(np.asarray(nodes)[:, :2], depths):
            f.write(f"{x:.6f} {y:.6f} {h:.6f}\n")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-finish-depths",
        description="Floor, cap and r-factor-smooth FVCOM depths; write a "
                    "TB-FVCOM-style <case>_dep_<variant>.dat.")
    p.add_argument("source", type=Path,
                   help="a refinement output directory, or an FVCOM *_grd.dat")
    p.add_argument("--hmin", type=float, required=True, help="minimum depth, m")
    p.add_argument("--hmax", type=float, default=300.0, help="maximum depth, m (default 300)")
    p.add_argument("--rfactor", type=float, default=0.2,
                   help="r-factor limit |dh|/(hi+hj) (default 0.2; 0 = no smoothing)")
    p.add_argument("--dep", type=Path, default=None,
                   help="FVCOM case only: the depth file to start from "
                        "(a name next to the grd file, or a path)")
    p.add_argument("--whole-mesh", action="store_true",
                   help="refinement directory only: also move the frozen base depths")
    p.add_argument("--tag", default=None,
                   help="variant name (default min<hmin>m[_rfac<r>]_cap<hmax>)")
    p.add_argument("--outdir", type=Path, default=None,
                   help="where to write (default: <dir>/fvcom_finished, or the "
                        "grd file's own directory)")
    p.add_argument("--method", choices=("equal", "limit"), default="equal",
                   help="equal (default): TB-FVCOM's smoother -- floor, smooth the "
                        "uncapped field with equal-and-opposite moves, then cap; "
                        "limit: the refinement's limiter, floor and cap first")
    p.add_argument("--rounds", type=int, default=2000,
                   help="smoothing sweeps allowed (default 2000)")
    p.add_argument("--allow-unconverged", action="store_true",
                   help="write the product even if the r-factor limit was not "
                        "reached (exit status stays 3); for diagnosis only")
    p.add_argument("--dry-run", action="store_true", help="report, write nothing")
    return p


def main(argv: list[str] | None = None) -> int:
    """Exit 0 on a run-ready product, 2 on bad input, 3 when not converged."""
    args = build_parser().parse_args(argv)
    # Only 0 (no smoothing) or a finite r in (0, 1) means anything: NaN and a
    # negative r used to be read as "no smoothing" and 1.2 or inf as a limit
    # every field already meets, all with exit 0 (review F8).
    if not (np.isfinite(args.rfactor) and 0.0 <= args.rfactor < 1.0):
        print(f"--rfactor must be 0 (no smoothing) or a finite number in (0, 1); "
              f"got {args.rfactor}", file=sys.stderr)
        return 2
    if args.rfactor > 0 and args.rounds < 1:
        print("--rounds must be at least 1 when smoothing", file=sys.stderr)
        return 2
    src = args.source.expanduser().resolve()
    refinement = src.is_dir()
    if refinement:
        fvdir = src / "fvcom"
        grd = next(fvdir.glob("*_grd.dat"), None) if fvdir.is_dir() else None
        if grd is None or not (src / "node_map.npy").exists():
            print(f"not a refinement output directory: {src}", file=sys.stderr)
            return 2
        name = grd.name[: -len("_grd.dat")]
        dep = fvdir / f"{name}_dep.dat"
    else:
        grd = src
        if not grd.name.endswith("_grd.dat") or not grd.exists():
            print(f"not an FVCOM *_grd.dat: {grd}", file=sys.stderr)
            return 2
        name = grd.name[: -len("_grd.dat")]
        dep = grd.parent / f"{name}_dep.dat" if args.dep is None else (
            args.dep if args.dep.is_absolute() or args.dep.exists()
            else grd.parent / args.dep)
    dep = Path(dep).expanduser().resolve()
    if not dep.exists():
        print(f"no depth file: {dep}", file=sys.stderr)
        return 2
    obc = grd.parent / f"{name}_obc.dat"
    mesh = read_fvcom_case(grd, dep, obc if obc.exists() else None)

    movable = np.ones(mesh.n_nodes, dtype=bool)
    if refinement and not args.whole_mesh:
        node_map = np.load(src / "node_map.npy")
        movable[node_map[node_map >= 0]] = False
    scope = "patch nodes only" if not movable.all() else "every node"
    rf = args.rfactor if args.rfactor > 0 else None
    tag = args.tag or variant_tag(args.hmin, args.hmax, rf)
    print(f"[finish] {name} from {dep.name}: {mesh.n_nodes:,} nodes, "
          f"{int(movable.sum()):,} movable ({scope}); variant {tag}")
    try:
        depths, rep = finish_depths(mesh.depths, mesh.elements, movable,
                                    hmin=args.hmin, hmax=args.hmax, rfactor=rf,
                                    rounds=args.rounds, method=args.method)
    except ValueError as exc:
        print(f"[finish] {exc}", file=sys.stderr)
        return 2
    r = rep["rfactor_report"]
    print(f"[finish] floor {args.hmin:g} m: {rep['n_floored']:,} node(s); "
          f"cap {args.hmax:g} m: {rep['n_capped']:,}")
    if rf is None:
        print("[finish] no r-factor smoothing (--rfactor 0)")
    else:
        print(f"[finish] r <= {rf:g}: {r['n_depths_changed']:,} depth(s) moved, "
              f"worst {r['max_depth_change_m']:.2f} m in {r['rounds']} sweep(s); "
              f"{'CONVERGED' if rep['converged_at_write_precision'] else 'NOT CONVERGED'} "
              f"at the file's precision, max r {rep['max_r_movable']:.6f} "
              f"(movable), {rep['max_r_frozen_pair']:.6f} (frozen pairs)")
    print(f"[finish] depths {rep['raw_min_m']:.2f}..{rep['raw_max_m']:.2f} -> "
          f"{rep['final_min_m']:.2f}..{rep['final_max_m']:.2f} m")
    if rep["max_frozen_change_m"] > 0:
        print("[finish] a frozen depth moved -- the refinement contract is broken",
              file=sys.stderr)
        return 1
    converged = bool(rep["converged_at_write_precision"])
    if not converged:
        # Not a run-ready product, and not to be staged as one: exit 3 and
        # write nothing, unless asked to for diagnosis (review F7).
        print(f"[finish] NOT CONVERGED: the r-factor limit is not met on "
              f"{rep.get('n_over_movable_at_tolerance', '?')} edge(s) with a movable "
              f"end (max r {rep.get('max_r_movable', float('nan')):.6f}). Raise "
              "--rounds, or look for a frozen depth the floor cannot reach"
              + ("" if args.allow_unconverged else "; nothing written"),
              file=sys.stderr)
        if not args.allow_unconverged:
            return 3
    if args.dry_run:
        return 0 if converged else 3

    outdir = (args.outdir or (src / "fvcom_finished" if refinement else grd.parent)).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    variant = outdir / f"{name}_dep_{tag}.dat"
    write_dep(variant, mesh.nodes, depths)
    rep.update({"source_depth": str(dep), "variant": tag, "scope": scope})
    (outdir / f"{name}_dep_{tag}.json").write_text(
        json.dumps(rep, indent=1, default=float) + "\n")
    print(f"[finish] wrote {variant}")
    if refinement:
        cor_path = src / "fvcom" / f"{name}_cor.dat"
        cor = _read_cor(cor_path) if cor_path.exists() else None
        written = export_fvcom_case(
            dataclasses.replace(mesh, depths=depths), outdir, name, cor=cor,
            obc_type=getattr(mesh, "obc_type", 1), obc_depth_control=False)
        print(f"[finish] wrote the case ({', '.join(sorted(written))}) in {outdir}")
    return 0 if converged else 3


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
