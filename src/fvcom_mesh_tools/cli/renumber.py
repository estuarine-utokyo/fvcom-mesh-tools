"""``fmesh-renumber`` CLI: renumber a finished FVCOM case for locality.

It works on a written case -- ``<prefix>_grd.dat`` and friends -- rather than
inside the mesh generators, because the numbering a generator produces is
meaningful to it. A local refinement keeps its retained elements first and in
base order, which is what makes the frozen zone verifiable; renumbering
belongs after that contract has been checked, not instead of it.

The sweep starts at the open boundary by default, which is what SMS does --
select the open-boundary nodestring, then renumber -- and what the production
base was numbered by: seeded there, RCM reproduces its bandwidth (78 against
the mesh's own 80) where an automatically chosen start gives 128. It also
leaves the open-boundary nodes as a contiguous block at the end of the
numbering, which is the other reason to do it that way.

Every file of the case moves together: the grid, the depths, the Coriolis
column, the open-boundary list and a sponge file when there is one. The
permutation is written beside them as ``<prefix>_renumber.json`` so that
anything else keyed by node id -- a node map, an observation index, a station
list -- can be carried across afterwards.

    fmesh-renumber outputs/banzu_nori_s5/fvcom/banzu_nori --outdir <dir>

Why it matters for FVCOM, and where it does not: FVCOM partitions the domain
itself at run time (``setup_domain.F`` -> ``DOMDEC`` -> METIS), so this does
not change the partition. It changes the order inside each rank, because
``genmap.F`` builds local numbering by walking the global ids in order and
keeping the ones that belong to the rank.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.renumber import METHODS, locality_stats, renumber_mesh


def _read_sponge(path: Path) -> list[tuple[int, float, float]]:
    """``<prefix>_spg.dat`` as ``(node0, radius, coefficient)`` triples."""
    rows = [ln.split() for ln in path.read_text().splitlines()[1:] if ln.strip()]
    return [(int(r[0]) - 1, float(r[1]), float(r[2])) for r in rows]


def _read_cor(path: Path) -> np.ndarray:
    rows = [ln.split() for ln in path.read_text().splitlines()[1:] if ln.strip()]
    return np.array([float(r[-1]) for r in rows], dtype=float)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-renumber",
        description="Renumber a finished FVCOM case for memory locality.")
    p.add_argument("case", type=Path,
                   help="case prefix, without _grd.dat")
    p.add_argument("--outdir", type=Path, default=None,
                   help="where to write the renumbered case "
                        "(default: <case>_renumbered beside the input)")
    p.add_argument("--casename", default=None,
                   help="output case name (default: the input's)")
    p.add_argument("--method", choices=METHODS, default="rcm",
                   help="node ordering (default: rcm)")
    p.add_argument("--seed", choices=("obc", "auto", "boundary"), default="obc",
                   help="where the RCM sweep starts (default: obc, which is "
                        "SMS's convention and the production base's)")
    p.add_argument("--elements", choices=("node", "centroid"), default="node",
                   help="element ordering (default: follow the nodes)")
    p.add_argument("--force", action="store_true",
                   help="renumber even when it would not improve the mesh")
    p.add_argument("--dry-run", action="store_true",
                   help="measure and report; write nothing")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    prefix = args.case
    grd, dep = Path(f"{prefix}_grd.dat"), Path(f"{prefix}_dep.dat")
    obc = Path(f"{prefix}_obc.dat")
    for q in (grd, dep):
        if not q.exists():
            print(f"not found: {q}", file=sys.stderr)
            return 2
    mesh = read_fvcom_case(grd, dep, obc if obc.exists() else None)
    before = locality_stats(mesh.n_nodes, mesh.elements)
    out_mesh, report = renumber_mesh(mesh, args.method, elements=args.elements,
                                     seeds=args.seed,
                                     only_if_better=not args.force)
    after = report["after"]
    print(f"[renumber] {args.method}(seed={args.seed})/{args.elements}: bandwidth "
          f"{before['bandwidth']:,} -> {after['bandwidth']:,}, mean edge span "
          f"{before['mean_edge_span']:.1f} -> {after['mean_edge_span']:.1f}, "
          f"p99 {before['p99_edge_span']:.1f} -> {after['p99_edge_span']:.1f}")
    if not report["applied"]:
        print(f"[renumber] not applied: {report['note']}")
        if not args.dry_run:
            return 1
    if args.dry_run:
        return 0

    outdir = args.outdir or prefix.parent / f"{prefix.name}_renumbered"
    outdir.mkdir(parents=True, exist_ok=True)
    name = args.casename or prefix.name
    perm = np.asarray(report["node_perm"], dtype=np.int64)

    cor_path = Path(f"{prefix}_cor.dat")
    cor = _read_cor(cor_path)[perm] if cor_path.exists() else None
    spg_path = Path(f"{prefix}_spg.dat")
    sponge = None
    if spg_path.exists():
        inv = np.empty(perm.shape[0], dtype=np.int64)
        inv[perm] = np.arange(perm.shape[0])
        sponge = [(int(inv[n]), r, c) for n, r, c in _read_sponge(spg_path)]

    written = export_fvcom_case(
        out_mesh, outdir, name, cor=cor, sponge=sponge,
        obc_type=getattr(mesh, "obc_type", 1), obc_depth_control=False)
    # The permutation is the only way back: anything keyed by node id has to
    # be carried across with it, and a renumbered case with no record of how
    # it was renumbered is a case nobody can compare with its own history.
    (outdir / f"{name}_renumber.json").write_text(json.dumps({
        "method": args.method, "seed": args.seed,
        "element_order": args.elements,
        "source": str(prefix.resolve()),
        "before": before, "after": after,
        "node_perm_old_ids": perm.tolist(),
        "element_perm_old_ids": np.asarray(
            report["element_perm"], dtype=np.int64).tolist(),
    }) + "\n")
    print(f"[renumber] wrote {', '.join(sorted(written))} in {outdir}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
