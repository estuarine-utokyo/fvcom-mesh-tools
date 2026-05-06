"""``fmesh-mesh-combine`` CLI: stitch independent ``fort.14`` files.

Three strategies, picked with ``--strategy``:

* ``disjoint`` - simple concat with index renumbering. Boundaries are
  preserved verbatim. No OCSMesh required. Best for non-overlapping
  basins (Tokyo Bay + Osaka Bay -> one model).
* ``overlap``  - OCSMesh-backed: carve foreground meshes into a
  background and remesh the seam. For nested-resolution scenarios.
* ``neighbor`` - OCSMesh-backed: snap meshes whose edges coincide.

Boundaries are dropped on the OCSMesh paths (``MeshData`` has no
boundary structure); pass the result through ``fmesh-buildmesh``-style
post-processing or ``omesh14-edit-bdy`` to re-establish them.

Example::

    fmesh-mesh-combine outputs/19_tokyo_bay_oceanmesh.14 \\
        outputs/20_osaka_bay_oceanmesh.14 \\
        outputs/21_kanto_kansai.14 \\
        --strategy disjoint
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_compose import STRATEGIES, combine


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-combine",
        description=(
            "Combine two or more fort.14 meshes into a single fort.14. "
            "Strategy 'disjoint' concatenates non-overlapping inputs "
            "with boundaries preserved (pure numpy). 'overlap' and "
            "'neighbor' use ocsmesh.ops.combine_mesh."
        ),
    )
    p.add_argument(
        "inputs", type=Path, nargs="+",
        help="Two or more input fort.14 files. The last one is treated as the output path.",
    )
    p.add_argument(
        "--strategy", choices=list(STRATEGIES), default="disjoint",
        help="Combination strategy (default: disjoint).",
    )
    p.add_argument(
        "--crs", type=str, default="EPSG:4326",
        help="CRS for OCSMesh-backed strategies (default: EPSG:4326).",
    )
    p.add_argument(
        "--buffer-size", type=float, default=0.0075,
        help="[overlap] OCSMesh buffer_size in CRS units (default: 0.0075).",
    )
    p.add_argument(
        "--buffer-domain", type=float, default=0.002,
        help="[overlap] OCSMesh buffer_domain in CRS units (default: 0.002).",
    )
    p.add_argument(
        "--min-int-ang", type=int, default=30,
        help="[overlap] Minimum interior angle for the seam remesh (default: 30 deg).",
    )
    p.add_argument(
        "--adjacent-layers", type=int, default=0,
        help="[overlap] OCSMesh adjacent_layers (default: 0).",
    )
    p.add_argument(
        "--no-clip-final", action="store_true",
        help="[overlap] Disable the final clip-to-domain step.",
    )
    p.add_argument(
        "--title", type=str, default=None,
        help="Title for the output fort.14 (default: combine of input titles).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if len(args.inputs) < 3:
        print(
            "fmesh-mesh-combine needs at least 2 input fort.14 files plus an output path "
            f"(got {len(args.inputs)} positional arguments).",
            file=sys.stderr,
        )
        return 2

    *input_paths, output = args.inputs
    for p in input_paths:
        if not p.exists():
            print(f"input not found: {p}", file=sys.stderr)
            return 2
    output.parent.mkdir(parents=True, exist_ok=True)

    meshes = [read_fort14(p) for p in input_paths]
    print(f"[mesh-combine] strategy={args.strategy}  inputs={len(meshes)}")
    for path, m in zip(input_paths, meshes):
        print(
            f"[mesh-combine]   {path.name}: NP={m.n_nodes:,} NE={m.n_elements:,} "
            f"open={len(m.open_boundaries)} land={len(m.land_boundaries)}"
        )

    kwargs: dict = {"title": args.title}
    if args.strategy == "overlap":
        kwargs.update(
            crs=args.crs,
            adjacent_layers=args.adjacent_layers,
            buffer_size=args.buffer_size,
            buffer_domain=args.buffer_domain,
            min_int_ang=args.min_int_ang,
            clip_final=not args.no_clip_final,
        )
    elif args.strategy == "neighbor":
        kwargs.update(crs=args.crs)
    # disjoint takes only `title`.

    combined = combine(args.strategy, meshes, **{k: v for k, v in kwargs.items() if v is not None})
    print(
        f"[mesh-combine] output: NP={combined.n_nodes:,} "
        f"NE={combined.n_elements:,} "
        f"open={len(combined.open_boundaries)} "
        f"land={len(combined.land_boundaries)}"
    )

    write_fort14(combined, output)
    print(f"[mesh-combine] wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
