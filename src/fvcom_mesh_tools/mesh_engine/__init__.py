"""Mesh-generation engine adapters.

Two backends are exposed via the ``fmesh-buildmesh --engine`` flag:

* ``oceanmesh`` -- DistMesh via :mod:`oceanmesh` (Roberts et al.). The
  Python port of OceanMesh2D; produces near-equilateral meshes
  (alpha mean ~0.96, ``frac<20deg`` ~0.03 % on Tokyo Bay) but is
  slower (PoC #18: ~26 min on Tokyo Bay).
* ``ocsmesh`` -- OCSMesh + gmsh. Faster (~40 s on Tokyo Bay) but
  lower quality (alpha ~0.85, ``frac<20deg`` ~1.13 %); useful for
  draft / iteration.

Both engines emit ``(points, cells)`` in EPSG:4326 lon/lat. The
caller (``fmesh-buildmesh``) is responsible for depth interpolation,
boundary classification, river inflow injection, perpfix, and
``fort.14`` writing - all mesher-agnostic.

Lazy imports keep the heavyweight dependencies optional: importing
:mod:`fvcom_mesh_tools.mesh_engine` does *not* import either engine
module.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

ENGINES = ("oceanmesh", "ocsmesh")


def build(
    engine: str,
    *,
    dem_path: Path,
    coastline_paths: list[Path],
    bbox: tuple[float, float, float, float],
    hmin_m: float,
    hmax_m: float,
    log: Callable[[str], None] = print,
    **engine_kwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """Dispatch to the requested mesh engine. Returns ``(points, cells)``.

    ``points`` is ``(NP, 2)`` lon/lat (EPSG:4326). ``cells`` is
    ``(NE, 3)`` zero-based vertex indices, CCW-oriented when output
    looks "normal" (otherwise the caller will flip).
    """
    if engine == "oceanmesh":
        from fvcom_mesh_tools.mesh_engine.oceanmesh import build as om_build
        return om_build(
            dem_path=dem_path,
            coastline_paths=coastline_paths,
            bbox=bbox,
            hmin_m=hmin_m,
            hmax_m=hmax_m,
            log=log,
            **engine_kwargs,
        )
    if engine == "ocsmesh":
        from fvcom_mesh_tools.mesh_engine.ocsmesh import build as ocs_build
        return ocs_build(
            dem_path=dem_path,
            coastline_paths=coastline_paths,
            bbox=bbox,
            hmin_m=hmin_m,
            hmax_m=hmax_m,
            log=log,
            **engine_kwargs,
        )
    raise ValueError(f"unknown engine: {engine!r}; choose from {ENGINES}")


__all__ = ["ENGINES", "build"]
