"""Base-mesh recipes: everything a whole-bay mesh depends on, in one file.

A base recipe (``recipes/base/*.yaml``) names the open boundary (an input
file of lon,lat nodes), the domain closure, the land window, the geometry
edits, every generation and finishing setting, and the depth product.
``notebooks/440_base_mesh.py`` builds the mesh from it and records what it
used; this module reads and checks the recipe.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = ["DEPTH_PRODUCTS", "REFERENCE_KEYS", "REQUIRED_SETTINGS", "compare_to_reference",
           "load_base_recipe", "read_open_boundary"]

#: Depth products the build knows.
DEPTH_PRODUCTS = ("m7001_production",)

#: Settings a recipe must write out: a default changed in the code must not
#: change a base mesh silently.
REQUIRED_SETTINGS = (
    "SR_H0", "SR_MAXEL", "SR_GRADE", "SR_DT", "SR_CRMIN", "SR_H_TARGET",
    "SR_DM_SCALE", "SR_MODE", "SR_ZBASE", "SR_ZW_LAT", "SR_ZE_LAT", "SR_DMAX",
    "SR_CFL_DT", "SR_CFL_CRMAX", "SR_CH_POLICY", "SR_CH_REFINE", "SR_WATERWAYS",
    "SR_NORMALIZE", "SR_OBC_LADDER", "SR_OBC_HSCALE", "SR_OBC_K", "SR_OBC_SKIP",
    "SR_ONE_WIDE", "SR_COAST_FIT", "SR_FIN_SEED", "SR_GEN_SEED",
)
_OPTIONAL_SETTINGS = ("SR_OBC_H0", "SR_OBC_H1", "SR_EDITS_EXCLUDE")
#: What a recipe's ``reference`` block records of the mesh it reproduces.
#: Byte identity is not required: the environment follows the latest
#: conda-forge releases and the oceanmesh fork evolves upward-compatibly
#: (owner, 2026-09-29), so a mesh that differs slightly but keeps these
#: within ``tolerance`` and passes every QA gate reproduces the reference.
REFERENCE_KEYS = ("fort14_sha256", "n_nodes", "n_elements", "n_obc_nodes", "wet_area_km2",
                  "tolerance")
_TOLERANCE_KEYS = ("n_nodes", "n_elements", "wet_area_km2")
_PATH_SETTINGS = ("SR_OUT", "SR_LAND", "SR_OBC_FILE", "SR_DOMAIN_FILE", "SR_EDITS_DIR",
                  "SR_SIZING")


def read_open_boundary(path) -> list[tuple[float, float]]:
    """The open boundary's nodes, ``[(lon, lat), ...]``, in order.

    One ``lon,lat`` pair per line; ``#`` lines and a ``lon,lat`` header are
    skipped.  At least two nodes, all finite.
    """
    import math

    pts = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        s = line.strip()
        if not s or s.startswith("#") or s.lower().startswith("lon"):
            continue
        parts = [v.strip() for v in s.split(",")]
        try:
            lon, lat = float(parts[0]), float(parts[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}:{n}: expected 'lon,lat', got {line!r}") from exc
        if not (math.isfinite(lon) and math.isfinite(lat)
                and -180 <= lon <= 360 and -90 <= lat <= 90):
            raise ValueError(f"{path}:{n}: not a lon,lat in degrees: {line!r}")
        pts.append((lon, lat))
    if len(pts) < 2:
        raise ValueError(f"{path}: an open boundary needs two nodes or more")
    return pts


def load_base_recipe(path) -> dict[str, Any]:
    """Read and check a base recipe; relative paths resolve against it.

    Returns the recipe with ``open_boundary``, ``domain`` and ``edits``
    resolved to absolute paths, ``settings`` as strings, and ``recipe_path``.
    """
    import json

    import yaml

    path = Path(path).resolve()
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: a base recipe is a mapping")
    missing = [k for k in ("name", "case", "open_boundary", "domain", "land", "edits",
                           "settings", "depths") if k not in raw]
    if missing:
        raise ValueError(f"{path}: missing {missing}")
    out = dict(raw)
    out["recipe_path"] = str(path)
    for key in ("open_boundary", "domain", "edits"):
        p = (path.parent / str(raw[key])).resolve()
        if not p.exists():
            raise ValueError(f"{path}: {key} {p} does not exist")
        out[key] = str(p)
    read_open_boundary(out["open_boundary"])
    dom = json.loads(Path(out["domain"]).read_text())
    for k in ("bbox", "closure_before", "closure_after", "closure_size_probe"):
        if k not in dom:
            raise ValueError(f"{out['domain']}: missing {k!r}")
    land = raw["land"]
    if not (isinstance(land, dict) and len(land.get("bbox", [])) == 4):
        raise ValueError(f"{path}: land.bbox must be [lon_min, lat_min, lon_max, lat_max]")
    settings = raw["settings"]
    if not isinstance(settings, dict):
        raise ValueError(f"{path}: settings is a mapping")
    unknown = sorted(set(settings) - set(REQUIRED_SETTINGS) - set(_OPTIONAL_SETTINGS))
    if unknown:
        # a path setting belongs to the run, and a misspelt knob would be
        # silently ignored by the generator
        raise ValueError(f"{path}: unknown settings {unknown}"
                         + (" (paths are set by the build)"
                            if set(unknown) & set(_PATH_SETTINGS) else ""))
    absent = [k for k in REQUIRED_SETTINGS if k not in settings]
    if absent:
        raise ValueError(f"{path}: settings must write out {absent}")
    out["settings"] = {k: str(v) for k, v in settings.items()}
    if raw["depths"] not in DEPTH_PRODUCTS:
        raise ValueError(f"{path}: depths must be one of {DEPTH_PRODUCTS}")
    ref = raw.get("reference")
    if ref is not None:
        if not isinstance(ref, dict) or sorted(ref) != sorted(REFERENCE_KEYS):
            raise ValueError(f"{path}: reference needs exactly {list(REFERENCE_KEYS)}")
        sha = ref["fort14_sha256"]
        if not (isinstance(sha, str) and len(sha) == 64):
            raise ValueError(f"{path}: reference.fort14_sha256 is a SHA-256 hex digest")
        tol = ref["tolerance"]
        if not (isinstance(tol, dict) and sorted(tol) == sorted(_TOLERANCE_KEYS)
                and all(isinstance(v, (int, float)) and 0 <= v < 1 for v in tol.values())):
            raise ValueError(f"{path}: reference.tolerance gives a relative bound in [0, 1) "
                             f"for each of {list(_TOLERANCE_KEYS)}")
        for k in ("n_nodes", "n_elements", "n_obc_nodes", "wet_area_km2"):
            if not (isinstance(ref[k], (int, float)) and ref[k] > 0):
                raise ValueError(f"{path}: reference.{k} must be positive")
    return out


def compare_to_reference(summary: dict[str, Any], reference: dict[str, Any],
                         qa_passed: bool) -> dict[str, Any]:
    """Does a built mesh reproduce the recipe's reference?

    ``summary`` holds the built mesh's ``fort14_sha256``, ``n_nodes``,
    ``n_elements``, ``n_obc_nodes`` and ``wet_area_km2``.  ``byte_identical``
    compares the hashes; ``reproduces`` asks for every QA gate passed, the
    same open-boundary node count (the boundary is an input), and node count,
    element count and wet area each within its relative tolerance.
    """
    tol = reference["tolerance"]
    rel = {k: (summary[k] - reference[k]) / reference[k] for k in _TOLERANCE_KEYS}
    within = {k: abs(rel[k]) <= tol[k] for k in _TOLERANCE_KEYS}
    same_obc = summary["n_obc_nodes"] == reference["n_obc_nodes"]
    return {
        "byte_identical": summary["fort14_sha256"] == reference["fort14_sha256"],
        "reproduces": bool(qa_passed and same_obc and all(within.values())),
        "qa_passed": bool(qa_passed),
        "same_open_boundary_nodes": bool(same_obc),
        "relative_difference": rel,
        "within_tolerance": within,
        "built": {k: summary[k] for k in ("n_nodes", "n_elements", "n_obc_nodes",
                                          "wet_area_km2")},
    }
