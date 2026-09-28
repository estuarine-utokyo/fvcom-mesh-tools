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

__all__ = ["DEPTH_PRODUCTS", "REQUIRED_SETTINGS", "load_base_recipe", "read_open_boundary"]

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
    ref = raw.get("reference_fort14_sha256")
    if ref is not None and not (isinstance(ref, str) and len(ref) == 64):
        raise ValueError(f"{path}: reference_fort14_sha256 is a SHA-256 hex digest")
    return out
