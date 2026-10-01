"""Extension recipes: a base mesh kept as it is, plus the sea out to a new boundary.

An extension recipe (``recipes/extend/*.yaml``) names the base case, the new
open boundary (a CSV of lon,lat nodes), the land window, the bathymetry
sources for sizing and for depths, every sizing setting and the depth rules.
``notebooks/445_extend_mesh.py`` builds the mesh from it; this module reads
and checks the recipe.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Any

from fvcom_mesh_tools.base_recipe import parse_open_boundary

__all__ = ["EXPECT_ENV", "REQUIRED_SETTINGS", "check_case_name", "check_expected",
           "load_extend_recipe"]

#: Every sizing and generation setting must be written out.
REQUIRED_SETTINGS = ("coast_h_m", "max_edge_m", "gradation", "cfl_dt_s", "cfl_cr",
                     "interface_band_m", "obc_band_m", "lattice_m", "dm_scale", "gen_seed",
                     "max_iter", "fin_seed")
_KEYS = ("name", "case", "base", "base_case", "open_boundary", "land", "bathymetry",
         "settings", "depths")
#: Settings that are counts or seeds: integral, not truncated (review F18).
_INTEGRAL = ("gen_seed", "fin_seed", "max_iter")


def check_case_name(name) -> str:
    """A case name is one plain file-name component (review F28).

    It prefixes every output file (``<case>_grd.dat``); a path in it could
    leave the reserved output directory.
    """
    if not (isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name)
            and name not in (".", "..")):
        raise ValueError(f"case name {name!r} must be one file-name component "
                         "(letters, digits, '_', '-', '.')")
    return name


def load_extend_recipe(path) -> dict[str, Any]:
    """Read and check an extension recipe; relative paths resolve against it."""
    import yaml

    from fvcom_mesh_tools.dem.sources import SOURCES

    path = Path(path).resolve()
    # the bytes parsed are the bytes hashed: a provenance hash taken later
    # could describe another file (review of the extend tools, round 5 F3)
    data = path.read_bytes()
    raw = yaml.safe_load(data.decode())
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: an extension recipe is a mapping")
    missing = [k for k in _KEYS if k not in raw]
    unknown = sorted(set(raw) - set(_KEYS))
    if missing or unknown:
        raise ValueError(f"{path}: missing {missing}, unknown {unknown}")
    out = dict(raw)
    out["recipe_path"] = str(path)
    out["recipe_sha256"] = hashlib.sha256(data).hexdigest()
    try:
        check_case_name(raw["case"])
        check_case_name(raw["base_case"])
    except ValueError as err:
        raise ValueError(f"{path}: {err}") from None
    for key in ("base", "open_boundary"):
        p = (path.parent / str(raw[key])).resolve()
        if not p.exists():
            raise ValueError(f"{path}: {key} {p} does not exist")
        out[key] = str(p)
    for kind in ("grd", "dep", "obc"):
        f = Path(out["base"]) / f"{raw['base_case']}_{kind}.dat"
        if not f.exists():
            raise ValueError(f"{path}: the base case has no {f.name}")
    # the boundary is read once: the bytes hashed are the bytes parsed, and
    # the stages use these coordinates rather than reading the file again
    # (review round 7 F8)
    obc_bytes = Path(out["open_boundary"]).read_bytes()
    out["open_boundary_lonlat"] = parse_open_boundary(obc_bytes.decode(), out["open_boundary"])
    # a boundary published by 444 carries its report beside it, with the CSV's
    # hash: a CSV that is not the one the report describes is refused (review
    # round 3 F5)
    marker = Path(out["open_boundary"] + ".RECOVER")
    if marker.exists():
        raise ValueError(f"{path}: a failed publication left {marker.name}; restore the "
                         f"boundary from its .prev files first")
    side = Path(out["open_boundary"]).with_suffix(".json")
    got = hashlib.sha256(obc_bytes).hexdigest()
    out["open_boundary_sha256"] = got
    if side.exists():
        import json

        want = json.loads(side.read_text()).get("csv_sha256")
        if want is not None and want != got:
            raise ValueError(f"{path}: {Path(out['open_boundary']).name} is not the boundary "
                             f"its report {side.name} describes (hash mismatch)")
    land = raw["land"]
    if not (isinstance(land, dict) and set(land) == {"bbox"} and len(land["bbox"]) == 4):
        raise ValueError(f"{path}: land is {{bbox: [lon_min, lat_min, lon_max, lat_max]}}")
    bb = land["bbox"]
    # finite, on the globe, and ordered (review round 2 F11)
    if not (all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
                for v in bb)
            and -180 <= bb[0] < bb[2] <= 360 and -90 <= bb[1] < bb[3] <= 90):
        raise ValueError(f"{path}: land.bbox {bb} must be finite lon_min < lon_max, "
                         "lat_min < lat_max on the globe")
    bathy = raw["bathymetry"]
    if not (isinstance(bathy, dict) and set(bathy) == {"sizing", "depths"}):
        raise ValueError(f"{path}: bathymetry has exactly 'sizing' and 'depths'")
    for use, names in bathy.items():
        if not (isinstance(names, list) and names):
            raise ValueError(f"{path}: bathymetry.{use} is a non-empty list")
        bad = [n for n in names if n not in SOURCES]
        if bad or len(set(names)) != len(names):
            raise ValueError(f"{path}: bathymetry.{use}: unknown or repeated {bad or names}; "
                             f"known {sorted(SOURCES)}")
    s = raw["settings"]
    if not isinstance(s, dict):
        raise ValueError(f"{path}: settings is a mapping")
    absent = [k for k in REQUIRED_SETTINGS if k not in s]
    extra = sorted(set(s) - set(REQUIRED_SETTINGS))
    if absent or extra:
        raise ValueError(f"{path}: settings must write out {absent}; unknown {extra}")
    for k in REQUIRED_SETTINGS:
        if not (isinstance(s[k], (int, float)) and not isinstance(s[k], bool)):
            raise ValueError(f"{path}: settings.{k} must be a number")
        if not math.isfinite(s[k]):
            raise ValueError(f"{path}: settings.{k} must be finite")
        if k in _INTEGRAL and s[k] != int(s[k]):
            raise ValueError(f"{path}: settings.{k} must be an integer")
        if k in ("gen_seed", "fin_seed") and not 0 <= s[k] < 2**32:
            raise ValueError(f"{path}: settings.{k} must be in [0, 2**32)")
        if k not in ("gen_seed", "fin_seed") and s[k] <= 0:
            raise ValueError(f"{path}: settings.{k} must be positive")
    if s["coast_h_m"] > s["max_edge_m"]:
        raise ValueError(f"{path}: coast_h_m exceeds max_edge_m")
    d = raw["depths"]
    if not (isinstance(d, dict) and set(d) == {"min_m", "max_m", "rfactor"}):
        raise ValueError(f"{path}: depths has exactly min_m, max_m, rfactor")
    for k in ("min_m", "max_m", "rfactor"):
        v = d[k]
        if k == "max_m" and v is None:
            continue
        if not (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)):
            raise ValueError(f"{path}: depths.{k} must be a finite number")
    if not d["min_m"] > 0 or (d["max_m"] is not None and not d["max_m"] > d["min_m"]):
        raise ValueError(f"{path}: depths need 0 < min_m < max_m (max_m may be null)")
    if not 0 < d["rfactor"] < 1:
        raise ValueError(f"{path}: depths.rfactor must be in (0, 1)")
    return out


#: environment variables a driver sets for its stages: the digests of the
#: recipe and open boundary it recorded (review round 5 F3)
EXPECT_ENV = {"recipe_sha256": "FMESH_EXPECT_RECIPE_SHA256",
              "open_boundary_sha256": "FMESH_EXPECT_OBC_SHA256"}


def check_expected(recipe: dict[str, Any]) -> None:
    """Refuse a recipe or open boundary that is not the one the driver
    recorded, when the driver says which (``EXPECT_ENV``)."""
    import os

    for key, var in EXPECT_ENV.items():
        want = os.environ.get(var)
        if want and want != recipe[key]:
            raise ValueError(f"{key.removesuffix('_sha256')} changed since the driver "
                             f"recorded it ({recipe[key][:12]} != {want[:12]})")

