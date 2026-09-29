"""Extension recipes: a base mesh kept as it is, plus the sea out to a new boundary.

An extension recipe (``recipes/extend/*.yaml``) names the base case, the new
open boundary (a CSV of lon,lat nodes), the land window, the bathymetry
sources for sizing and for depths, every sizing setting and the depth rules.
``notebooks/445_extend_mesh.py`` builds the mesh from it; this module reads
and checks the recipe.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fvcom_mesh_tools.base_recipe import read_open_boundary

__all__ = ["REQUIRED_SETTINGS", "load_extend_recipe"]

#: Every sizing and generation setting must be written out.
REQUIRED_SETTINGS = ("coast_h_m", "max_edge_m", "gradation", "cfl_dt_s", "cfl_cr",
                     "interface_band_m", "obc_band_m", "lattice_m", "dm_scale", "gen_seed",
                     "max_iter", "fin_seed")
_KEYS = ("name", "case", "base", "base_case", "open_boundary", "land", "bathymetry",
         "settings", "depths")


def load_extend_recipe(path) -> dict[str, Any]:
    """Read and check an extension recipe; relative paths resolve against it."""
    import yaml

    from fvcom_mesh_tools.dem.sources import SOURCES

    path = Path(path).resolve()
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: an extension recipe is a mapping")
    missing = [k for k in _KEYS if k not in raw]
    unknown = sorted(set(raw) - set(_KEYS))
    if missing or unknown:
        raise ValueError(f"{path}: missing {missing}, unknown {unknown}")
    out = dict(raw)
    out["recipe_path"] = str(path)
    for key in ("base", "open_boundary"):
        p = (path.parent / str(raw[key])).resolve()
        if not p.exists():
            raise ValueError(f"{path}: {key} {p} does not exist")
        out[key] = str(p)
    for kind in ("grd", "dep", "obc"):
        f = Path(out["base"]) / f"{raw['base_case']}_{kind}.dat"
        if not f.exists():
            raise ValueError(f"{path}: the base case has no {f.name}")
    read_open_boundary(out["open_boundary"])
    land = raw["land"]
    if not (isinstance(land, dict) and set(land) == {"bbox"} and len(land["bbox"]) == 4):
        raise ValueError(f"{path}: land is {{bbox: [lon_min, lat_min, lon_max, lat_max]}}")
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
        if k not in ("gen_seed", "fin_seed") and s[k] <= 0:
            raise ValueError(f"{path}: settings.{k} must be positive")
    if s["coast_h_m"] > s["max_edge_m"]:
        raise ValueError(f"{path}: coast_h_m exceeds max_edge_m")
    d = raw["depths"]
    if not (isinstance(d, dict) and set(d) == {"min_m", "max_m", "rfactor"}):
        raise ValueError(f"{path}: depths has exactly min_m, max_m, rfactor")
    if not d["min_m"] > 0 or (d["max_m"] is not None and not d["max_m"] > d["min_m"]):
        raise ValueError(f"{path}: depths need 0 < min_m < max_m (max_m may be null)")
    if not 0 < d["rfactor"] < 1:
        raise ValueError(f"{path}: depths.rfactor must be in (0, 1)")
    return out
