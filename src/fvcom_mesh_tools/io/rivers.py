"""Load river-mouth point coordinates from a vector or text file.

Accepted formats:

* GeoJSON / shapefile / any vector file readable by ``geopandas`` whose
  geometries are points (or multipoints; flattened).
* Plain CSV with at least one of the column-name pairs
  ``(lon, lat)`` / ``(longitude, latitude)`` / ``(x, y)``. Column
  matching is case-insensitive.

All inputs are reprojected to EPSG:4326 (lon/lat) on load.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _from_csv(path: Path) -> np.ndarray:
    pts: list[tuple[float, float]] = []
    with path.open() as f:
        header_line = f.readline().rstrip("\n")
        cols = [c.strip().lower() for c in header_line.split(",")]
        # Find lon/lat column positions.
        lon_keys = ("lon", "longitude", "x")
        lat_keys = ("lat", "latitude", "y")
        lon_idx = next((i for i, c in enumerate(cols) if c in lon_keys), None)
        lat_idx = next((i for i, c in enumerate(cols) if c in lat_keys), None)
        if lon_idx is None or lat_idx is None:
            raise ValueError(
                f"{path}: header must contain a lon/lat column pair "
                f"(got columns {cols!r})."
            )
        for line in f:
            line = line.rstrip("\n").strip()
            if not line:
                continue
            row = [c.strip() for c in line.split(",")]
            pts.append((float(row[lon_idx]), float(row[lat_idx])))
    return np.asarray(pts, dtype=np.float64)


def _from_vector(path: Path) -> np.ndarray:
    try:
        import geopandas as gpd
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "loading vector river files requires geopandas; "
            "install via `mamba install -c conda-forge geopandas`."
        ) from exc

    gdf = gpd.read_file(path)
    if gdf.empty:
        return np.empty((0, 2), dtype=np.float64)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    else:
        gdf = gdf.to_crs("EPSG:4326")

    pts: list[tuple[float, float]] = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        gt = geom.geom_type
        if gt == "Point":
            pts.append((geom.x, geom.y))
        elif gt == "MultiPoint":
            for sub in geom.geoms:
                pts.append((sub.x, sub.y))
        else:
            raise ValueError(f"{path}: unsupported geometry type {gt}")
    return np.asarray(pts, dtype=np.float64)


def load_river_points(paths: Iterable[Path | str]) -> np.ndarray:
    """Read river-mouth points from one or more files into ``(N, 2)``
    ``(lon, lat)`` array. Points are concatenated in the order of
    ``paths``.

    Raises ``FileNotFoundError`` if any path does not exist, and
    ``ValueError`` for file formats that cannot be interpreted.
    """
    chunks: list[np.ndarray] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(f"river points file not found: {path}")
        suffix = path.suffix.lower()
        if suffix == ".csv":
            chunks.append(_from_csv(path))
        else:
            chunks.append(_from_vector(path))
    if not chunks:
        return np.empty((0, 2), dtype=np.float64)
    return np.vstack(chunks)
