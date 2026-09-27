"""What a mesh was made from, so that it can be made again.

A refinement is reproducible from its seed only if everything else is the
same too: the code of this package and of the oceanmesh fork, the input
files, and the libraries that do the floating-point work.  Rebuilding the
five port recipes gave byte-identical meshes (2026-09-25), but the report
recorded only the seed, the recipe and the base mesh's hash -- not the code
or the OSM coastline it was cut from.  :func:`collect` records the rest.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import subprocess
from pathlib import Path
from typing import Any

__all__ = ["collect", "file_sha256", "git_state"]

# the libraries whose arithmetic or geometry decides the mesh
LIBRARIES = ("numpy", "scipy", "shapely", "geopandas", "pyproj", "rasterio",
             "matplotlib", "netCDF4", "oceanmesh")


def file_sha256(path) -> str | None:
    """SHA-256 of a file, or None when it cannot be read."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()
    except OSError:
        return None


def git_state(path) -> dict[str, Any] | None:
    """The commit a source tree is at, and whether it has local changes.

    ``dirty`` lists changed and untracked files: a mesh made from
    uncommitted code cannot be remade from the commit alone, and saying so is
    the point.  ``path_tracked`` says whether ``path`` itself is in the
    commit (None for a directory).
    Returns None outside a git work tree or without git.
    """
    root = Path(path).resolve()
    if root.is_file():
        root = root.parent

    target = Path(path).resolve()

    def run(*args, check=True):
        return subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                              text=True, timeout=30, check=check)

    try:
        top = run("rev-parse", "--show-toplevel").stdout.strip()
        commit = run("rev-parse", "HEAD").stdout.strip()
        # Porcelain lines are "XY path"; the leading status column may be a
        # space, so the output is not stripped before slicing.  Untracked
        # files count: a new, uncommitted source file is exactly what the
        # commit alone cannot give back (review, round 1).
        changed = [ln[3:] for ln in run("status", "--porcelain").stdout.splitlines()
                   if ln.strip()]
        tracked = None
        if target.is_file():
            tracked = run("ls-files", "--error-unmatch", str(target),
                          check=False).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return None
    return {"root": top, "commit": commit, "dirty": changed, "path_tracked": tracked}


def _version(name: str) -> str | None:
    """The installed distribution's version, without importing it.

    Importing would load ``oceanmesh``, which is GPL and may not be imported
    from this Apache-2.0 package (CLAUDE.md; review, round 1).  For an
    editable install the metadata is from install time; the commit in
    ``code`` is what identifies such a tree.
    """
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def collect(*, code: dict[str, Any] | None = None,
            files: dict[str, Any] | None = None,
            libraries=LIBRARIES) -> dict[str, Any]:
    """A provenance record for ``report.json``.

    ``code`` maps a name to a path inside a git tree (``{"fvcom_mesh_tools":
    __file__, "oceanmesh": oceanmesh.__file__}``); ``files`` maps a name to a
    file, or to a list of files hashed together in order (a shapefile and its
    sidecars).  Missing files hash to None rather than failing the run.
    """
    out: dict[str, Any] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "libraries": {name: _version(name) for name in libraries},
        "code": {name: git_state(p) for name, p in (code or {}).items()},
        "files": {},
    }
    for name, p in (files or {}).items():
        if isinstance(p, (list, tuple)):
            parts = [file_sha256(q) for q in p]
            out["files"][name] = {
                "paths": [str(q) for q in p],
                "sha256": None if any(x is None for x in parts)
                else hashlib.sha256("".join(parts).encode()).hexdigest()}
        else:
            out["files"][name] = {"path": str(p), "sha256": file_sha256(p)}
    return out
