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

__all__ = ["changed_files", "changed_inventory", "code_state", "collect", "dataset_files",
           "file_sha256", "git_state"]

# the libraries whose arithmetic or geometry decides the mesh
LIBRARIES = ("numpy", "scipy", "shapely", "geopandas", "pyproj", "rasterio",
             "matplotlib", "netCDF4", "oceanmesh", "fvcom-mesh-tools")

# the suffix of the marker dataset_files adds when sidecars cannot be listed
UNLISTED = ".<sidecars-unlisted>"


def file_sha256(path) -> str | None:
    """SHA-256 of a file, or None when it cannot be read."""
    if str(path).endswith(UNLISTED):
        return None                     # a dataset whose sidecars are unknown
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()
    except OSError:
        return None


def dataset_files(path) -> list[Path]:
    """The files a dataset is read from: a shapefile's sidecars too.

    A shapefile's geometry is in ``.shp`` and ``.shx``, its attributes (which
    a ``where:`` selects on) in ``.dbf``, its CRS in ``.prj`` and its text
    encoding in ``.cpg``; any single-file format (GeoJSON, GeoPackage) is
    itself (review, round 2).

    When the directory cannot be listed, which sidecars exist is unknown --
    not the same as "none" -- and the list ends with a marker path,
    ``<stem>.<sidecars-unlisted>``, that names this and hashes to None, so
    the dataset's digest in :func:`collect` is null rather than a digest of
    the ``.shp`` alone (review, round 8).
    """
    p = Path(path)
    if p.suffix.lower() != ".shp":
        return [p]
    # matched without case: GDAL reads LAND.SHP with LAND.PRJ (review, round 3)
    side = {".shx", ".dbf", ".prj", ".cpg"}
    try:
        sibs = sorted(q for q in p.parent.iterdir()
                      if q.stem == p.stem and q.suffix.lower() in side)
    except OSError:
        return [p, p.with_name(p.stem + UNLISTED)]
    return [p] + sibs



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
        # commit alone cannot give back (review, round 1).  Untracked files
        # are asked for explicitly: status.showUntrackedFiles=no hid them
        # (review, round 16).
        changed = parse_porcelain_z(run("status", "--porcelain", "-z",
                                        "--untracked-files=all").stdout)
        tracked = None
        if target.is_file():
            tracked = run("ls-files", "--error-unmatch", str(target),
                          check=False).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return None
    return {"root": top, "commit": commit, "dirty": changed, "path_tracked": tracked}


def parse_porcelain_z(out: str) -> list[str]:
    """Every path ``git status --porcelain -z`` names, both sides of a rename.

    NUL-separated records ``XY path``; a rename or copy (X or Y in ``R``/``C``)
    is followed by one more record, the source path. Line-based parsing kept
    ``old -> new`` as one path, so a file renamed INTO a tree did not mark it
    dirty, and names with spaces or quotes were mangled (review of the extend
    tools, round 2 F10).
    """
    recs = out.split("\0")
    paths, i = [], 0
    while i < len(recs):
        rec = recs[i]
        i += 1
        if len(rec) < 4:
            continue
        paths.append(rec[3:])
        if "R" in rec[:2] or "C" in rec[:2]:
            if i < len(recs) and recs[i]:
                paths.append(recs[i])
            i += 1
    return paths


def code_state(name: str, path) -> dict[str, Any]:
    """What identifies the code at ``path``, in or out of git.

    Always the resolved path and the installed distribution's version (if
    ``name`` is one); the git state when the path is in a work tree; and a
    SHA-256 over the Python sources beside it, sorted by name, whenever the
    commit does not identify them -- outside git (a package installed there
    had no identity at all, review round 4) and for a path the enclosing
    repository does not track, such as an ignored ``build/`` or ``.venv``
    inside a checkout (round 21).  ``commit_identifies_code`` says which.
    """
    p = Path(path).resolve()
    out: dict[str, Any] = {"path": str(p), "distribution": _version(name),
                           "git": git_state(p)}
    out["commit_identifies_code"] = bool(out["git"] and out["git"].get("path_tracked"))
    if out["commit_identifies_code"]:
        # a tracked tree with uncommitted changes under it is not what the
        # commit says (review of the extend tools, round 1, F17)
        top = Path(out["git"]["root"])
        here = p.parent if p.is_file() else p
        dirty_here = [d for d in out["git"].get("dirty") or []
                      if (top / d).resolve() == here or here in (top / d).resolve().parents]
        if dirty_here:
            out["commit_identifies_code"] = False
            out["dirty_under_path"] = dirty_here
    if not out["commit_identifies_code"]:
        root = p.parent if p.is_file() else p
        # each file is read once, and one that cannot be read makes the
        # identity incomplete, said so, rather than a digest that looks whole
        # (review, round 5)
        # Path.rglob swallows a directory it cannot scan, and an unreadable
        # tree hashed as empty (review, round 6): walk with the errors kept
        # Symlinked directories are followed -- Python imports through them,
        # and a subpackage behind one went unhashed (review, round 13) -- in
        # sorted order and under EVERY logical path, since each alias is an
        # importable path; a link to one of its own ancestors is a cycle and
        # ends the walk there (review, round 14: first-alias-wins made the
        # digest depend on traversal order and ignore a removed alias).
        h, failed = hashlib.sha256(), []
        files, links = [], []
        for d, dirs, names in root.walk(follow_symlinks=True, on_error=lambda e: failed.append(
                str(e.filename) if e.filename else str(root))):
            try:
                real = d.resolve()
                ancestors = {(root / a).resolve() for a in d.relative_to(root).parents}
            except OSError:
                failed.append(str(d))
                dirs[:] = []
                continue
            if real in ancestors:
                # a link back to an ancestor is an importable path too: it is
                # recorded, as the ancestor it names inside the tree, though
                # it is not walked again (review, round 15)
                up = next(str(a) for a in d.relative_to(root).parents
                          if (root / a).resolve() == real)
                links.append(f"{d.relative_to(root)} -> {up}")
                dirs[:] = []
                continue
            dirs.sort()
            files.extend(d / n for n in names if n.endswith(".py"))
        for f in sorted(files):
            digest = file_sha256(f)
            if digest is None:
                failed.append(str(f))
                continue
            h.update(str(f.relative_to(root)).encode())
            h.update(digest.encode())
        for ln in sorted(links):
            h.update(ln.encode())
        out["source_sha256"] = None if failed else h.hexdigest()
        out["source_unreadable"] = failed
    return out


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
        "code": {name: code_state(name, p) for name, p in (code or {}).items()},
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


def changed_files(prov: dict[str, Any], names) -> list[str]:
    """The entries ``names`` of ``prov["files"]`` whose content is not what
    :func:`collect` recorded: changed, removed or created since.

    A build reads its inputs after its provenance is taken; checking them
    again at the end binds the record to what was consumed (review of the
    extend tools, round 4 F8).
    """
    spec = {}
    for k in names:
        v = prov["files"][k]
        spec[k] = v["paths"] if "paths" in v else v["path"]
    now = collect(files=spec, libraries=())["files"]
    return [k for k in spec if now[k]["sha256"] != prov["files"][k]["sha256"]]


def changed_inventory(prov: dict[str, Any], files: dict[str, Any]) -> list[str]:
    """The entries of ``files`` (name -> path or list of paths, listed again
    now) whose membership differs from what ``prov["files"]`` recorded, or
    that it did not record. A dataset whose files are found by a search can
    gain a file the first hash never saw (review of the extend tools, round 7
    F9); :func:`changed_files` then checks the contents."""
    out = []
    for k, p in files.items():
        rec = prov["files"].get(k)
        now = [str(q) for q in p] if isinstance(p, (list, tuple)) else str(p)
        was = None if rec is None else rec.get("paths", rec.get("path"))
        if now != was:
            out.append(k)
    return out

