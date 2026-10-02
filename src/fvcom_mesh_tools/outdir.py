"""Reserve an output directory for one run, atomically.

Checking that a directory is empty and then creating it lets two runs in
(review of the extend tools, rounds 1 and 2). ``reserve`` creates the
directory if needed, refuses one that holds anything but its own marker,
and takes ownership by creating ``.reserved`` with ``O_EXCL``: of two runs
racing for the same directory exactly one gets it.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["MARKER", "TOKEN_ENV", "claim", "new_token", "reserve"]

MARKER = ".reserved"
#: the environment variable a driver passes its reservation token in
TOKEN_ENV = "FMESH_OUT_TOKEN"


def new_token() -> str:
    """A fresh reservation token."""
    import secrets

    return secrets.token_hex(16)


def reserve(path, token: str | None = None) -> Path:
    """Create or take the empty directory ``path`` for this run; return it
    resolved. ``token``, if given, is written into the marker, so that the
    run's stages can prove they belong to it (:func:`claim`)."""
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir()
    except FileExistsError:
        if not path.is_dir():
            raise SystemExit(f"{path} exists and is not a directory") from None
        if any(p.name != MARKER for p in path.iterdir()):
            raise SystemExit(f"{path} is not empty; give a fresh directory") from None
    try:
        fd = os.open(path / MARKER, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise SystemExit(f"{path} is reserved by another run") from None
    os.write(fd, f"pid {os.getpid()}\n".encode()
             + (f"token {token}\n".encode() if token else b""))
    os.close(fd)
    return path


def claim(path, owner=None) -> Path:
    """The output of one stage of a run: ``path``, resolved.

    Under a driver -- ``TOKEN_ENV`` set -- the reservation marker of
    ``owner`` (``path`` itself by default) must carry that token: the stage
    writes only into a directory its own run reserved (review of the extend
    tools, round 19 F2). Run alone, the stage reserves ``path`` afresh, so an
    earlier output, finished or not, is refused.
    """
    path = Path(path).resolve()
    token = os.environ.get(TOKEN_ENV)
    if not token:
        return reserve(path)
    marker = Path(owner).resolve() / MARKER if owner is not None else path / MARKER
    try:
        held = marker.read_text().splitlines()
    except OSError:
        raise SystemExit(f"{marker} is missing: {path} is not reserved for this run") from None
    if f"token {token}" not in held:
        raise SystemExit(f"{path} is reserved by another run")
    path.mkdir(parents=True, exist_ok=True)
    return path

