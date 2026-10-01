"""Writing values into an FVCOM run namelist, and moving a staged case.

FVCOM reads ``INPUT_DIR`` and ``OUTPUT_DIR`` into ``CHARACTER(LEN=80)``
(``mod_main.F``): a longer path is cut without a word and the run reads or
writes elsewhere. Strings are Fortran list-directed values: an apostrophe
inside one is doubled, and a path is kept to printable ASCII (review of the
extend tools, rounds 4-6).
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

__all__ = ["FVCOM_DIR_MAX", "check_fvcom_dirs", "fortran_string", "relocate_case",
           "set_value"]

FVCOM_DIR_MAX = 80


def fortran_string(text: str) -> str:
    """``text`` as a quoted Fortran string; printable ASCII only."""
    if not (text.isascii() and text.isprintable()):
        raise ValueError(f"not printable ASCII, refused in a namelist: {text!r}")
    return "'" + text.replace("'", "''") + "'"


def check_fvcom_dirs(*dirs) -> list[str]:
    """The run directories as FVCOM will read them (absolute, with a trailing
    slash); ValueError for one that is not printable ASCII or is longer than
    ``FVCOM_DIR_MAX`` bytes."""
    paths = [f"{Path(d).resolve()}/" for d in dirs]
    bad = [p for p in paths if not (p.isascii() and p.isprintable())]
    if bad:
        raise ValueError(f"FVCOM run directories must be printable ASCII: {bad}")
    long = [p for p in paths if len(p.encode("ascii")) > FVCOM_DIR_MAX]
    if long:
        raise ValueError(f"FVCOM keeps {FVCOM_DIR_MAX} bytes of a run directory; "
                         f"use a shorter root: {long}")
    return paths


_VALUE = re.compile(r"""'(?:[^']|'')*'|"(?:[^"]|"")*"|[^,\s/]+""")


def _mask_quoted(text: str) -> str:
    """``text`` with quoted contents blanked (same length), so a search for
    ``KEY =`` cannot land inside a string."""
    out, quote, i = list(text), None, 0
    while i < len(text):
        ch = text[i]
        if quote:
            if ch == quote and i + 1 < len(text) and text[i + 1] == quote:
                out[i] = out[i + 1] = " "
                i += 2
                continue
            if ch == quote:
                quote = None
            elif ch != "\n":
                out[i] = " "
        elif ch in "'\"":
            quote = ch
        i += 1
    return "".join(out)


def set_value(text: str, key: str, value: str) -> str:
    """Replace the value of the one assignment ``key = ...`` with ``value``.

    Only that value is replaced: other assignments on the same line stay
    (review round 9 F7), and ``key =`` inside a quoted string is not an
    assignment. ``value`` is written as given (quote strings with
    :func:`fortran_string`).
    """
    hits = [m for m in re.finditer(rf"(?<![\w%]){re.escape(key)}\s*=\s*",
                                   _mask_quoted(text), re.IGNORECASE)]
    if len(hits) != 1:
        raise ValueError(f"namelist key {key}: expected once, found {len(hits)}")
    start = hits[0].end()
    m = _VALUE.match(text, start)
    end = m.end() if m else start
    return text[:start] + value + text[end:]


def relocate_case(src: Path, dst: Path, nml: str = "m2_run.nml",
                  end_date: str | None = None) -> Path:
    """Copy a staged case to ``dst`` (its output left behind), point its
    namelist at ``dst/input`` and ``dst/output``, and optionally set
    END_DATE.

    Everything is checked and the namelist rendered before anything is
    written; ``dst`` may not overlap ``src`` (after resolving links), and the
    copy is built beside ``dst`` and moved into place (review of the extend
    tools, round 7 F6).
    """
    src, dst = Path(src).resolve(), Path(dst).resolve()
    if dst == src or src in dst.parents or dst in src.parents:
        raise ValueError(f"the copy {dst} overlaps the case {src}")
    check_fvcom_dirs(dst / "input", dst / "output")
    text = (src / nml).read_text()
    text = set_value(text, "INPUT_DIR", fortran_string(f"{dst / 'input'}/"))
    text = set_value(text, "OUTPUT_DIR", fortran_string(f"{dst / 'output'}/"))
    if end_date is not None:
        text = set_value(text, "END_DATE", fortran_string(end_date))
    dst.parent.mkdir(parents=True, exist_ok=True)
    import tempfile

    tmp = Path(tempfile.mkdtemp(dir=dst.parent, prefix=f".{dst.name}."))
    # removed only after a known outcome; an exit in the middle of putting
    # the previous copy back keeps it in tmp/previous (review round 9 F3)
    state = "staging"
    try:
        work = tmp / "case"
        shutil.copytree(src, work, ignore=shutil.ignore_patterns("output"))
        (work / "output").mkdir(exist_ok=True)
        (work / nml).write_text(text)
        # the previous destination is moved aside, not deleted, until the new
        # one is in place, and moved back if that fails (review round 8 F3)
        prev = tmp / "previous"
        if dst.exists():
            dst.rename(prev)
        state = "publishing"
        try:
            work.rename(dst)
        except BaseException:
            state = "restoring"
            if prev.exists():
                try:
                    prev.rename(dst)
                except OSError:
                    raise OSError(f"could not put {dst} back; it is in {prev}") from None
            state = "restored"
            raise
        state = "done"
    finally:
        if state in ("staging", "done", "restored"):
            shutil.rmtree(tmp, ignore_errors=True)
    return dst / nml
