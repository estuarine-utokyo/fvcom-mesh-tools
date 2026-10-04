"""``fmesh-check-run``: did an FVCOM run finish, and is its output usable?

A zero exit code is not success for FVCOM: some STOP paths return 0, and a
log free of the usual error words proves nothing about the output. A run
passes only if **all** of these hold:

1. its log ends the run (``TADA``) and names no fatal condition;
2. the history output -- ``<casename>_0001.nc``, ``_0002.nc``, ... numbered
   without a gap; restart and other files are not history -- carries
   ``zeta``, ``ua``, ``va`` and ``Times``; its times increase record by
   record with no gap over 1.5 output intervals (``NC_OUT_INTERVAL``, in
   seconds/minutes/hours/days or cycles of the internal step; a declaration
   that cannot be read fails); it starts at ``NC_FIRST_OUT`` (else
   ``START_DATE``) and reaches ``END_DATE``, each within one interval (one
   hour when none is declared);
3. ``zeta``, ``ua`` and ``va`` are finite in every record.

With ``--marker PATH`` the verdict is also written as JSON to ``PATH`` --
only on success -- so that a later batch job can require it: NQSV starts a
dependent job when its predecessor ENDS, whatever the outcome, and the
marker is how the chain knows the difference (review F5, F6).

    fmesh-check-run RUN_DIR [--log fvcom.log] [--nml m2_run.nml] [--marker RUN_DIR/RUN_OK]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

_FATAL = re.compile(r"fatal|non[ -]?finite|floating.*exception|segmentation|nan detected",
                    re.IGNORECASE)


def _parse_time(text: str) -> datetime:
    text = text.strip().replace("T", " ").rstrip("Z")
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    raise ValueError(f"cannot read a time from {text!r}")


def _strip_comments(text: str) -> str:
    """A namelist without its ``!`` comments -- those outside quotes only.

    A commented-out setting read as the active one made incomplete history
    pass: ``! NC_OUT_INTERVAL = 'days = 1'`` above the real hourly interval
    widened the gap gate to a day (review 4, U1).
    """
    out = []
    for line in text.splitlines():
        quote, cut = None, len(line)
        for i, ch in enumerate(line):
            if quote:
                if ch == quote:
                    quote = None
            elif ch in "'\"":
                quote = ch
            elif ch == "!":
                cut = i
                break
        out.append(line[:cut])
    return "\n".join(out)


def _mask_quoted(text: str) -> str:
    """``text`` with the inside of every quoted value blanked, same length,
    so that a search for ``KEY =`` cannot land inside a string."""
    out, quote, i = list(text), None, 0
    while i < len(text):
        ch = text[i]
        if quote:
            if ch == quote and i + 1 < len(text) and text[i + 1] == quote:
                out[i] = out[i + 1] = " "          # a doubled delimiter
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


def _nml_value(text: str, key: str) -> str | None:
    """A namelist value, single- or double-quoted or bare; None when absent.

    Inside a quoted value a doubled delimiter is one character, as Fortran
    writes it (``'run''s'`` is ``run's``; review round 7 F3). Only real
    assignments count -- ``KEY =`` inside another quoted value does not --
    and a key assigned twice with different values is refused (round 8 F8).
    """
    value = re.compile(r"\s*=\s*(?:'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"|([^,\s/]+))")
    found = []
    for m in re.finditer(rf"(?<![\w%]){re.escape(key)}(?=\s*=)", _mask_quoted(text),
                         re.IGNORECASE):
        v = value.match(text, m.end())
        if v is None:
            continue
        single, double, bare = v.groups()
        if single is not None:
            found.append(single.replace("''", "'").strip())
        elif double is not None:
            found.append(double.replace('""', '"').strip())
        else:
            found.append(bare.strip())
    if len(set(found)) > 1:
        raise ValueError(f"{key} is set more than once, to {found}")
    return found[0] if found else None


def _grid_counts(path: Path) -> dict[str, int] | None:
    """``{"node": NP, "nele": NE}`` from an FVCOM ``_grd.dat`` header, or None
    when the file cannot be read as one."""
    try:
        with open(path) as fh:
            head = [next(fh) for _ in range(2)]
        np_, ne_ = (int(line.split("=")[1]) for line in head)
    except (OSError, StopIteration, IndexError, ValueError):
        return None
    if np_ < 3 or ne_ < 1:
        return None
    return {"node": np_, "nele": ne_}


def _staged_mesh(grid_path: Path, dep_path: Path | None):
    """(nodes, elements, depths or None) of the staged case, or a reason."""
    from fvcom_mesh_tools.io.fvcom_native import read_grd

    try:
        nodes, elements = read_grd(grid_path)
    except (OSError, ValueError) as exc:
        return f"the staged grid {grid_path} cannot be read: {exc}"
    depths = None
    if dep_path is not None:
        from fvcom_mesh_tools.io.fvcom_native import read_dep

        # the strict reader: header, count, rows and finite values, as
        # FVCOM's READ_DEPTH requires (review round 21 F4)
        try:
            dxy, depths = read_dep(dep_path)
        except (OSError, ValueError, IndexError) as exc:
            return f"the staged depths {dep_path} cannot be read: {exc}"
        if dxy.shape != nodes[:, :2].shape or not np.allclose(dxy, nodes[:, :2],
                                                              rtol=0, atol=1e-6):
            return f"the staged depths {dep_path} do not match the staged grid"
    return nodes, elements, depths


def _close(var, got, want) -> bool:
    """``got`` (a history variable) equals ``want`` (staged doubles) to the
    precision it is stored in: float32 output keeps ~7 digits (round 20 F5)."""
    if var.dtype == np.float32:
        tol = np.maximum(1e-3, 2 * np.spacing(np.abs(want).astype(np.float32)).astype(float))
    else:
        tol = 1e-3
    return got.shape == want.shape and bool(np.all(np.abs(got - want) <= tol))


def _grid_identity(ds, staged, name: str) -> list[str]:
    """Reasons the history ``ds`` is not on the staged mesh: node
    coordinates, each element's nodes (in any order) and, when the namelist
    names the depth file, the bathymetry h (review rounds 19 F7, 20 F2-F5).
    Masked values are refused, not unmasked."""
    nodes, elements, depths = staged
    out = []
    for v, want in (("x", nodes[:, 0]), ("y", nodes[:, 1]),
                    *((("h", depths),) if depths is not None else ())):
        if v not in ds.variables:
            out.append(f"{name}: no {v} to tie the history to the staged case")
            continue
        raw = ds[v][:]
        if np.ma.is_masked(raw):
            out.append(f"{name}: {v} has masked values")
            continue
        got = np.asarray(np.ma.getdata(raw), float)
        if not (np.isfinite(got).all() and _close(ds[v], got, want)):
            out.append(f"{name}: {v} is not the staged case's"
                       + (" (another bathymetry)" if v == "h" else ""))
    if "nv" not in ds.variables:
        out.append(f"{name}: no nv to tie the history to the staged mesh")
    else:
        raw = ds["nv"][:]
        if np.ma.is_masked(raw):
            out.append(f"{name}: nv has masked values")
        else:
            nv = np.asarray(np.ma.getdata(raw))
            if nv.dtype.kind not in "iu" or nv.ndim != 2 or nv.shape[0] != 3:
                out.append(f"{name}: nv is not an integer (3, nele) array")
            elif nv.T.shape != elements.shape or not np.array_equal(
                    np.sort(nv.T - 1, axis=1), np.sort(elements, axis=1)):
                out.append(f"{name}: connectivity is not the staged mesh's")
    return out


def _fortran_float(text: str) -> float:
    """A Fortran real as a namelist writes it: ``1.``, ``1.5d0``, ``10``."""
    t = text.strip().lower().replace("d", "e")
    return float(t + "0" if t.endswith(".") else t)


def _interval(value: str, text: str) -> timedelta:
    """``NC_OUT_INTERVAL``: seconds/minutes/hours/days, or cycles of the
    internal step (EXTSTEP_SECONDS x ISPLIT), as FVCOM's mod_set_time reads it."""
    m = re.fullmatch(r"\s*(seconds|minutes|hours|days|cycles)\s*=\s*([0-9.eE+-]+)\s*",
                     value, re.IGNORECASE)
    if not m:
        raise ValueError(f"cannot read NC_OUT_INTERVAL = {value!r}")
    unit, n = m.group(1).lower(), float(m.group(2))
    if unit == "cycles":
        dte, isplit = _nml_value(text, "EXTSTEP_SECONDS"), _nml_value(text, "ISPLIT")
        if dte is None or isplit is None:
            raise ValueError("NC_OUT_INTERVAL in cycles needs EXTSTEP_SECONDS and ISPLIT")
        # a finite positive step and a positive whole ISPLIT, checked before
        # any arithmetic (review round 22 F4)
        step, split = _fortran_float(dte), _fortran_float(isplit)
        if not (math.isfinite(step) and step > 0):
            raise ValueError(f"EXTSTEP_SECONDS = {dte!r} is not a finite positive step")
        if not (math.isfinite(split) and split >= 1 and split == int(split)):
            raise ValueError(f"ISPLIT = {isplit!r} is not a positive whole number")
        seconds = n * step * int(split)
    else:
        seconds = n * {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}[unit]
    # finite and positive, or the gap and coverage checks mean nothing
    # (review round 21 F9)
    if not (math.isfinite(seconds) and seconds > 0):
        raise ValueError(f"NC_OUT_INTERVAL = {value!r} is not a finite positive interval")
    try:
        out = timedelta(seconds=seconds)
    except OverflowError:
        raise ValueError(f"NC_OUT_INTERVAL = {value!r} is too long") from None
    if out <= timedelta(0):                 # below a microsecond
        raise ValueError(f"NC_OUT_INTERVAL = {value!r} rounds to zero")
    return out


def check_run(run_dir, *, log="fvcom.log", nml="m2_run.nml", casename=None) -> dict:
    """The verdict on one run directory, with the reasons for any failure.

    The history output is ``<casename>_NNNN.nc`` -- the casename defaults to
    the namelist's name without ``_run.nml`` -- numbered from 0001 without a
    gap. Restart and other files are not history, and cannot stand in for it
    (review 3, T1).
    """
    import netCDF4

    run = Path(run_dir).resolve()
    reasons: list[str] = []
    info: dict = {"run_dir": str(run)}
    case = casename or re.sub(r"_run\.nml$", "", nml)

    log_path = run / log
    if not log_path.exists():
        reasons.append(f"no log {log_path.name}")
    else:
        text = log_path.read_text(errors="replace")
        if "TADA" not in text:
            reasons.append("the log does not end the run (no TADA)")
        bad = sorted({m.group(0).lower() for m in _FATAL.finditer(text)})
        if bad:
            reasons.append(f"the log reports: {', '.join(bad)}")

    end = first = interval = None
    nml_path = run / nml
    nml_text = ""
    start_parsed = None
    vals: dict[str, str | None] = {}
    if not nml_path.exists():
        reasons.append(f"no namelist {nml}")
    else:
        nml_text = _strip_comments(nml_path.read_text())
        # a key set twice, differently, is a failure (review round 8 F8)
        for key in ("END_DATE", "NC_FIRST_OUT", "START_DATE", "NC_OUT_INTERVAL",
                    "OUTPUT_DIR", "GRID_FILE", "INPUT_DIR", "DEPTH_FILE"):
            try:
                vals[key] = _nml_value(nml_text, key)
            except ValueError as exc:
                reasons.append(str(exc))
                vals[key] = None
        # an unreadable date is a failure reason, not a crash (round 21 F9)
        def _date(key):
            try:
                return _parse_time(vals[key]) if vals.get(key) else None
            except (ValueError, OverflowError) as exc:
                reasons.append(f"{key} = {vals[key]!r} cannot be read ({exc})")
                return None

        end = _date("END_DATE")
        # START_DATE on its own, always: a valid NC_FIRST_OUT used to hide a
        # bad one (review round 22 F2)
        start_parsed = _date("START_DATE")
        first = _date("NC_FIRST_OUT") or start_parsed
        v = vals["NC_OUT_INTERVAL"]
        if v is not None:
            # present but unreadable fails closed: it used to switch the gap
            # check off (review 3, T3)
            try:
                interval = _interval(v, nml_text)
                info["output_interval_s"] = interval.total_seconds()
            except ValueError as exc:
                reasons.append(str(exc))
    # the history is where the namelist sends it (review round 6 F3); a
    # namelist without OUTPUT_DIR falls back to <run>/output
    outdir = run / "output"
    od = vals.get("OUTPUT_DIR")
    if od:
        outdir = Path(od) if Path(od).is_absolute() else run / od
    info["output_dir"] = str(outdir)
    # the staged mesh's size, when the namelist names a readable grid file:
    # the history must be on it (review of the extend tools, round 4 F14)
    staged = None
    staged_grid = None
    grid, indir = vals.get("GRID_FILE"), vals.get("INPUT_DIR")
    if grid and indir:
        gpath = Path(indir) / grid if Path(indir).is_absolute() else run / indir / grid
        staged = _grid_counts(gpath)
        if staged is not None:
            dfile = vals.get("DEPTH_FILE")
            dpath = None
            if dfile:
                dpath = Path(indir) / dfile if Path(indir).is_absolute() else run / indir / dfile
            staged_grid = _staged_mesh(gpath, dpath)
            if isinstance(staged_grid, str):
                reasons.append(staged_grid)
                staged_grid = None
        # a grid named but unreadable is a failure, not a skipped check
        # (review round 5 F6)
        if staged is None:
            reasons.append(f"the namelist's grid {gpath} cannot be read for its counts")
    # a run must integrate over a positive interval: END_DATE after
    # START_DATE (output may start as late as END_DATE; review round 12 F5)
    start = start_parsed if nml_path.exists() else None   # parsed once, above
    if nml_path.exists() and start is None and not vals.get("START_DATE"):
        reasons.append(f"no START_DATE found in {nml}")
    if end is not None and start is not None and end <= start:
        reasons.append(f"END_DATE {end} is not after the start {start}")
    # FVCOM refuses a first output outside the run (round 23 F5)
    nfo = vals.get("NC_FIRST_OUT")
    if nfo and start is not None and end is not None and first is not None and not (
            start <= first <= end):
        reasons.append(f"NC_FIRST_OUT {first} is outside [{start}, {end}]")
    if end is None:
        reasons.append(f"no END_DATE found in {nml}")
    else:
        info["end_date"] = end.isoformat(sep=" ")
    tol = interval if interval is not None else timedelta(hours=1)

    # the history stacks, in number order, from 0001 without a gap (T2)
    pat = re.compile(rf"^{re.escape(case)}_(\d{{4}})\.nc$")
    history = sorted((int(m.group(1)), f) for f in outdir.glob("*.nc")
                     if (m := pat.match(f.name)))
    if not history:
        reasons.append(f"no history output ({case}_0001.nc ...) in {outdir}")
    elif [k for k, _ in history] != list(range(1, len(history) + 1)):
        reasons.append("the history stacks are not numbered 0001.. without a gap: "
                       + ", ".join(f.name for _, f in history))
    stamps: list = []
    sizes: dict[str, int] = {}
    for _, f in history:
        try:
            with netCDF4.Dataset(f) as ds:
                missing = [v for v in ("zeta", "ua", "va", "Times") if v not in ds.variables]
                if missing:
                    reasons.append(f"{f.name}: the history output lacks {', '.join(missing)}")
                    continue
                # Times indexed by the fields' record dimension (round 24 F6)
                tv = ds["Times"]
                if tv.ndim != 2 or tv.dimensions[0] != "time" or tv.dtype.kind not in "SU":
                    reasons.append(f"{f.name}: Times is {tv.dimensions} {tv.dtype}, not "
                                   f"(time, DateStrLen) characters")
                    continue
                raw_times = ds["Times"][:]
                # a masked character is unknown, not its fill (round 23 F4)
                if np.ma.is_masked(raw_times):
                    reasons.append(f"{f.name}: Times has masked characters")
                    continue
                times = [_parse_time(str(t)) for t in
                         np.atleast_1d(netCDF4.chartostring(raw_times))]
                if not times:
                    reasons.append(f"{f.name}: no records")
                    continue
                stamps += times
                # the history must be on the staged mesh itself, not one of
                # the same size (review round 19 F7): node coordinates, and
                # each element's nodes in any order (FVCOM may reverse them)
                if staged_grid is not None:              # every stack (round 20 F2)
                    reasons += _grid_identity(ds, staged_grid, f.name)
                for var, dim in (("zeta", "node"), ("ua", "nele"), ("va", "nele")):
                    # the layout FVCOM writes; a square transposed array
                    # has the right shape and the wrong meaning (round 13 F4)
                    if ds[var].dimensions != ("time", dim):
                        reasons.append(f"{f.name}: {var} has dimensions "
                                       f"{ds[var].dimensions}, not ('time', '{dim}')")
                        continue
                    a = np.ma.filled(ds[var][:], np.nan)
                    # time by space, with space not empty and the same in
                    # every stack and on the staged mesh (round 4 F14)
                    if a.ndim != 2 or a.shape[1] == 0:
                        reasons.append(f"{f.name}: {var} has shape {a.shape}, not "
                                       f"(time, {dim}) with {dim} > 0")
                        continue
                    if sizes.setdefault(dim, a.shape[1]) != a.shape[1]:
                        reasons.append(f"{f.name}: {var} has {a.shape[1]} {dim}, other "
                                       f"output {sizes[dim]}")
                    if staged is not None and a.shape[1] != staged[dim]:
                        reasons.append(f"{f.name}: {var} has {a.shape[1]} {dim}, the staged "
                                       f"mesh {staged[dim]}")
                    if a.shape[0] != len(times):
                        reasons.append(f"{f.name}: {var} has {a.shape[0]} records for "
                                       f"{len(times)} times")
                    elif not np.isfinite(a).all():
                        reasons.append(f"{f.name}: {var} is not finite everywhere")
                # every other floating field written per record -- 3-D
                # velocities, scalars, turbulence -- must be finite too; a
                # failed 3-D solution can sit under healthy barotropic fields
                # (review round 10 F4). One record at a time, to bound memory.
                for name, v in ds.variables.items():
                    if (name in ("zeta", "ua", "va") or "time" not in v.dimensions
                            or v.dtype.kind != "f"):
                        continue
                    # FVCOM writes time first; a field with time elsewhere is
                    # refused, not skipped (review round 12 F4)
                    if v.dimensions[0] != "time":
                        reasons.append(f"{f.name}: {name} has time at position "
                                       f"{v.dimensions.index('time')} of {v.dimensions}")
                        continue
                    if v.shape[0] != len(times):
                        reasons.append(f"{f.name}: {name} has {v.shape[0]} records for "
                                       f"{len(times)} times")
                        continue
                    # an empty record is finite only vacuously (round 11 F4)
                    if any(n == 0 for n in v.shape[1:]):
                        reasons.append(f"{f.name}: {name} has shape {v.shape}, an empty "
                                       f"dimension")
                        continue
                    for dim, n in zip(v.dimensions[1:], v.shape[1:]):
                        if dim in ("node", "nele") and n != sizes.get(dim, n):
                            reasons.append(f"{f.name}: {name} has {n} {dim}, the history "
                                           f"{sizes[dim]}")
                    for k in range(v.shape[0]):
                        if not np.isfinite(np.ma.filled(v[k], np.nan)).all():
                            reasons.append(f"{f.name}: {name} is not finite at record {k}")
                            break
        except Exception as exc:        # unreadable, truncated, not NetCDF
            reasons.append(f"{f.name}: cannot be read ({exc.__class__.__name__}: {exc})")
    info["n_records"] = len(stamps)
    info["first_output"] = stamps[0].isoformat(sep=" ") if stamps else None
    info["last_output"] = stamps[-1].isoformat(sep=" ") if stamps else None
    if stamps:
        # in file order, never re-sorted: a reversed or overlapping stack is
        # a defect to report, not to repair (T2)
        # microseconds, not seconds: sub-second records are distinct
        # (review round 14 F7)
        steps = np.diff(np.array(stamps, dtype="datetime64[us]")).astype(float) / 1e6
        # one record has no cadence to judge; its coverage still is (U2)
        if steps.size and (steps <= 0).any():
            reasons.append("the history times do not increase record by record")
        elif steps.size and interval is not None and steps.max() > 1.5 * tol.total_seconds():
            reasons.append(f"the output has a gap of {steps.max() / 3600:.2f} h against "
                           f"an interval of {tol.total_seconds() / 3600:.2f} h")
        if first is not None and abs(stamps[0] - first) > tol:
            reasons.append(f"the output starts at {stamps[0]}, not at {first}")
        if end is not None and stamps[-1] < end - tol:
            reasons.append(f"the output stops at {stamps[-1]} before END_DATE {end}")
        # nor run past it: records after END_DATE are another run's (round 21 F1)
        if end is not None and stamps[-1] > end + tol:
            reasons.append(f"the output runs to {stamps[-1]}, past END_DATE {end}")
        # some output after the start: an initial record alone is no
        # integration, however short the run (review round 13 F5)
        if start is not None and stamps[-1] <= start:
            reasons.append(f"no output after the start {start}")
    info["ok"] = not reasons
    info["reasons"] = reasons
    return info


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="fmesh-check-run",
                                description="Decide whether an FVCOM run finished "
                                            "with usable output.")
    p.add_argument("run_dir", type=Path)
    p.add_argument("--log", default="fvcom.log", help="log file in RUN_DIR")
    p.add_argument("--nml", default="m2_run.nml", help="namelist in RUN_DIR")
    p.add_argument("--casename", default=None,
                   help="history files are <casename>_NNNN.nc (default: the namelist's "
                        "name without _run.nml)")
    p.add_argument("--marker", type=Path, default=None,
                   help="write the verdict here, on success only")
    return p


def _is_verdict(path: Path) -> bool:
    """True for a small JSON file of the shape ``check_run`` returns."""
    import os
    import stat

    try:
        # a regular file only, opened without following a link or blocking on a
        # FIFO, and read to at most the size limit (review round 54 F1)
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        try:
            st = os.fstat(fd)
            if not stat.S_ISREG(st.st_mode) or st.st_size > 1_000_000:
                return False
            with os.fdopen(fd, "rb", closefd=False) as fh:
                data = fh.read(1_000_001)
        finally:
            os.close(fd)
        if len(data) > 1_000_000:
            return False
        doc = json.loads(data.decode())
    except (OSError, ValueError):
        return False
    return isinstance(doc, dict) and isinstance(doc.get("ok"), bool) and "reasons" in doc


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    # A marker from an earlier attempt must not survive this one's failure
    # (review 2, R1): it goes first, and comes back only on success.
    if args.marker is not None:
        # only a verdict this tool wrote may be removed: whatever else the marker
        # path names (an input, a history file, a link to one) is left alone, so the
        # guard does not depend on knowing every input (review rounds 51 F3, 52 F3, 53 F1)
        mk = args.marker
        if mk.is_symlink() or (os.path.lexists(mk) and not _is_verdict(mk)):
            print(f"[check-run] --marker {mk} exists and is not a verdict written by this "
                  "tool; refusing to remove it", file=sys.stderr)
            return 2
        args.marker.unlink(missing_ok=True)
    info = check_run(args.run_dir, log=args.log, nml=args.nml, casename=args.casename)
    verdict = "OK" if info["ok"] else "FAILED"
    print(f"[check-run] {args.run_dir}: {verdict}; {info['n_records']} record(s), "
          f"last {info['last_output']}, END_DATE {info.get('end_date')}")
    for r in info["reasons"]:
        print(f"[check-run]   {r}")
    if info["ok"] and args.marker is not None:
        args.marker.write_text(json.dumps(info, indent=1) + "\n")
    return 0 if info["ok"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
