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
import re
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


def _nml_value(text: str, key: str) -> str | None:
    """A namelist value, single- or double-quoted or bare; None when absent."""
    m = re.search(rf"\b{key}\s*=\s*(?:'([^']*)'|\"([^\"]*)\"|([^,\s/]+))", text,
                  re.IGNORECASE)
    if not m:
        return None
    return next(g for g in m.groups() if g is not None).strip()


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
        return timedelta(seconds=n * _fortran_float(dte) * int(_fortran_float(isplit)))
    return timedelta(seconds=n * {"seconds": 1, "minutes": 60, "hours": 3600,
                                  "days": 86400}[unit])


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
    if not nml_path.exists():
        reasons.append(f"no namelist {nml}")
    else:
        nml_text = nml_path.read_text()
        v = _nml_value(nml_text, "END_DATE")
        end = _parse_time(v) if v else None
        v = _nml_value(nml_text, "NC_FIRST_OUT") or _nml_value(nml_text, "START_DATE")
        first = _parse_time(v) if v else None
        v = _nml_value(nml_text, "NC_OUT_INTERVAL")
        if v is not None:
            # present but unreadable fails closed: it used to switch the gap
            # check off (review 3, T3)
            try:
                interval = _interval(v, nml_text)
                info["output_interval_s"] = interval.total_seconds()
            except ValueError as exc:
                reasons.append(str(exc))
    if end is None:
        reasons.append(f"no END_DATE found in {nml}")
    else:
        info["end_date"] = end.isoformat(sep=" ")
    tol = interval if interval is not None else timedelta(hours=1)

    # the history stacks, in number order, from 0001 without a gap (T2)
    pat = re.compile(rf"^{re.escape(case)}_(\d{{4}})\.nc$")
    history = sorted((int(m.group(1)), f) for f in (run / "output").glob("*.nc")
                     if (m := pat.match(f.name)))
    if not history:
        reasons.append(f"no history output ({case}_0001.nc ...)")
    elif [k for k, _ in history] != list(range(1, len(history) + 1)):
        reasons.append("the history stacks are not numbered 0001.. without a gap: "
                       + ", ".join(f.name for _, f in history))
    stamps: list = []
    for _, f in history:
        try:
            with netCDF4.Dataset(f) as ds:
                missing = [v for v in ("zeta", "ua", "va", "Times") if v not in ds.variables]
                if missing:
                    reasons.append(f"{f.name}: the history output lacks {', '.join(missing)}")
                    continue
                times = [_parse_time(str(t)) for t in
                         np.atleast_1d(netCDF4.chartostring(ds["Times"][:]))]
                if not times:
                    reasons.append(f"{f.name}: no records")
                    continue
                stamps += times
                for var in ("zeta", "ua", "va"):
                    a = np.ma.filled(ds[var][:], np.nan)
                    if a.shape[0] != len(times):
                        reasons.append(f"{f.name}: {var} has {a.shape[0]} records for "
                                       f"{len(times)} times")
                    elif not np.isfinite(a).all():
                        reasons.append(f"{f.name}: {var} is not finite everywhere")
        except Exception as exc:        # unreadable, truncated, not NetCDF
            reasons.append(f"{f.name}: cannot be read ({exc.__class__.__name__}: {exc})")
    info["n_records"] = len(stamps)
    info["first_output"] = stamps[0].isoformat(sep=" ") if stamps else None
    info["last_output"] = stamps[-1].isoformat(sep=" ") if stamps else None
    if stamps:
        # in file order, never re-sorted: a reversed or overlapping stack is
        # a defect to report, not to repair (T2)
        steps = np.diff(np.array(stamps, dtype="datetime64[s]")).astype(float)
        if (steps <= 0).any():
            reasons.append("the history times do not increase record by record")
        elif interval is not None and steps.max() > 1.5 * tol.total_seconds():
            reasons.append(f"the output has a gap of {steps.max() / 3600:.2f} h against "
                           f"an interval of {tol.total_seconds() / 3600:.2f} h")
        if first is not None and abs(stamps[0] - first) > tol:
            reasons.append(f"the output starts at {stamps[0]}, not at {first}")
        if end is not None and stamps[-1] < end - tol:
            reasons.append(f"the output stops at {stamps[-1]} before END_DATE {end}")
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


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    # A marker from an earlier attempt must not survive this one's failure
    # (review 2, R1): it goes first, and comes back only on success.
    if args.marker is not None:
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
