"""``fmesh-check-run``: did an FVCOM run finish, and is its output usable?

A zero exit code is not success for FVCOM: some STOP paths return 0, and a
log free of the usual error words proves nothing about the output. A run
passes only if **all** of these hold:

1. its log ends the run (``TADA``) and names no fatal condition;
2. every output NetCDF opens, and the last output time reaches the
   namelist's ``END_DATE`` (within one output interval);
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


def check_run(run_dir, *, log="fvcom.log", nml="m2_run.nml") -> dict:
    """The verdict on one run directory, with the reasons for any failure."""
    import netCDF4

    run = Path(run_dir).resolve()
    reasons: list[str] = []
    info: dict = {"run_dir": str(run)}

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

    end = None
    nml_path = run / nml
    if nml_path.exists():
        m = re.search(r"END_DATE\s*=\s*'([^']+)'", nml_path.read_text())
        if m:
            end = _parse_time(m.group(1))
            info["end_date"] = end.isoformat(sep=" ")
    if end is None:
        reasons.append(f"no END_DATE found in {nml}")

    files = sorted((run / "output").glob("*.nc"))
    if not files:
        reasons.append("no output NetCDF")
    last, step = None, None
    n_records = 0
    for f in files:
        try:
            with netCDF4.Dataset(f) as ds:
                times = [str(t) for t in np.atleast_1d(
                    netCDF4.chartostring(ds["Times"][:]))] \
                    if "Times" in ds.variables else []
                if times:
                    stamps = [_parse_time(t) for t in times]
                    last = max(stamps) if last is None else max(last, max(stamps))
                    if len(stamps) > 1:
                        step = stamps[-1] - stamps[-2]
                n_records += len(times)
                for var in ("zeta", "ua", "va"):
                    if var in ds.variables:
                        if not np.isfinite(np.ma.filled(ds[var][:], np.nan)).all():
                            reasons.append(f"{f.name}: {var} is not finite everywhere")
        except Exception as exc:        # unreadable, truncated, not NetCDF
            reasons.append(f"{f.name}: cannot be read ({exc.__class__.__name__}: {exc})")
    info["n_records"] = n_records
    info["last_output"] = last.isoformat(sep=" ") if last else None
    if files and last is None:
        reasons.append("the output has no Times records")
    if end is not None and last is not None:
        slack = step if step is not None else timedelta(hours=1)
        if last < end - slack:
            reasons.append(f"the output stops at {last} before END_DATE {end}")
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
    p.add_argument("--marker", type=Path, default=None,
                   help="write the verdict here, on success only")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    info = check_run(args.run_dir, log=args.log, nml=args.nml)
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
