#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=4
#PJM -L elapse=00:15:00
#PJM -j
#PJM -X
#PJM -N fmesh_longpath
#============================================================================
# Does FVCOM keep a long run-directory path? Runs the same tiny case (the
# 2-minute zero-forcing smoke in outputs/fvcom_smoke) with the OLD and the NEW
# FVCOM binary, in a short and in a >200-byte directory, with absolute
# INPUT_DIR / OUTPUT_DIR, and compares the results.
#
#   pjsub -x NEW=<new fvcom>,OLD=<old fvcom>,FVCOM_LIBS=<lib prefix> \
#       jobs/genkai/fvcom_longpath_test.sh
# Work dir: scratch/longpath/ (git-ignored). Log: logs/fvcom_longpath.<jobid>.log
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
: "${NEW:?set NEW}"; : "${OLD:?set OLD}"; : "${FVCOM_LIBS:?set FVCOM_LIBS}"
. jobs/genkai/common.sh fvcom_longpath 1
set +u; conda deactivate; set -u
module load intel/2025.1.3 impi/2021.15 netcdf/4.9.2 netcdf-fortran/4.6.1 hdf5/1.14.4
export OMP_NUM_THREADS=1 INSTALLDIR=$FVCOM_LIBS
export LD_LIBRARY_PATH="$FVCOM_LIBS/lib:$FVCOM_LIBS/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
SRC=$PWD/outputs/fvcom_smoke
ROOT=$PWD/scratch/longpath
rm -rf "$ROOT"; mkdir -p "$ROOT"
LONGSEG=$(printf 'a_deliberately_long_directory_name_%.0s' 1 2 3 4 5)
run_case () {   # $1 label  $2 binary  $3 directory
    local d=$3
    mkdir -p "$d/output"; cp -r "$SRC/input" "$SRC"/*_run.nml "$d/"
    sed -i "s#^ *INPUT_DIR .*#INPUT_DIR = '$d/input/',#; s#^ *OUTPUT_DIR .*#OUTPUT_DIR = '$d/output/',#" "$d"/*_run.nml
    # history, average, surface and restart output (their identifiers share the directory prefix), so that the long OUTPUT_DIR is written to as well as read
    sed -i "s#^ *NC_ON .*#NC_ON = T,#; s#^ *NC_OUT_INTERVAL .*#NC_OUT_INTERVAL = 'seconds=60.0',#; s#^ *NCAV_ON .*#NCAV_ON = T,#; s#^ *NCAV_OUT_INTERVAL .*#NCAV_OUT_INTERVAL = 'seconds=60.0',#; s#^ *NCSF_ON .*#NCSF_ON = T,#; s#^ *NCSF_OUT_INTERVAL .*#NCSF_OUT_INTERVAL = 'seconds=60.0',#; s#^ *RST_ON .*#RST_ON = T,#; s#^ *RST_FIRST_OUT .*#RST_FIRST_OUT = '2020-01-01 00:01:00',#; s#^ *RST_OUT_INTERVAL .*#RST_OUT_INTERVAL = 'seconds=60.0',#" "$d"/*_run.nml
    echo "[$1] dir length $(( ${#d} + 1 )) bytes"
    local extra=()
    case "$1" in new_*) extra=("--logfile=$d/fvcom_own.log");; esac   # an absolute log file name is read from the command line
    ( cd "$d" && mpiexec -np 4 "$2" --casename=tokyo_bay_v1_smoke "${extra[@]}" > fvcom.log 2>&1 ) && rc=0 || rc=$?
    case "$1" in new_*) [ -s "$d/fvcom_own.log" ] || { echo "[$1] --logfile was not written to $d"; return 1; };; esac
    echo "[$1] exit=$rc; netcdf files: $(ls "$d/output" | grep -c '\.nc$')"
    echo "[$1] $(grep -c -i 'fatal' "$d/fvcom.log") fatal lines; $(grep -c TADA "$d/fvcom.log") TADA"
    return $rc
}
run_case old_short "$OLD" "$ROOT/o"
run_case new_short "$NEW" "$ROOT/n"
run_case new_long  "$NEW" "$ROOT/$LONGSEG/$LONGSEG/case"
# the old binary is expected to FAIL in the long directory (it cuts the path at 80 bytes)
if run_case old_long "$OLD" "$ROOT/$LONGSEG/$LONGSEG/oldcase"; then echo "[old_long] UNEXPECTEDLY ran"; exit 1; fi
module purge; unset LD_LIBRARY_PATH
. "${FMESH_CONDA_ROOT:?}/etc/profile.d/conda.sh"; set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
python - "$ROOT" "$LONGSEG" <<'PY'
import sys, glob, os, numpy as np, netCDF4 as nc
root, seg = sys.argv[1:]
cases = {"old_short": f"{root}/o", "new_short": f"{root}/n", "new_long": f"{root}/{seg}/{seg}/case"}
def files(d):
    return {os.path.basename(f): f for f in glob.glob(f"{d}/output/*.nc")}
ref = files(cases["old_short"])
bad = []
expected = {"tokyo_bay_v1_smoke_0001.nc", "tokyo_bay_v1_smoke_avg_0001.nc", "tokyo_bay_v1_smoke_surface_0001.nc", "tokyo_bay_v1_smoke_restart_0001.nc"}
if not expected <= set(ref):
    bad.append(f"reference run lacks {sorted(expected - set(ref))}")
for lab in ("new_short", "new_long"):
    got = files(cases[lab])
    if set(got) != set(ref):
        bad.append(f"{lab}: files {sorted(got)} != {sorted(ref)}"); continue
    worst = 0.0
    for name, f in ref.items():
        a, b = nc.Dataset(f), nc.Dataset(got[name])
        if set(a.variables) != set(b.variables) or set(a.dimensions) != set(b.dimensions):
            bad.append(f"{lab}/{name}: variables or dimensions differ"); continue
        for d in a.dimensions:
            if a.dimensions[d].isunlimited() != b.dimensions[d].isunlimited():
                bad.append(f"{lab}/{name}: dimension {d} unlimited status differs")
            if len(a.dimensions[d]) != len(b.dimensions[d]):
                bad.append(f"{lab}/{name}: dimension {d} has {len(a.dimensions[d])} vs {len(b.dimensions[d])}")
        for v in a.variables:
            if v == "file_date":     # the creation time of the file: differs by design
                continue
            x, y = a[v][:], b[v][:]
            if a[v].dimensions != b[v].dimensions or np.shape(x) != np.shape(y) or a[v].dtype != b[v].dtype:
                bad.append(f"{lab}/{name}/{v}: shape or type differs"); continue
            if np.ma.getmaskarray(x).tolist() != np.ma.getmaskarray(y).tolist():
                bad.append(f"{lab}/{name}/{v}: masks differ"); continue
            x, y = np.ma.getdata(x), np.ma.getdata(y)
            same = np.array_equal(x, y, equal_nan=x.dtype.kind in "fc")   # a NaN on one side only differs
            if same and x.dtype.kind in "fc" and not np.all(np.isfinite(x)) and name.endswith("_0001.nc"):
                bad.append(f"{lab}/{name}/{v}: non-finite values in both files")
            fin = np.isfinite(x - y) if x.dtype.kind in "fc" else np.ones(x.shape, bool)
            diff = 0.0 if same else (float(np.max(np.abs((x - y)[fin]))) if fin.any() else float("inf"))
            if not same:
                bad.append(f"{lab}/{name}/{v}: differs (max {diff}, {int(np.sum(x != y))} points)")
            worst = max(worst, diff)
    print(lab, "files", sorted(got), "max |diff| over all variables of all files:", worst)
    if worst != 0.0:
        bad.append(f"{lab}: differs from old_short by {worst}")
print("RESULT:", "FAIL " + "; ".join(bad) if bad else "OK")
sys.exit(1 if bad else 0)
PY
