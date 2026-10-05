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
    mkdir -p "$d"; cp -r "$SRC/input" "$SRC/output" "$SRC"/*_run.nml "$d/"
    sed -i "s#^ *INPUT_DIR .*#INPUT_DIR = '$d/input/',#; s#^ *OUTPUT_DIR .*#OUTPUT_DIR = '$d/output/',#" "$d"/*_run.nml
    # history and restart output, so that the long OUTPUT_DIR is written to as well as read
    sed -i "s#^ *NC_ON .*#NC_ON = T,#; s#^ *NC_OUT_INTERVAL .*#NC_OUT_INTERVAL = 'seconds=60.0',#; s#^ *RST_ON .*#RST_ON = T,#; s#^ *RST_FIRST_OUT .*#RST_FIRST_OUT = '2020-01-01 00:01:00',#; s#^ *RST_OUT_INTERVAL .*#RST_OUT_INTERVAL = 'seconds=60.0',#" "$d"/*_run.nml
    echo "[$1] dir length $(( ${#d} + 1 )) bytes"
    ( cd "$d" && mpiexec -np 4 "$2" --casename=tokyo_bay_v1_smoke > fvcom.log 2>&1 ) && rc=0 || rc=$?
    echo "[$1] exit=$rc; netcdf files: $(ls "$d/output" | grep -c '\.nc$')"
    echo "[$1] $(grep -c -i 'fatal' "$d/fvcom.log") fatal lines; $(grep -c TADA "$d/fvcom.log") TADA"
}
run_case old_short "$OLD" "$ROOT/o"
run_case new_short "$NEW" "$ROOT/n"
run_case new_long  "$NEW" "$ROOT/$LONGSEG/$LONGSEG/case"
run_case old_long  "$OLD" "$ROOT/$LONGSEG/$LONGSEG/oldcase" || true
module purge; unset LD_LIBRARY_PATH
. "${FMESH_CONDA_ROOT:?}/etc/profile.d/conda.sh"; set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
python - "$ROOT" "$LONGSEG" <<'PY'
import sys, glob, numpy as np, netCDF4 as nc
root, seg = sys.argv[1:]
def first(d):
    f = sorted(glob.glob(f"{d}/output/*.nc")); return f[0] if f else None
ref = first(f"{root}/o")
for lab, d in (("new_short", f"{root}/n"), ("new_long", f"{root}/{seg}/{seg}/case")):
    f = first(d)
    if not (ref and f): print(lab, "NO OUTPUT"); continue
    a, b = nc.Dataset(ref), nc.Dataset(f)
    diffs = {v: float(np.max(np.abs(a[v][:] - b[v][:]))) for v in ("zeta", "u", "v", "temp", "salinity") if v in a.variables}
    print(lab, "max |diff| vs old_short:", diffs)
PY
