#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L node=1
#PJM -L elapse=02:00:00
#PJM -j
#PJM -X
#PJM -N fmesh_448smoke
#============================================================================
# GENKAI twin of jobs/octopus/448_extend_smoke.sh: does a wide mesh RUN?
# Stage one case (notebooks/448_extend_smoke.py: uniform M2 on the open
# boundary, stability only), integrate it for FMESH_DAYS, and judge it as 423
# does: exit code, the log grepped for fatal signatures, and fmesh-check-run
# (TADA, output reaching END_DATE, finite fields).
#
#   pjsub -x FMESH_CASE=<case prefix>,FMESH_RUN_ROOT=<new dir>,FVCOM_LIBS=<lib prefix> \
#       jobs/genkai/448_extend_smoke.sh
# Optional: FMESH_DAYS (2), FMESH_RANKS (120 = one whole node: ~225 elements per rank on
#   27,000 elements, the 'scale up' regime of the rank rule; the shared pool had no free
#   node when this was written, so the job asks for a whole exclusive node),
#   FMESH_GAUGE (MERA), FMESH_FVCOM.
# FMESH_RUN_ROOT: 448 refuses a used root; SMOKE_OK is never removed (review r2 F21).
# Keep it inside this repository (scratch/, git-ignored).
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
: "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}"
: "${FMESH_CASE:?set FMESH_CASE}"
: "${FVCOM_LIBS:?set FVCOM_LIBS (the FVCOM library install prefix)}"
. jobs/genkai/common.sh 448_extend_smoke 1
RUN_ROOT=$FMESH_RUN_ROOT
RANKS=${FMESH_RANKS:-120}
FVCOM=$(fmesh_fvcom)
echo "[448] fvcom $FVCOM $(sha256sum "$FVCOM" | cut -c1-16) ranks=$RANKS"
python notebooks/448_extend_smoke.py --case "$FMESH_CASE" --root "$RUN_ROOT" \
    --days "${FMESH_DAYS:-2}" --gauge "${FMESH_GAUGE:-MERA}"
CASE_DIR=$RUN_ROOT/extended
set +u; conda deactivate; set -u
module purge
module load intel/2025.1.3 impi/2021.15 netcdf/4.9.2 netcdf-fortran/4.6.1 hdf5/1.14.4
export INSTALLDIR=$FVCOM_LIBS OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
fail=0
t0=$(date +%s)
if ( cd "$CASE_DIR" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 > fvcom.log 2>&1 ); then
    rc=0; else rc=$?; fi
echo "[448] exit=$rc seconds=$(( $(date +%s) - t0 ))"
if grep -Ei 'fatal|non[ -]?finite|floating exception|segmentation|nan detected' \
        "$CASE_DIR/fvcom.log" >/dev/null; then
    echo "[448] UNHEALTHY log"; tail -25 "$CASE_DIR/fvcom.log"; fail=1
fi
[ "$rc" -eq 0 ] || { tail -25 "$CASE_DIR/fvcom.log"; fail=1; }
module purge
unset LD_LIBRARY_PATH
. "${FMESH_CONDA_ROOT:-$WORK_DIR/miniforge3}/etc/profile.d/conda.sh"
set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
python -m fvcom_mesh_tools.cli.check_run "$CASE_DIR" || fail=1
[ "$fail" -eq 0 ] && date -Is > "$RUN_ROOT/SMOKE_OK"
echo "end=$(date -Is) fail=$fail"
exit $fail
