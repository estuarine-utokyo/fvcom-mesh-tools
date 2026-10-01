#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=120GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_423smoke
#PBS -j o
#PBS -o logs/423_m2_smoke.pbs.log
#PBS -r n
# Does the mesh RUN?  Both staged cases, shortened, one after another.
#
# A zero exit code is not success for a Fortran solver, so the log is grepped
# and the output records are counted as well.  The question this answers is
# whether the case integrates at all -- not whether the answer is right,
# which is 413's job on the full run.
#
# Required: FMESH_RUN_ROOT.  Optional: FMESH_DAYS (2), FMESH_RANKS (64).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
# One mutating stage per run root (review round 10 F3): the root lock is
# taken atomically before any marker is invalidated or case staged, and no
# 412 run may hold a case lock meanwhile (412 checks for this lock in turn).
STAGE_LOCK=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/.staging
mkdir -p -- "$FMESH_RUN_ROOT"
mkdir -- "$STAGE_LOCK" 2>/dev/null \
    || { echo "another stage holds $STAGE_LOCK (remove it if that job is gone)"; exit 2; }
trap 'rmdir -- "$STAGE_LOCK"' EXIT
for c in "$FMESH_RUN_ROOT"/*/.running; do
    [ -e "$c" ] && { echo "a run is active: $c"; exit 2; }
done
# INVALIDATE first, before anything that can fail (review 2, R1)
rm -f "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/SMOKE_OK" \
      "$FMESH_RUN_ROOT/base/RUN_OK" "$FMESH_RUN_ROOT/refined/RUN_OK"
. jobs/octopus/common.sh 423_m2_smoke 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
DAYS=${FMESH_DAYS:-2}
RANKS=${FMESH_RANKS:-64}
FVCOM=$(fmesh_fvcom)
SMOKE=$RUN_ROOT/smoke
[ -f "$RUN_ROOT/STAGED" ] || { echo "not staged: $RUN_ROOT (no STAGED marker)"; exit 2; }
python - "$RUN_ROOT" "$SMOKE" "$DAYS" <<'PY'
import re, sys
from pathlib import Path
from fvcom_mesh_tools.io.fvcom_namelist import end_after, relocate_case
root, smoke, days = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
for case in ("base", "refined"):
    nml = (root / case / "m2_run.nml").read_text()
    start = re.search(r"START_DATE\s*=\s*'([^']+)'", nml).group(1)
    end = end_after(start, days)           # finite, positive (review round 12 F5)
    # the moved directories are checked against FVCOM's 80 bytes before
    # anything is written (review round 6 F11)
    relocate_case(root / case, smoke / case, end_date=end)
    print(f"[423] {case}: end -> {end}")
PY
set +u; conda deactivate; set -u
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r $init ]]; then source "$init"; break; fi
    done
fi
module purge
module load BaseCPU/2026
module load hdf5/1.14.6 netcdf-c/4.9.3 netcdf-fortran/4.6.2
INSTALLDIR=$FVCOM_LIBS
export INSTALLDIR OMP_NUM_THREADS=1 PROJ_DATA=/usr/share/proj
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
fail=0
for case in base refined; do
    t0=$(date +%s.%N)
    if ( cd "$SMOKE/$case" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 \
            > fvcom.log 2>&1 ); then rc=0; else rc=$?; fi
    t1=$(date +%s.%N)
    printf '[423] %-8s exit=%s seconds=%.1f\n' "$case" "$rc" "$(echo "$t1 - $t0"|bc)"
    if grep -Ei 'fatal|non[ -]?finite|floating exception|segmentation|nan detected' \
            "$SMOKE/$case/fvcom.log" >/dev/null; then
        echo "[423] $case: UNHEALTHY log"; tail -25 "$SMOKE/$case/fvcom.log"; fail=1
    fi
    [ "$rc" -eq 0 ] || { tail -25 "$SMOKE/$case/fvcom.log"; fail=1; }
done
# A zero exit code and a quiet log are not a finished run: some STOP paths
# return 0.  fmesh-check-run requires TADA, output reaching END_DATE and
# finite fields (review F6); it needs the conda Python back.
module purge
unset LD_LIBRARY_PATH
set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
for case in base refined; do
    python -m fvcom_mesh_tools.cli.check_run "$SMOKE/$case" || fail=1
done
[ "$fail" -eq 0 ] && date -Is > "$RUN_ROOT/SMOKE_OK"
echo "end=$(date -Is) fail=$fail"
exit $fail
