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
. jobs/octopus/common.sh 423_m2_smoke 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
DAYS=${FMESH_DAYS:-2}
RANKS=${FMESH_RANKS:-64}
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
SMOKE=$RUN_ROOT/smoke
[ -f "$RUN_ROOT/STAGED" ] || { echo "not staged: $RUN_ROOT (no STAGED marker)"; exit 2; }
rm -f "$RUN_ROOT/SMOKE_OK"
python - "$RUN_ROOT" "$SMOKE" "$DAYS" <<'PY'
import re, shutil, sys
from datetime import datetime, timedelta
from pathlib import Path
root, smoke, days = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
for case in ("base", "refined"):
    dst = smoke / case
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(root / case, dst, ignore=shutil.ignore_patterns("output"))
    (dst / "output").mkdir(exist_ok=True)
    nml = (dst / "m2_run.nml").read_text()
    start = datetime.fromisoformat(re.search(r"START_DATE\s*=\s*'([^']+)'", nml).group(1))
    end = (start + timedelta(days=days)).strftime("%Y-%m-%d %H:%M:%S")
    nml = re.sub(r"(END_DATE\s*=\s*')[^']+(')", rf"\g<1>{end}\g<2>", nml)
    nml = nml.replace(str(root / case), str(dst))
    (dst / "m2_run.nml").write_text(nml)
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
INSTALLDIR=/octfs/work/G16445/share/local/fvcom/libs/install-oneapi-2025.3.1
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
set +u; conda activate "${FMESH_ENV:-oceanmesh-bench}"; set -u
for case in base refined; do
    python -m fvcom_mesh_tools.cli.check_run "$SMOKE/$case" || fail=1
done
[ "$fail" -eq 0 ] && date -Is > "$RUN_ROOT/SMOKE_OK"
echo "end=$(date -Is) fail=$fail"
exit $fail
