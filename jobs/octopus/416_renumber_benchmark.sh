#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=120GB
#PBS -l elapstim_req=06:00:00
#PBS -N fmesh_416bench
#PBS -j o
#PBS -o logs/416_renumber_benchmark.pbs.log
#PBS -r n
# Does renumbering make FVCOM faster? Measured, one variable at a time.
#
# The first attempt did not answer the question: the two cases ran
# CONCURRENTLY on one shared node, once each, and the renumbered one came
# back 14 % SLOWER. Two things were wrong with that as an experiment.
#
#   * Concurrency. Two 64-rank jobs on the same node share its memory
#     bandwidth, so each one's time depends on what the other is doing.
#     Here the cases run one after another, in the same job, alternating,
#     so every repeat sees the same machine.
#   * The partition. FVCOM runs METIS itself, on the mesh as given, and
#     METIS is not invariant under renumbering -- so a 64-rank comparison
#     measures the memory layout AND whatever partition METIS happened to
#     choose. A single-rank run has no partition at all (setup_domain.F
#     skips EL_PID entirely for NPROCS = 1), so it measures the layout
#     alone. Both are run, and reported separately, because both questions
#     are real: the second is the claim in the request, the first is what
#     the owner's production configuration would actually experience.
#
# The run is shortened to FMESH_BENCH_DAYS (default 2) -- a timing test does
# not need twenty days of tide -- and repeated FMESH_BENCH_REPEATS times.
#
# Required: FMESH_RUN_ROOT holding the two staged cases (from 414).
# Optional: FMESH_BENCH_RANKS ("1:64"), FMESH_BENCH_DAYS, FMESH_BENCH_REPEATS.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh "416_renumber_benchmark" 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
# Colon-separated, not space-separated: qsub -v splits on whitespace.
RANKS_LIST=${FMESH_BENCH_RANKS:-"1:64"}
RANKS_LIST=${RANKS_LIST//:/ }
DAYS=${FMESH_BENCH_DAYS:-2}
REPEATS=${FMESH_BENCH_REPEATS:-3}
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
for case in base refined; do
    [ -f "$RUN_ROOT/$case/m2_run.nml" ] || { echo "not staged: $RUN_ROOT/$case"; exit 2; }
done

# A shortened copy of each case, so the staged twenty-day run is untouched.
BENCH=$RUN_ROOT/bench
python - "$RUN_ROOT" "$BENCH" "$DAYS" <<'PY'
import re, shutil, sys
from datetime import datetime, timedelta
from pathlib import Path
root, bench, days = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
for case in ("base", "refined"):
    dst = bench / case
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(root / case, dst, ignore=shutil.ignore_patterns("output"))
    (dst / "output").mkdir(exist_ok=True)
    nml = (dst / "m2_run.nml").read_text()
    start = datetime.fromisoformat(
        re.search(r"START_DATE\s*=\s*'([^']+)'", nml).group(1))
    end = (start + timedelta(days=days)).strftime("%Y-%m-%d %H:%M:%S")
    nml = re.sub(r"(END_DATE\s*=\s*')[^']+(')", rf"\g<1>{end}\g<2>", nml)
    nml = nml.replace(str(root / case), str(dst))
    (dst / "m2_run.nml").write_text(nml)
    print(f"[416] {case}: end -> {end}")
PY

set +u
conda deactivate
set -u
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

echo "ranks_list=$RANKS_LIST days=$DAYS repeats=$REPEATS host=$(hostname -s)"
for ranks in $RANKS_LIST; do
    for rep in $(seq 1 "$REPEATS"); do
        # Alternate the order within a repeat so a systematic warm-up or
        # drift in the machine cannot favour one case.
        order="base refined"
        [ $((rep % 2)) -eq 0 ] && order="refined base"
        for case in $order; do
            rm -f "$BENCH/$case/output/"*.nc
            t0=$(date +%s.%N)
            ( cd "$BENCH/$case" && mpiexec -np "$ranks" "$FVCOM" \
                --casename=m2 > "fvcom_${ranks}_${rep}.log" 2>&1 ) || {
                echo "FAILED ranks=$ranks rep=$rep case=$case"
                tail -15 "$BENCH/$case/fvcom_${ranks}_${rep}.log"
                exit 1
            }
            t1=$(date +%s.%N)
            if grep -Ei 'fatal|non[ -]?finite|segmentation|nan detected' \
                    "$BENCH/$case/fvcom_${ranks}_${rep}.log" >/dev/null; then
                echo "UNHEALTHY ranks=$ranks rep=$rep case=$case"; exit 1
            fi
            printf 'BENCH ranks=%s rep=%s case=%-8s seconds=%.2f\n' \
                "$ranks" "$rep" "$case" "$(echo "$t1 - $t0" | bc)"
        done
    done
done
echo "end=$(date -Is)"
