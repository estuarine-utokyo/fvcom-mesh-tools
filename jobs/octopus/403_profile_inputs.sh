#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=16
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_403
#PBS -j o
#PBS -o logs/403_profile_inputs.pbs.log
#PBS -r n
# Profile the mesh-generation chain to find where the wall time goes.
# The run log of job 402 shows the INPUT stage alone takes ~1450 s of the
# ~1800 s a variant needs, so the profile targets that stage first.
# Submit from the repository root:
#   qsub jobs/octopus/403_profile_inputs.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 403_profile_inputs 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/profile_403.${JOBID}
mkdir -p "$WORK"
for d in notebooks recipes; do rsync -a --delete "$REPO/$d/" "$WORK/$d/"; done
mkdir -p "$WORK/outputs/figures"
rsync -a "$REPO/outputs/tb_varres_3r/" "$WORK/outputs/tb_varres_3r/"
rsync -a "$REPO/outputs/sample_repro/" "$WORK/outputs/sample_repro/"
rm -f "$WORK/outputs/sample_repro/sample_repro_final.14"

cd "$WORK"
export SR_NORMALIZE=on SR_OBC_H1=1680 FMESH_OVERWRITE=1
export SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=forbid
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
python -m cProfile -o "$WORK/325.prof" notebooks/325_sample_repro.py > "$WORK/325.log" 2>&1

python - "$WORK/325.prof" <<'PY'
import pstats, sys
s = pstats.Stats(sys.argv[1])
print("=== top 40 by cumulative time ===")
s.sort_stats("cumulative").print_stats(40)
print("=== top 40 by total (self) time ===")
s.sort_stats("tottime").print_stats(40)
PY
cp -p "$WORK/325.prof" "$REPO/outputs/325_${JOBID}.prof"
echo "end=$(date -Is)"
