#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=48
#PBS -l memsz_job=144GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_409
#PBS -j o
#PBS -o logs/409_verify.pbs.log
#PBS -r n
# Final verification of this round: the certified chain with the coastline
# fit off and on, under the repaired 364 width measure, the two-measure
# time-step guard and the hoisted land buffer.  Expected: both variants
# pass QA and the 364 gate, the fit closes the shoreline offset, neither
# time step moves, and the input stage runs in ~250 s instead of ~1,430 s.
# Submit from the repository root:
#   qsub jobs/octopus/409_verify.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 409_verify 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/verify_409.${JOBID}
COLLECT="$REPO/outputs/verify_409.${JOBID}"
mkdir -p "$WORK" "$COLLECT"

run_variant() {
    local name="$1"
    local envs="$2"
    local dir="$WORK/$name"
    mkdir -p "$dir"
    for d in notebooks recipes; do
        rsync -a --delete "$REPO/$d/" "$dir/$d/"
    done
    mkdir -p "$dir/outputs/figures"
    rsync -a "$REPO/outputs/tb_varres_3r/" "$dir/outputs/tb_varres_3r/"
    rsync -a "$REPO/outputs/sample_repro/" "$dir/outputs/sample_repro/"
    # These are fresh task-owned copies; avoid reporting stale results if
    # generation or finishing fails.
    rm -f "$dir/outputs/sample_repro/sample_repro_final.14" \
          "$dir/outputs/sample_repro/sample_repro_final_qa.json" \
          "$dir/outputs/sample_repro/kept_channel_paths.json"
    (
        cd "$dir"
        export SR_NORMALIZE=on SR_OBC_H1=1680 FMESH_OVERWRITE=1
        export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
        # shellcheck disable=SC2163
        for kv in $envs; do export "$kv"; done
        echo "=== variant $name ($envs) ==="
        _t0=$(date +%s)
        python notebooks/325_sample_repro.py
        echo "[job] 325 wall = $(( $(date +%s) - _t0 )) s"
        python notebooks/331_finish2.py
        fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 || true
        python notebooks/342_connectivity_check.py || true
        python notebooks/346_one_wide_flags.py || true
        python notebooks/364_over_resolution.py || true
        python notebooks/389_coverage_gaps.py || true
        python notebooks/394_spec_match.py outputs/sample_repro/sample_repro_final.14 "$name" || true
        python notebooks/404_coast_fit.py outputs/sample_repro/sample_repro_final.14 \
            outputs/sample_repro || true
    ) > "$dir/run.log" 2>&1
    echo "variant $name finished rc=$? $(date -Is)"
}

pids=()
run_variant nofit "SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=forbid SR_COAST_FIT=off" & pids+=($!)
run_variant fit   "SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=forbid SR_COAST_FIT=on" & pids+=($!)
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
echo "=== variants done (status=$status) $(date -Is) ==="

for v in nofit fit; do
    mkdir -p "$COLLECT/$v"
    cp -p "$WORK/$v/run.log" "$COLLECT/$v/" 2>/dev/null || true
    cp -p "$WORK/$v"/outputs/figures/404_*.png "$COLLECT/$v/" 2>/dev/null || true
    for f in sample_repro_final.14 sample_repro_final_qa.json one_wide_cells.json \
             over_resolution.json coverage_gaps.json waterways.json \
             channel_policy.json coast_fit.json coast_offsets.json land_breaches.json; do
        [ -f "$WORK/$v/outputs/sample_repro/$f" ] && cp -p "$WORK/$v/outputs/sample_repro/$f" "$COLLECT/$v/"
    done
done
echo "=== summary ==="
python notebooks/392_sizing_summary.py "$COLLECT" | tee "$COLLECT/summary.txt"
python notebooks/393_sizing_compare.py "$COLLECT" || true
echo "end=$(date -Is) status=$status"
exit "$status"
