#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=48
#PBS -l memsz_job=144GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_408
#PBS -j o
#PBS -o logs/408_natural.pbs.log
#PBS -r n
# Compare the three channel-width policies under one coastline fit.
# forbid  pushes banks into land until two rows fit (certified)
# natural carves kept port/dead-end records at their natural width
# allow   also lowers the keep bar (known to fail QA; kept as the control)
# This run also measures the detect_waterways speed-up: the input stage
# was 1,427 s of 1,487 s wall before the loop-invariant buffer was hoisted.
# Submit from the repository root:
#   qsub jobs/octopus/408_natural.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 408_natural 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/natural_408.${JOBID}
COLLECT="$REPO/outputs/natural_408.${JOBID}"
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
run_variant forbid  "SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=forbid SR_COAST_FIT=on" & pids+=($!)
run_variant natural "SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=natural SR_COAST_FIT=on" & pids+=($!)
run_variant allow   "SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=allow SR_COAST_FIT=on" & pids+=($!)
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
echo "=== variants done (status=$status) $(date -Is) ==="

for v in forbid natural allow; do
    mkdir -p "$COLLECT/$v"
    cp -p "$WORK/$v/run.log" "$COLLECT/$v/" 2>/dev/null || true
    cp -p "$WORK/$v"/outputs/figures/404_*.png "$COLLECT/$v/" 2>/dev/null || true
    for f in sample_repro_final.14 sample_repro_final_qa.json one_wide_cells.json \
             over_resolution.json coverage_gaps.json waterways.json \
             channel_policy.json coast_fit.json land_breaches.json; do
        [ -f "$WORK/$v/outputs/sample_repro/$f" ] && cp -p "$WORK/$v/outputs/sample_repro/$f" "$COLLECT/$v/"
    done
done
echo "=== summary ==="
python notebooks/392_sizing_summary.py "$COLLECT" | tee "$COLLECT/summary.txt"
python notebooks/393_sizing_compare.py "$COLLECT" || true
echo "end=$(date -Is) status=$status"
exit "$status"
