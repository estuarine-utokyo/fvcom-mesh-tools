#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=192GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_387
#PBS -j o
#PBS -o logs/387_edit_ablation_each.pbs.log
#PBS -r n
#============================================================================
# Per-edit ablation: what does EACH of the four human-judgment edits still do?
#
# Four variants, each excluding exactly one edit, run IN PARALLEL in private
# copies of the input tree under scratch (notebooks/325 writes to fixed
# repository-relative paths, so the runs cannot share a working directory).
# The package itself is imported from the editable install, so the copies
# only need data, notebooks and recipes.
#
# Results are collected into outputs/ablation_387/<variant>/ ; the repository
# working tree and the certified outputs are not touched.
#
# Usage (from the repository root): qsub jobs/octopus/387_edit_ablation_each.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 387_edit_ablation_each 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/abl_387.${JOBID}
COLLECT="$REPO/outputs/ablation_387"
VARIANTS="edit_001_haneda_d_runway edit_003_west_edge_crack edit_004_ow05_harbor edit_005_ow05_urayasu"
mkdir -p "$WORK" "$COLLECT"

run_variant() {   # run_variant EXCLUDED_EDIT
    local excl="$1"
    # NOTE: a second assignment in the same `local` statement does not see
    # `excl` yet (bash expands the whole statement first), which under
    # `set -u` aborted every variant in job 115212.
    local dir="$WORK/excl_${excl}"
    mkdir -p "$dir"
    # Only what the scripts read: notebooks, recipes, the OSM land shapefile
    # and the sample-repro inputs (sample_original.14 and friends).
    for d in notebooks recipes; do
        rsync -a --delete "$REPO/$d/" "$dir/$d/"
    done
    mkdir -p "$dir/outputs"
    rsync -a "$REPO/outputs/tb_varres_3r/" "$dir/outputs/tb_varres_3r/"
    rsync -a "$REPO/outputs/sample_repro/" "$dir/outputs/sample_repro/"
    mkdir -p "$dir/outputs/figures"
    (
        cd "$dir"
        export SR_NORMALIZE=on SR_OBC_H1=1680 SR_EDITS_EXCLUDE="$excl"
        export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
        echo "=== variant: excluding $excl ==="
        python notebooks/325_sample_repro.py
        python notebooks/331_finish2.py
        fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 || true
        python notebooks/342_connectivity_check.py || true
        python notebooks/346_one_wide_flags.py || true
        python notebooks/364_over_resolution.py || true
        python notebooks/372_axis_check.py || true
    ) > "$dir/run.log" 2>&1
    echo "variant excl_${excl} finished rc=$? $(date -Is)"
}

echo "=== launching ${VARIANTS} in parallel ==="
pids=()
for v in $VARIANTS; do
    run_variant "$v" &
    pids+=($!)
done
status=0
for p in "${pids[@]}"; do
    wait "$p" || status=1
done
echo "=== all variants done (status=$status) $(date -Is) ==="

for v in $VARIANTS; do
    d="$WORK/excl_${v}"
    mkdir -p "$COLLECT/excl_${v}"
    cp -p "$d/run.log" "$COLLECT/excl_${v}/" 2>/dev/null || true
    for f in sample_repro_final.14 sample_repro_final_qa.json one_wide_cells.json \
             ow_registry.json axis_check.json normalize.json; do
        [ -f "$d/outputs/sample_repro/$f" ] && cp -p "$d/outputs/sample_repro/$f" "$COLLECT/excl_${v}/"
    done
done

echo "=== summary ==="
python notebooks/387_ablation_summary.py "$COLLECT" | tee "$COLLECT/summary.txt"
echo "end=$(date -Is) status=$status"
exit "$status"
