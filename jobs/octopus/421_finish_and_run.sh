#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:40:00
#PBS -N fmesh_421finish
#PBS -j o
#PBS -o logs/421_finish_and_run.pbs.log
#PBS -r n
# Finish a hires refinement's depths and stage an M2 pair against the base.
# Required: FMESH_OUT (a refinement output dir), FMESH_RUN_ROOT.
# Optional: FMESH_HMIN (3), FMESH_HMAX (300), FMESH_RFACTOR (0.2),
#           FMESH_METHOD (equal = TB-FVCOM's smoother; limit = the refinement's).
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
# INVALIDATE first -- this stage's marker and every later one -- before
# anything that can fail, including the environment set-up: a marker left
# from an earlier attempt in a reused run root would let the next stage
# start on this attempt's failure (review 2, R1).
rm -f "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/STAGED" "$FMESH_RUN_ROOT/SMOKE_OK" \
      "$FMESH_RUN_ROOT/base/RUN_OK" "$FMESH_RUN_ROOT/refined/RUN_OK"
. jobs/octopus/common.sh 421_finish_and_run 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
OUTDIR=${FMESH_OUT:?set FMESH_OUT}
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
# NQSV starts this when the refinement ENDS, whatever its outcome; only an
# accepted product may be finished and staged (review F5).
[ -f "$OUTDIR/ACCEPTED" ] || { echo "not accepted: $OUTDIR (no ACCEPTED marker)"; exit 2; }
# exits 3 -- and stops this job -- when the r-factor limit is not reached
python -m fvcom_mesh_tools.cli.finish_depths "$OUTDIR" \
    --hmin "${FMESH_HMIN:-3}" --hmax "${FMESH_HMAX:-300}" \
    --rfactor "${FMESH_RFACTOR:-0.2}" --method "${FMESH_METHOD:-equal}"
FIN=$(ls "$OUTDIR"/fvcom_finished/*_grd.dat)
FIN=${FIN%_grd.dat}
echo "finished case = $FIN"
# The base is the one the refinement started from, as its report names it:
# the same grid, open boundary and depth product. goto2023 was hard-coded
# here, and a refinement of this project's own base (TokyoBayTool) was then
# refused for an open boundary that is not its own. The depth file is the
# one the refinement inherited -- goto2023's grid directory holds two, and
# staging the other would compare bathymetry products as well as meshes.
BASEDIR=$RUN_ROOT/base_case
mkdir -p "$BASEDIR"
read -r B_GRD B_DEP B_OBC < <(python -c "
import json, sys
r = json.load(open(sys.argv[1]))
print(r['base_mesh'], r['base_depth'], r['base_obc'])" "$OUTDIR/report.json")
for f in "$B_GRD" "$B_DEP" "$B_OBC"; do
    [ -f "$f" ] || { echo "the base named in the report is missing: $f"; exit 2; }
done
echo "base = $B_GRD"
cp "$B_GRD" "$BASEDIR/TokyoBayB_grd.dat"
cp "$B_OBC" "$BASEDIR/TokyoBayB_obc.dat"
B_COR=${B_GRD%_grd.dat}_cor.dat
[ -f "$B_COR" ] && cp "$B_COR" "$BASEDIR/TokyoBayB_cor.dat"
cp "$B_DEP" "$BASEDIR/TokyoBayB_dep.dat"
python notebooks/414_refine_m2_prep.py --root "$RUN_ROOT" \
    --base "$BASEDIR/TokyoBayB" --refined "$FIN"
date -Is > "$RUN_ROOT/STAGED"
echo "end=$(date -Is)"
