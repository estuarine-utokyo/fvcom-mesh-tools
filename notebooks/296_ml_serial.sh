#!/bin/bash
# serialize MATLAB jobs after Ex7 finishes (license limit -87)
set -u
until [[ -z $(pjstat 2>/dev/null | awk '$1==6180046') ]]; do sleep 300; done
cd "${HOME}/Github/fvcom-mesh-tools"
J5=$(pjsub notebooks/294_matlab_ex5.pjsub 2>&1 | grep -oE "[0-9]{7}")
echo "ML5 resubmitted $J5"
until [[ -z $(pjstat 2>/dev/null | awk -v j="$J5" '$1==j') ]]; do sleep 180; done
grep -E "GOLDEN|Error" "$(ls -t fvcm-mlx5.*.out | head -1)" | tail -1
