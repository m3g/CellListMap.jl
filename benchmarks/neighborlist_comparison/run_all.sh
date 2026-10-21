#!/bin/bash
# Usage: ./run_all.sh <python> <julia_project>
# All runs pinned to the same 8 physical (P-)cores.
PY=$1; JLPROJ=$2
CORES=0,2,4,6,8,10,12,14
cd "$(dirname "$0")"
mkdir -p results
for t in 1 8; do
    OMP_NUM_THREADS=$t taskset -c $CORES $PY bench_python.py $t results/python_${t}t.json
    taskset -c $CORES julia --project=$JLPROJ -t $t bench_julia.jl results/julia_${t}t.json
done
