# Neighbor list comparison

Benchmarks of `CellListMap.neighborlist` against NearestNeighbors.jl, PointNeighbors.jl,
scipy, vesin, freud and matscipy, used for the plots in the main README.

- `common_sizes.json`: density, cutoff, box side and system sizes for both cases.
- `bench_julia.jl`, `bench_python.py`: time each package. The results are merged into
  `results/<lang>_<n>t.json` after each package; an optional third argument restricts the packages to run.
  For each package and case, larger systems are skipped once the predicted time of the next size
  exceeds 60 s.
- `run_all.sh <python> <julia_project>`: runs everything, pinned to 8 physical cores
  (`taskset -c 0,2,...,14`; adjust `CORES` to your machine).
- `plot.jl`: writes `docs/src/assets/b_nl_cd.png` and `docs/src/assets/b_nl_cv.png`.

Environments used:

```bash
python -m venv pyenv
./pyenv/bin/pip install numpy scipy freud-analysis vesin matscipy
julia --project=jlenv -e 'using Pkg; Pkg.develop(path="../.."); Pkg.add(["NearestNeighbors","PointNeighbors","StaticArrays","JSON3","Plots"])'
julia --project=jlenv plot.jl
```
