# usage: julia -t 8 --project prototype/bench.jl [proto|clm] N
include(joinpath(@__DIR__, "FracCellList.jl"))
using .FracCellList, CellListMap, StaticArrays, Random
which, N = ARGS[1], parse(Int, ARGS[2])
L = (N / 0.1)^(1 / 3)
x = [SVector{3,Float64}(L .* rand(MersenneTwister(1 + i), 3)) for i in 1:N]
nrep = N > 500_000 ? 5 : 20
best(f) = (f(); minimum(@elapsed(f()) for _ in 1:nrep) * 1e3)
r(t) = lpad(round(t, digits = 2), 8)
for cell in (:cubic, :tric), par in (false, true)
    M = cell == :cubic ? [L 0 0; 0 L 0; 0 0 L] : [L 0.3L 0.2L; 0 L 0.1L; 0 0 L]
    if which == "clm"
        uc = cell == :cubic ? [L, L, L] : M
        px = CellListMap.ParticleSystemPositions(x); box = CellListMap.Box(uc, 5.0)
        cl = CellListMap.CellList(px, box; parallel = par); aux = par ? CellListMap.AuxThreaded(cl) : nothing
        tb = best(() -> (px.updated[] = true; CellListMap.UpdateCellList!(px, box, cl, aux; parallel = par)))
        sys = ParticleSystem(; xpositions = x, cutoff = 5.0, unitcell = uc, output = 0, parallel = par)
        tm = best(() -> pairwise!((p, n) -> n + 1, sys))
        npairs = pairwise!((p, n) -> n + 1, sys)
        tn = best(() -> neighborlist(; xpositions = x, cutoff = 5.0, unitcell = uc, parallel = par))
        println(rpad("clm $cell par=$par", 26), "build", r(tb), "  map", r(tm), "  neighborlist", r(tn), "  pairs ", npairs)
    else
        for k in (1,)
            grid = FracCellList.Grid(SMatrix{3,3}(M), 5.0; cells_per_cutoff = k)
            cl = FracCellList.CellList(x, grid; parallel = par)
            tb = best(() -> FracCellList.build!(cl, x; parallel = par))
            count(i, j, xi, xj, d2, S, n) = n + 1
            tm = best(() -> FracCellList.map_pairs(count, 0, cl; parallel = par))
            npairs = FracCellList.map_pairs(count, 0, cl; parallel = par)
            tn = best(() -> FracCellList.neighborlist(FracCellList.build!(cl, x; parallel = par); parallel = par))
            println(rpad("proto k=$k $cell par=$par", 26), "build", r(tb), "  map", r(tm), "  neighborlist", r(tn), "  pairs ", npairs, "  offsets ", length(grid.offsets))
        end
    end
    GC.gc()
end
