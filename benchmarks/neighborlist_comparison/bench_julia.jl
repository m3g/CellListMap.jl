# Neighbor list benchmarks of Julia packages (periodic cubic boxes).
# Usage: julia -t <nthreads> bench_julia.jl <output.json> [package,package,...]
# Results are merged into <output.json> after each package.
using CellListMap, NearestNeighbors, PointNeighbors, StaticArrays, JSON3, Random, Pkg

const cfg = JSON3.read(read(joinpath(@__DIR__, "common_sizes.json"), String))
const cutoff = Float64(cfg.cutoff)
const TMAX = 60.0
const parallel = Threads.nthreads() > 1

pkgversion_str(m) = string(pkgversion(m))

# minimum-image distance
@inline function pdist(a, b, L)
    d = a - b
    d = d - L * round.(d / L)
    return sqrt(sum(abs2, d))
end

function run_celllistmap(x, L)
    nb = neighborlist(; xpositions=x, cutoff, unitcell=SVector(L, L, L), parallel)
    return length(nb)
end

function nn_chunk(kdtree, x, L, irange)
    ptree = PeriodicTree(kdtree, SVector(0.0, 0.0, 0.0), SVector(L, L, L))
    out = Tuple{Int,Int,Float64}[]
    idxs = Int[]
    for i in irange
        empty!(idxs)
        inrange!(idxs, ptree, x[i], cutoff)
        for j in idxs
            j > i && push!(out, (i, j, pdist(x[i], x[j], L)))
        end
    end
    return out
end

function run_nearestneighbors(x, L)
    kdtree = KDTree(x)
    if !parallel
        return length(nn_chunk(kdtree, x, L, eachindex(x)))
    end
    chunks = Iterators.partition(eachindex(x), cld(length(x), 8 * Threads.nthreads()))
    tasks = [Threads.@spawn nn_chunk(kdtree, x, L, c) for c in chunks]
    return length(reduce(vcat, fetch.(tasks)))
end

function run_pointneighbors(x, L)
    xm = reinterpret(reshape, Float64, x)
    # FullGridCellList is much faster than the default DictionaryCellList for a known domain
    cell_list = FullGridCellList(; min_corner=[0.0, 0.0, 0.0], max_corner=[L, L, L], search_radius=cutoff,
        max_points_per_cell=ceil(Int, 2 * length(x) * (cutoff / L)^3) + 64)
    nhs = GridNeighborhoodSearch{3}(; search_radius=cutoff, n_points=length(x), cell_list,
        periodic_box=PeriodicBox(min_corner=[0.0, 0.0, 0.0], max_corner=[L, L, L]))
    backend = parallel ? ThreadsStaticBackend() : SerialBackend()
    initialize!(nhs, xm, xm; parallelization_backend=backend)
    lists = [Tuple{Int,Int,Float64}[] for _ in 1:Threads.maxthreadid()]
    foreach_point_neighbor(xm, xm, nhs; parallelization_backend=backend) do i, j, _, d
        i < j && push!(lists[Threads.threadid()], (i, j, d))
    end
    return length(reduce(vcat, lists))
end

const competitors = [
    "CellListMap" => (run_celllistmap, pkgversion_str(CellListMap)),
    "NearestNeighbors" => (run_nearestneighbors, pkgversion_str(NearestNeighbors)),
    "PointNeighbors" => (run_pointneighbors, pkgversion_str(PointNeighbors)),
]

function timeit(f, x, L)
    best = Inf; n = 0
    for _ in 1:3
        GC.gc()
        t = @elapsed (n = f(x, L))
        best = min(best, t)
        t > 5.0 && break
    end
    return best, n
end

function main(outfile, selected)
    rng = Xoshiro(321)
    results = isfile(outfile) ? JSON3.read(read(outfile, String), Dict{String,Any}) : Dict{String,Any}()
    for (name, (f, ver)) in competitors
        name in selected || continue
        f(30.0 .* rand(rng, SVector{3,Float64}, 3000), 30.0) # compile (box side must exceed 2*cutoff)
        res = Dict{String,Any}("version" => ver, "cd" => [], "cv" => [])
        for case in ("cd", "cv")
            Ns = cfg[case*"_N"]
            for (k, N) in enumerate(Ns)
                L = case == "cd" ? (N / cfg.density)^(1 / 3) : Float64(cfg.cv_side)
                x = L .* rand(rng, SVector{3,Float64}, N)
                t, npairs = timeit(f, x, L)
                println("$name $(Threads.nthreads())t $case N=$N t=$(round(t; digits=3)) npairs=$npairs"); flush(stdout)
                push!(res[case], [N, t, npairs])
                p = case == "cd" ? 1 : 2
                k < length(Ns) && t * (Ns[k+1] / N)^p > TMAX && break
            end
        end
        results[name] = res
        open(outfile, "w") do io
            JSON3.pretty(io, results)
        end
    end
end

main(ARGS[1], length(ARGS) > 1 ? split(ARGS[2], ",") : first.(competitors))
