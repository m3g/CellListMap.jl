include(joinpath(@__DIR__, "FracCellList.jl"))
using .FracCellList, StaticArrays, LinearAlgebra, Random, Test
canonical(i, j, S) = (i < j || (i == j && S > zero(S))) ? (i, j, S) : (j, i, -S)
function proto_pairs(x, M, rc, parallel, k)
    N = length(x[1])
    grid = FracCellList.Grid(SMatrix{N,N}(M), rc; cells_per_cutoff = k)
    cl = FracCellList.CellList(x, grid; parallel)
    T = Tuple{Int,Int,SVector{N,Int},Float64}
    out = FracCellList.map_pairs((i, j, xi, xj, d2, S, l) -> begin
            # consistency of the image coordinates with the shift
            @assert norm((x[j] + M * S - x[i]) - (xj - xi)) < 1e-8
            push!(l, (i, j, S, sqrt(d2))); l
        end, T[], cl; parallel, reducer = append!, copy_output = l -> similar(l, 0))
    sort([(canonical(o[1], o[2], o[3])..., o[4]) for o in out])
end
function brute(x, M, rc)
    N = length(x[1]); n = 2 + ceil(Int, rc / minimum(FracCellList.face_distances(SMatrix{N,N}(M))))
    out = Tuple{Int,Int,SVector{N,Int},Float64}[]
    for i in eachindex(x), j in i:length(x)
        S0 = round.(Int, M \ (x[i] - x[j]))
        for S in Iterators.product(ntuple(_ -> -n:n, N)...)
            Sv = S0 + SVector{N,Int}(S)
            j == i && !(Sv > zero(Sv)) && continue
            d = norm(x[j] + M * Sv - x[i])
            d <= rc && push!(out, (i, j, Sv, d))
        end
    end
    sort(out)
end
rng = MersenneTwister(11)
@testset "fractional grid vs brute force" begin
for N in (2, 3), kind in (:ortho, :tric), size in (:small, :large), parallel in (false, true), _ in 1:15
    rc = 0.5 + 2 * rand(rng)
    scale = size == :small ? 1.0 : 6.0
    M = Matrix(Diagonal(scale .* (0.5 .+ 4 .* rand(rng, N))))
    kind == :tric && (M .+= scale .* 1.5 .* (rand(rng, N, N) .- 0.5))
    x = [SVector{N}(scale .* (-5 .+ 10 .* rand(rng, N))) for _ in 1:(size == :small ? 10 : 60)]
    k = rand(rng, 1:3)
    a = proto_pairs(x, M, rc, parallel, k); b = brute(x, M, rc)
    ok = length(a) == length(b) && all(p -> p[1][1:3] == p[2][1:3] && isapprox(p[1][4], p[2][4]; atol = 1e-8), zip(a, b))
    ok || println("FAIL N=$N $kind $size par=$parallel k=$k rc=$rc: $(length(a)) vs $(length(b))")
    @test ok
end
end
# neighborlist: same pairs as map_pairs, serial and parallel
x = [SVector{3}(20 .* rand(rng, 3)) for _ in 1:3000]; M = SMatrix{3,3}([20.0 3 -2; 0 20 4; 0 0 20])
cl = FracCellList.CellList(x, FracCellList.Grid(M, 2.5))
ref = sort(FracCellList.map_pairs((i, j, xi, xj, d2, S, l) -> (push!(l, (min(i,j), max(i,j))); l), Tuple{Int,Int}[], cl; parallel = false))
for par in (false, true)
    nb = FracCellList.neighborlist(cl; parallel = par)
    @test sort([(min(p[1],p[2]), max(p[1],p[2])) for p in nb]) == ref
end
println("neighborlist ok: ", length(ref), " pairs")
