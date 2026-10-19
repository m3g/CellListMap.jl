@testitem "NeighborPair shifts" begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra: Diagonal, norm
    using Random: MersenneTwister

    # The shift satisfies y[j] + unitcell * shift - x[i] == pair.y - pair.x, for
    # coordinates spread over several periodic images (including negative ones).
    rng = MersenneTwister(3)
    for T in (Float64, Float32), N in (2, 3), cellkind in (:ortho, :tric, :none), cross in (false, true), parallel in (false, true)
        L = T(20)
        unitcell = cellkind == :ortho ? fill(L, N) :
            cellkind == :tric ? (N == 3 ? T[L 3 -2; 0 L 4; 0 0 L] : T[L 5; -3 L]) : nothing
        scale = cellkind == :none ? T(1 / 4) : T(1)
        x = [scale * SVector{N, T}(T(-60) .+ T(150) .* rand(rng, T, N)) for _ in 1:2000]
        y = [scale * SVector{N, T}(T(-60) .+ T(150) .* rand(rng, T, N)) for _ in 1:1000]
        yy = cross ? y : x
        M = unitcell isa AbstractVector ? Diagonal(unitcell) : unitcell
        kw = (; cutoff = T(4), unitcell, parallel)
        sys_err = cross ? ParticleSystem(; xpositions = x, ypositions = y, output = 0.0, kw...) :
            ParticleSystem(; xpositions = x, output = 0.0, kw...)
        err = pairwise!(sys_err) do p, e
            r = isnothing(M) ? yy[p.j] - x[p.i] : yy[p.j] + M * p.shift - x[p.i]
            e + norm(r - (p.y - p.x))
        end
        sys_count = cross ? ParticleSystem(; xpositions = x, ypositions = y, output = 0, kw...) :
            ParticleSystem(; xpositions = x, output = 0, kw...)
        n_shifted = pairwise!((p, n) -> n + any(!iszero, p.shift), sys_count)
        n_pairs = pairwise!((p, n) -> n + 1, sys_count)
        @test err / n_pairs < (T == Float64 ? 1.0e-10 : 1.0e-4)
        if cellkind == :none
            @test n_shifted == 0
        else
            @test n_shifted > 0
        end
    end

    # The (i, j, shift) set equals the one obtained by brute force over periodic images
    x = [SVector{3}(-30 .+ 80 .* rand(rng, 3)) for _ in 1:300]
    M = [20.0 3 -2; 0 20 4; 0 0 20]
    sys = ParticleSystem(; xpositions = x, cutoff = 4.0, unitcell = M, output = Tuple{Int, Int, SVector{3, Int}}[], parallel = false)
    pairwise!(sys) do p, list
        push!(list, p.i < p.j ? (p.i, p.j, p.shift) : (p.j, p.i, -p.shift))
    end
    brute = Tuple{Int, Int, SVector{3, Int}}[]
    for i in eachindex(x), j in (i + 1):length(x), S in Iterators.product(-6:6, -6:6, -6:6)
        norm(x[j] + M * SVector(S) - x[i]) <= 4.0 && push!(brute, (i, j, SVector(S)))
    end
    @test length(brute) > 0
    @test sort(sys.output) == sort(brute)

    # Coordinates far from the unit cell
    x = [SVector(1.0e6 + 1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)]
    sys = ParticleSystem(; xpositions = x, cutoff = 4.0, unitcell = [10.0, 10.0, 10.0], output = Tuple{Int, Int, SVector{3, Int}}[], parallel = false)
    pairwise!((p, list) -> push!(list, (p.i, p.j, p.shift)), sys)
    @test length(sys.output) == 1
    i, j, S = only(sys.output)
    @test norm(x[j] + 10.0 * S - x[i]) ≈ 1.0
end
