@testitem "Unit cells smaller than the cutoff" begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra: Diagonal, norm
    using Random: MersenneTwister

    struct PairList{T}
        v::Vector{T}
    end
    CellListMap.copy_output(p::PairList) = PairList(copy(p.v))
    CellListMap.reset_output!(p::PairList) = (empty!(p.v); p)
    CellListMap.reducer(a::PairList, b::PairList) = (append!(a.v, b.v); a)

    # Representation of a pair which is the same for (i, j, S) and (j, i, -S)
    canonical(i, j, S) = (i < j || (i == j && S > zero(S))) ? (i, j, S) : (j, i, -S)

    function celllist_pairs(x, y, unitcell, cutoff, parallel, lcell)
        N = length(x[1])
        out = PairList(Tuple{Int, Int, SVector{N, Int}, Float64}[])
        kw = (; cutoff, unitcell, output = out, parallel, lcell)
        sys = isnothing(y) ? ParticleSystem(; xpositions = x, kw...) : ParticleSystem(; xpositions = x, ypositions = y, kw...)
        pairs = pairwise!((p, l) -> (push!(l.v, (p.i, p.j, p.shift, p.d)); l), sys).v
        return isnothing(y) ? sort([(canonical(p[1], p[2], p[3])..., p[4]) for p in pairs]) : sort(pairs)
    end

    # Brute force over periodic images: shifts around the one that brings j close to i
    function brute_force_pairs(x, y, M, cutoff)
        N = length(x[1])
        n = 2 + ceil(Int, cutoff / minimum(CellListMap._face_distances(SMatrix{N, N}(M))))
        out = Tuple{Int, Int, SVector{N, Int}, Float64}[]
        yy = isnothing(y) ? x : y
        for i in eachindex(x), j in eachindex(yy)
            isnothing(y) && j < i && continue
            S0 = round.(Int, M \ (x[i] - yy[j]))
            for S in Iterators.product(ntuple(_ -> -n:n, N)...)
                Sv = S0 + SVector{N, Int}(S)
                isnothing(y) && j == i && !(Sv > zero(Sv)) && continue
                d = norm(yy[j] + M * Sv - x[i])
                d <= cutoff && push!(out, (i, j, Sv, d))
            end
        end
        return sort(out)
    end

    rng = MersenneTwister(7)
    for N in (2, 3), kind in (:ortho, :tric), cross in (false, true), parallel in (false, true), _ in 1:10
        cutoff = 0.5 + 4 * rand(rng)
        lcell = rand(rng, 1:3)
        if kind == :ortho
            sides = 0.4 .+ 4 .* rand(rng, N)
            unitcell = sides
            M = Matrix(Diagonal(sides))
        else
            M = Matrix(Diagonal(0.6 .+ 4 .* rand(rng, N))) .+ 2 .* (rand(rng, N, N) .- 0.5)
            unitcell = M
        end
        x = [SVector{N}(-5 .+ 10 .* rand(rng, N)) for _ in 1:8]
        y = cross ? [SVector{N}(-5 .+ 10 .* rand(rng, N)) for _ in 1:6] : nothing
        cl = celllist_pairs(x, y, unitcell, cutoff, parallel, lcell)
        bf = brute_force_pairs(x, y, M, cutoff)
        @test length(cl) == length(bf)
        @test all(p -> p[1][1:3] == p[2][1:3] && isapprox(p[1][4], p[2][4]; atol = 1.0e-8), zip(cl, bf))
    end

    # Single particle within the cutoff of its own images
    nb = neighborlist(; positions = [SVector(0.5, 0.5, 0.5)], cutoff = 1.5, unitcell = [1.0, 1.0, 1.0])
    @test length(nb) == 9
    @test count(p -> p[3] ≈ 1.0, nb) == 3
    @test count(p -> p[3] ≈ sqrt(2), nb) == 6

    # Degenerate unit cells are rejected
    @test_throws ArgumentError ParticleSystem(; positions = rand(SVector{3, Float64}, 10), cutoff = 0.1, unitcell = [1.0 2.0 0.0; 1.0 2.0 0.0; 0.0 0.0 1.0], output = 0.0)
end
