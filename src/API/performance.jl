@testitem "performance" begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra

    d = (46.36)^3 / 10_000 # water atomic density
    N = 10_000
    V = N * d
    L = V^(1/3)
    cutoff = 12.0 # typical for MD

    x = L .* rand(SVector{3,Float64}, N)
    y = L .* rand(SVector{3,Float64}, N)

    f(pair, out) = out += pair.d

    function s1(f, x, cutoff, parallel) 
        sys = ParticleSystem(
            xpositions=x,
            cutoff=cutoff,
            output=0.0,
            parallel=parallel,
        )
        pairwise!(f, sys)
    end

    function s2(f, x, y, cutoff, parallel) 
        sys = ParticleSystem(
            xpositions=x,
            ypositions=y,
            cutoff=cutoff,
            output=0.0,
            parallel=parallel,
        )
        pairwise!(f, sys)
    end

    function snaive(x, cutoff)
        s = 0.0
        for i in eachindex(x)
            for j in i+1:length(x)
                d = norm(x[j] - x[i])
                if d <= cutoff
                    s += d
                end
            end
        end
        return s
    end






end