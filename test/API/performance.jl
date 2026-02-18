@testitem "performance" begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools

    # Detect eventual gross peformance degradatation

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

    function snaive1(x, cutoff)
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

    function snaive2(x, y, cutoff)
        s = 0.0
        for i in eachindex(x)
            for j in eachindex(y)
                d = norm(y[j] - x[i])
                if d <= cutoff
                    s += d
                end
            end
        end
        return s
    end

    d = (46.36)^3 / 10_000 # water atomic density
    N = 10_000
    V = N * d
    L = V^(1/3)
    cutoff = 12.0 # typical for MD

    x = L .* rand(SVector{3,Float64}, N)
    y = L .* rand(SVector{3,Float64}, N)

    @test s1(f, x, cutoff, false) ≈ snaive1(x, cutoff)
    @test s2(f, x, y, cutoff, false) ≈ snaive2(x, y, cutoff)

    ts = @benchmark s1($f, $x, $cutoff, false) evals=1 samples=3
    tn = @benchmark snaive1($x, $cutoff) evals=1 samples=3
    @test minimum(tn.times) / minimum(ts.times) > 4.2   

    ts = @benchmark s2($f, $x, $y, $cutoff, false) evals=1 samples=3
    tn = @benchmark snaive2($x, $y, $cutoff) evals=1 samples=3
    @test minimum(tn.times) / minimum(ts.times) > 4.8

end