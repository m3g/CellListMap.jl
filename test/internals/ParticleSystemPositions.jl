@testitem "ParticleSystemPositions interface" begin
    using StaticArrays
    using CellListMap
    using ShowMethodTesting
    const PSP = CellListMap.ParticleSystemPositions

    # Construction from vector of SVectors
    x_svecs = rand(SVector{3,Float64}, 10)
    p = PSP(x_svecs)
    @test p isa PSP{3,Float64}
    @test length(p) == 10
    @test p.updated[] == true
    @test p[1] == x_svecs[1]

    # getindex
    @test p[1] isa SVector{3,Float64}

    # setindex! flags updated
    p.updated[] = false
    p[1] = SVector(1.0, 2.0, 3.0)
    @test p[1] == SVector(1.0, 2.0, 3.0)
    @test p.updated[] == true

    # length, size, ndims
    @test length(p) == 10
    @test size(p) == (10,)
    @test ndims(p) == 1
    @test ndims(typeof(p)) == 1

    # axes, keys, eachindex
    @test axes(p) == (Base.OneTo(10),)
    @test keys(p) == 1:10
    @test collect(eachindex(p)) == collect(1:10)

    # firstindex, lastindex, first, last
    @test firstindex(p) == 1
    @test lastindex(p) == 10
    @test first(p) == p[1]
    @test last(p) == p[10]

    # iterate
    let count = 0
        for v in p
            count += 1
            @test v isa SVector{3,Float64}
        end
        @test count == length(p)
    end

    # copy preserves data and updated flag
    p.updated[] = false
    pc = copy(p)
    @test length(pc) == length(p)
    @test pc[1] == p[1]
    @test pc.updated[] == false
    # copy is independent
    pc[1] = SVector(9.0, 9.0, 9.0)
    @test p[1] != pc[1]

    # similar
    ps = similar(p)
    @test length(ps) == length(p)
    @test ps.updated[] == true

    # push! flags updated and adds element
    p5 = PSP(rand(SVector{3,Float64}, 3))
    p5.updated[] = false
    push!(p5, SVector(4.0, 5.0, 6.0))
    @test length(p5) == 4
    @test p5[4] == SVector(4.0, 5.0, 6.0)
    @test p5.updated[] == true

    # append! flags updated and adds elements
    p6 = PSP(rand(SVector{3,Float64}, 3))
    p6.updated[] = false
    append!(p6, [SVector(1.0, 1.0, 1.0), SVector(2.0, 2.0, 2.0)])
    @test length(p6) == 5
    @test p6[4] == SVector(1.0, 1.0, 1.0)
    @test p6[5] == SVector(2.0, 2.0, 2.0)
    @test p6.updated[] == true

    # resize! flags updated and changes length
    p7 = PSP(rand(SVector{3,Float64}, 5))
    p7.updated[] = false
    resize!(p7, 8)
    @test length(p7) == 8
    @test p7.updated[] == true

    # empty! flags updated and clears
    p8 = PSP(rand(SVector{3,Float64}, 5))
    p8.updated[] = false
    empty!(p8)
    @test length(p8) == 0
    @test p8.updated[] == true

    # In-place broadcasting flags updated
    p9 = PSP(rand(SVector{3,Float64}, 5))
    p9.updated[] = false
    new_vals = [SVector(0.0, 0.0, 0.0) for _ in 1:5]
    p9 .= new_vals
    @test p9.updated[] == true
    @test all(p9[i] == SVector(0.0, 0.0, 0.0) for i in 1:5)

    # Out-of-place broadcasting returns plain Array
    p10 = PSP(rand(SVector{3,Float64}, 5))
    result = p10 .+ Ref(SVector(1.0, 1.0, 1.0))
    @test result isa Vector{SVector{3,Float64}}
    @test result[1] ≈ p10[1] + SVector(1.0, 1.0, 1.0)

    # copyto! from another ParticleSystemPositions flags updated
    p11 = PSP(rand(SVector{3,Float64}, 5))
    p12 = PSP(rand(SVector{3,Float64}, 5))
    p11.updated[] = false
    copyto!(p11, p12)
    @test p11.updated[] == true
    @test all(p11[i] == p12[i] for i in 1:5)

    # view shares the updated reference flag
    p13 = PSP(rand(SVector{3,Float64}, 10))
    p13.updated[] = false
    v = @view p13[1:5]
    v[1] = SVector(1.0, 1.0, 1.0)
    @test length(v) == 5
    @test p13[1] == v[1]
    @test p13.updated[] = true

    # show does not error
    @test isapprox(parse_show(p; repl=Dict("StaticArraysCore." => "")), """
      ParticleSystemPositions, updated: false, with 10-element Vector{SVector{3, Float64}}:
        [1.0, 2.0, 3.0]
        [0.939318968863584, 0.5097640417408886, 0.19846095185478763]
        [0.5368806196615254, 0.8508139226462262, 0.4786069714983523]
        [0.4127933306918298, 0.05227010510480201, 0.04247943213111527]
    """; float_match = (x,y) -> true)

    # error on construction from non-SVector element type
    @test_throws MethodError PSP(Vector{Float64}[])
    @test_throws MethodError PSP([rand(3) for _ in 1:5])

end

@testitem "ParticleSystemPositions integration with ParticleSystem" begin
    using StaticArrays
    using CellListMap

    # Flag is reset after pairwise! for ParticleSystem1
    x = rand(SVector{3,Float64}, 100)
    sys = ParticleSystem(positions=x, cutoff=0.1, unitcell=[1, 1, 1], output=0.0)
    @test sys.xpositions.updated[] == false
    sys.xpositions[1] = rand(SVector{3,Float64})
    @test sys.xpositions.updated[] == true
    pairwise!((pair, r) -> r += pair.d2, sys)
    @test sys.xpositions.updated[] == false

    # Flag is reset after pairwise! for ParticleSystem2
    x = rand(SVector{3,Float64}, 100)
    y = rand(SVector{3,Float64}, 50)
    sys2 = ParticleSystem(xpositions=x, ypositions=y, cutoff=0.1, unitcell=[1, 1, 1], output=0.0)
    @test sys2.xpositions.updated[] == false
    @test sys2.ypositions.updated[] == false
    # Modifying only y triggers only y flag
    sys2.ypositions[1] = SVector(0.5, 0.5, 0.5)
    @test sys2.xpositions.updated[] == false
    @test sys2.ypositions.updated[] == true
    pairwise!((pair, r) -> r += pair.d2, sys2)
    @test sys2.xpositions.updated[] == false
    @test sys2.ypositions.updated[] == false

    # update! with unitcell sets updated flag
    sys.xpositions.updated[] = false
    update!(sys; unitcell=[2, 2, 2])
    @test sys.xpositions.updated[] == true

    # update! with cutoff sets updated flag
    sys.xpositions.updated[] = false
    update!(sys; cutoff=0.2)
    @test sys.xpositions.updated[] == true

    # update! with cutoff sets both flags for ParticleSystem2
    sys2.xpositions.updated[] = false
    sys2.ypositions.updated[] = false
    update!(sys2; cutoff=0.2)
    @test sys2.xpositions.updated[] == true
    @test sys2.ypositions.updated[] == true

    # push! works through ParticleSystem and triggers recomputation
    x = rand(SVector{3,Float64}, 100)
    sys = ParticleSystem(positions=x, cutoff=0.1, unitcell=[1, 1, 1], output=0.0)
    pairwise!((pair, r) -> r += pair.d2, sys)
    @test sys.xpositions.updated[] == false
    push!(sys.xpositions, rand(SVector{3,Float64}))
    @test length(sys.xpositions) == 101
    @test sys.xpositions.updated[] == true
end
