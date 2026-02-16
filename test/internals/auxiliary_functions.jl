@testitem "reduce fallback error" begin
    using CellListMap
    @test_throws ArgumentError CellListMap.reduce("hello", [1, 2, 3])
end

@testitem "pairwise with show_progress" begin
    using CellListMap
    using StaticArrays
    x = rand(SVector{3, Float64}, 100)
    box = CellListMap.Box([1, 1, 1], 0.1)
    const PSP = CellListMap.ParticleSystemPositions
    cl = CellListMap.CellList(PSP(x), box)
    r = CellListMap._pairwise!((pair, r) -> r + pair.d2, 0.0, box, cl; show_progress = true)
    @test r >= 0.0
end

@testitem "get_dim" begin
    using StaticArrays
    using CellListMap: get_dim

    @test get_dim([1,2], [[1,2], [3,4]]) == 2
    @test get_dim([1,2], [[1,2], [3,4]], [[1,2],[3,4]]) == 2
    @test get_dim([1,2,3], [[1,2,3], [3,4,5]]) == 3
    @test get_dim([1,2,3], [[1,2,3], [3,4,5]], [[1,2,3],[3,4,5]]) == 3
    @test get_dim([1,2], Vector{Float64}[]) == 2
    @test get_dim([], Vector{Float64}[[1,2]]) == 2
    @test get_dim(nothing, Vector{Float64}[[1,2]]) == 2
    @test get_dim([1,2], SVector{2,Float64}[[1,2], [3,4]]) == 2

    @test_throws "vectors, or matrices" get_dim([1,2], [])
    @test_throws "Incompatible dimensions" get_dim([1,2], [[1,2,3]])
    @test_throws "2 or 3" get_dim([1,2,3,4], [[1,2,3,4]])
    @test_throws "Could not infer" get_dim(nothing, [[]])

    sys = ParticleSystem(
        xpositions=rand(3,10),
        unitcell=[1,1,1],
        cutoff=0.1,
        output=0.0
    )
    @test get_dim(sys) == 3
    @test get_dim(sys.xpositions) == 3
end

