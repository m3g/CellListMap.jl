@testitem "reduce fallback error" begin
    using CellListMap
    @test_throws ArgumentError CellListMap.reduce("hello", [1, 2, 3])
end

@testitem "pairwise with show_progress" begin
    using CellListMap
    using StaticArrays
    x = rand(SVector{3,Float64}, 100)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(x, box)
    r = pairwise!((pair, r) -> r + pair.d2, 0.0, box, cl; show_progress=true)
    @test r >= 0.0
end