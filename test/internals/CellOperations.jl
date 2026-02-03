@testitem "wrap_to_first" begin
    x = [15.0, 13.0]
    unit_cell_matrix = [10.0 0.0; 0.0 10.0]
    @test CellListMap.wrap_to_first(x, unit_cell_matrix) ≈ [5.0, 3.0]
end

@testitem "wrap_relative_to" begin
    using StaticArrays
    using Unitful
    for (x, y, xy, yx) in [
            ([15.0, 13.0], [4.0, 2.0], [5.0, 3.0], [14.0, 12.0]),
            ([-7.0, -6.0], [1.0, 2.0], [3.0, 4.0], [-9.0, -8.0]),
        ]
        # triclinic cells: unitcell is a matrix
        unit_cell_matrix = SMatrix{2, 2}(10.0, 0.0, 0.0, 10.0)
        @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, unit_cell_matrix) ≈ yx
        unit_cell_matrix = SMatrix{2, 2}(-10.0, 0.0, 0.0, -10.0)
        @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, unit_cell_matrix) ≈ yx
        # orthorhombic cells: sides is a vector
        sides = [10.0, 10.0]
        @test CellListMap.wrap_relative_to(x, y, sides) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, sides) ≈ yx
        sides = [-10.0, -10.0]
        @test CellListMap.wrap_relative_to(x, y, sides) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, sides) ≈ yx
    end
    # test unit propagations
    x = [15.0, 13.0]u"nm"
    y = [4.0, 2.0]u"nm"
    unit_cell_matrix = SMatrix{2, 2}(10.0, 0.0, 0.0, 10.0)u"nm"
    @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ [5.0, 3.0]u"nm"
    unit_cell_matrix = [10.0, 10.0]u"nm"
    @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ [5.0, 3.0]u"nm"
end

@testitem "translation image" begin
    using StaticArrays
    x = SVector{2}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    unitcell = [10.0 0.0; 0.0 10.0]
    @test CellListMap.translation_image(x, unitcell, [0, 0]) ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    @test CellListMap.translation_image(x, unitcell, [1, 1]) ≈ [[11.0, 11.0], [12.0, 12.0], [13.0, 13.0]]
    @test CellListMap.translation_image(x, unitcell, [1, 2]) ≈ [[11.0, 21.0], [12.0, 22.0], [13.0, 23.0]]
    @test CellListMap.translation_image(x, unitcell, [2, 1]) ≈ [[21.0, 11.0], [22.0, 12.0], [23.0, 13.0]]
    @test CellListMap.translation_image(x, unitcell, [-1, -1]) ≈ [[-9.0, -9.0], [-8.0, -8.0], [-7.0, -7.0]]
end

@testitem "replicate system" begin
    using StaticArrays
    unitcell = [10.0 0.0; 0.0 10.0]
    x = SVector{2, Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    CellListMap.replicate_system!(x, unitcell, (1, 0))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [11.0, 1.0], [12.0, 2.0], [13.0, 3.0]]
    x = SVector{2, Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    CellListMap.replicate_system!(x, unitcell, (1, 1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [11.0, 11.0], [12.0, 12.0], [13.0, 13.0]]
    x = SVector{2, Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    CellListMap.replicate_system!(x, unitcell, (0, 1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [1.0, 11.0], [2.0, 12.0], [3.0, 13.0]]
    x = SVector{2, Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    CellListMap.replicate_system!(x, unitcell, (-1, -1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [-9.0, -9.0], [-8.0, -8.0], [-7.0, -7.0]]
    # with a matrix input
    x = [1.0 2.0 3.0; 1.0 2.0 3.0]
    y = CellListMap.replicate_system(x, unitcell, (1, 0))
    @test y ≈ [1.0 2.0 3.0 11.0 12.0 13.0; 1.0 2.0 3.0 1.0 2.0 3.0]
    # throw error if dimensions do not match
    x = SVector{2, Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    @test_throws DimensionMismatch CellListMap.replicate_system!(x, unitcell, (1, 0, 1))
end

@testitem "cell cartesian/linear indices" begin
    using StaticArrays
    @test CellListMap.cell_cartesian_indices(SVector(3, 3), 1) == CartesianIndex(1, 1)
    @test CellListMap.cell_cartesian_indices(SVector(3, 3), 2) == CartesianIndex(2, 1)
    @test CellListMap.cell_cartesian_indices(SVector(3, 3), 4) == CartesianIndex(1, 2)
    @test CellListMap.cell_cartesian_indices(SVector(3, 3), 9) == CartesianIndex(3, 3)
    @test CellListMap.cell_linear_index(SVector(3, 3), CartesianIndex(1, 1)) == 1
    @test CellListMap.cell_linear_index(SVector(3, 3), CartesianIndex(2, 1)) == 2
    @test CellListMap.cell_linear_index(SVector(3, 3), CartesianIndex(1, 2)) == 4
    @test CellListMap.cell_linear_index(SVector(3, 3), CartesianIndex(3, 3)) == 9
end

@testitem "_minmax" setup = [AllocTest] begin
    using BenchmarkTools
    using StaticArrays
    import CellListMap: _minmax
    x = [[0.0, 0.5, 1.0], [0.5, 1.0, 0.0], [1.0, 0.0, 0.5]]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    x = [SVector(0.0, 0.5, 1.0), SVector(0.5, 1.0, 0.0), SVector(1.0, 0.0, 0.5)]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    a = @ballocated _minmax($x) evals = 1 samples = 1
    @test a == Allocs(0)
end

@testitem "limits - invalid coordinates" begin
    using CellListMap
    using StaticArrays: SVector
    x = rand(SVector{3, Float64}, 100)
    x[50] = SVector(1.0, NaN, 1.0)
    y = rand(SVector{3, Float64}, 100)
    @test_throws ArgumentError CellListMap.limits(x)
    @test_throws ArgumentError CellListMap.limits(x, y)
    @test_throws ArgumentError CellListMap.limits(y, x)
    x = rand(3, 100)
    x[2, 50] = NaN
    y = rand(3, 100)
    @test_throws ArgumentError CellListMap.limits(x)
    @test_throws ArgumentError CellListMap.limits(x, y)
    @test_throws ArgumentError CellListMap.limits(y, x)
end

@testitem "align_cell" setup = [TestingNeighborLists] begin
    import CellListMap: align_cell
    using StaticArrays
    using LinearAlgebra

    l = sqrt(2) / 2

    m = @SMatrix[1.0 0.0; 0.0 1.0]
    @test align_cell(m) == (one(m), one(m))

    m = @SMatrix[l 0; l 1]
    mt, R = align_cell(m)
    @test mt ≈ [1 l; 0 l]
    @test R ≈ [l l; -l l]

    m = @SMatrix[-l 0; l 1]
    mt, R = align_cell(m)
    @test mt ≈ [1 l; 0 -l]
    @test R ≈ [-l l; -l -l]

    #! format: off
    m = @SMatrix[ 
        1  0  0
        0  1  0 
        0  0  1
    ]
    #! format: on
    @test align_cell(m) == (one(m), one(m))

    #! format: off
    m = @SMatrix[ 
        3  0  0
        0  2  0 
        0  0  1
    ]
    #! format: on
    for _ in 1:5
        local R
        R = random_rotation()
        mr = R * m
        ma, Ra = align_cell(mr)
        @test ma[:, 1] ≈ m[:, 1]
        @test cross([1, 0, 0], cross(ma[:, 2], ma[:, 3])) ≈ zeros(3) atol = 1.0e-10
    end

    # throw error if not 2D or 3D
    @test_throws ArgumentError align_cell(rand(4, 4))

end

@testitem "cell_limits" begin
    import CellListMap: cell_limits, align_cell
    m = [1 0; 0 1]
    @test cell_limits(m) == ([0.0, 0.0], [1.0, 1.0])

    m = [10 5; 5 10]
    @test cell_limits(m) == ([0.0, 0.0], [15.0, 15.0])

    mr, _ = align_cell(m)
    cl = cell_limits(mr)
    @test cl[1] ≈ [0.0, 0.0] atol = 1.0e-10
    @test cl[2] ≈ [20.12461179749811, 6.708203932499369]

    m = [10 5; 0 10]
    @test cell_limits(m) == ([0.0, 0.0], [15.0, 10.0])

    mr, _ = align_cell(m)
    cl = cell_limits(mr)
    @test cl[1] ≈ [0.0, -8.94427190999916]
    @test cl[2] ≈ [15.652475842498529, 0.0]

    m = [1 0 0; 0 1 0; 0 0 1]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])

    m = [1 0 0; 0 2 0; 0 0 1]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 2.0, 1.0])

    mr, _ = align_cell(m)
    @test cell_limits(mr) == ([0.0, -1.0, 0.0], [2.0, 0.0, 1.0])

    m = [1 0 0; 0 2 0; 0 0 3]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])

    mr, _ = align_cell(m)
    @test cell_limits(mr) == ([0.0, 0.0, -1.0], [3.0, 2.0, 0.0])
end

@testitem "cell vertices" begin
    import CellListMap: cell_vertices
    m = [1 0; 0 1]
    @test all(
        isapprox.(cell_vertices(m), [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], atol = 1.0e-6)
    )
    m = [1 0.5; 0 1]
    @test all(
        isapprox.(cell_vertices(m), [[0.0, 0.0], [1.0, 0.0], [1.5, 1.0], [0.5, 1.0]], atol = 1.0e-6)
    )
    m = [1 0 0; 0 1 0; 0 0 1]
    @test all(
        isapprox.(
            cell_vertices(m),
            [
                [0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [1, 0, 1],
                [0, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
            atol = 1.0e-6
        )
    )
    m = [1 0.5 0; 0 1 0; 0 0 1]
    @test all(
        isapprox.(
            cell_vertices(m),
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.5, 1.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.5, 1.0, 1.0],
                [1.5, 1.0, 1.0],
            ],
            atol = 1.0e-6
        )
    )
    # throw error if not 2D or 3D
    @test_throws ArgumentError cell_vertices(rand(4, 4))
end
