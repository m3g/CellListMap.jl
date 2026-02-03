@testitem "CellListMap.unitcelltype" begin
    using CellListMap
    @test CellListMap.unitcelltype(CellListMap.Box([1, 1, 1], 0.1)) == CellListMap.OrthorhombicCell
    @test CellListMap.unitcelltype(CellListMap.Box([1 0 0; 0 1 0; 0 0 1], 0.1)) == CellListMap.TriclinicCell
    x = rand(3, 100)
    @test CellListMap.unitcelltype(CellListMap.Box(CellListMap.limits(x), 0.1)) == CellListMap.NonPeriodicCell
end

@testitem "promote types" begin
    import CellListMap: _promote_types
    T = _promote_types([1, 1, 1], 0.1)
    @test T == Float64
    T = _promote_types([1, 1, 1], 0.1f0)
    @test T == Float32
    T = _promote_types([1.0, 1, 1], 0.1f0)
    @test T == Float64
    T = _promote_types([1.0f0, 1, 1], 0.1f0)
    @test T == Float32
    T = _promote_types([10, 10, 10.0], 1)
    @test T == Float64
    T = _promote_types([10.0f0, 10, 10], 1)
    @test T == Float32
    T = _promote_types([10, 10, 10], 1)
    @test T == Float64
end

@testitem "cell_matrix_from_sides: units" begin
    using StaticArrays
    using Unitful
    import CellListMap: cell_matrix_from_sides

    sides = SVector(1.0, 1.0, 1.0)u"nm"
    cell_matrix = cell_matrix_from_sides(sides)
    @test eltype(sides) == eltype(cell_matrix)

    sides = (1.0u"nm", 1.0u"nm", 1.0u"nm")
    cell_matrix = cell_matrix_from_sides(sides)
    @test eltype(sides) == eltype(cell_matrix)
end

@testitem "Update box with Limits" begin
    using CellListMap
    r = [[1.0, 1.0, 1.0]]
    system = InPlaceNeighborList(x = r, cutoff = 3.0, parallel = false)
    list = neighborlist!(system)
    update!(system, r)
    @test list == Tuple{Int64, Int64, Float64}[]
    r = [[1.0, 1.0, 1.0], [10.0, 1.0, 1.0], [3.0, 1.0, 1.0]]
    update!(system, r)
    list = neighborlist!(system)
    @test list == [(1, 3, 2.0)]
    r = [[7, 10, 10], [18, 10, 10]]
    system = InPlaceNeighborList(x = r, cutoff = 3.0, parallel = false)
    list = neighborlist!(system)
    @test list == Tuple{Int64, Int64, Float64}[]
    r = [[10.89658911843461, 3.709237933444153, 10.0], [13.894156281793144, 11.054172259416013, 10.0]]
    update!(system, r)
    list = neighborlist!(system)
    @test list == Tuple{Int64, Int64, Float64}[]
end

@testitem "Stable Box update" setup = [AllocTest] begin
    using CellListMap
    using StaticArrays
    using BenchmarkTools
    using LinearAlgebra: diag

    # test box construction with different types
    box = CellListMap.Box([1, 1, 1], 0.1f0)
    @test typeof(box.cutoff) == Float32
    @test eltype(box.input_unit_cell.matrix) == Float32
    box = CellListMap.Box([1, 1, 1], 0.1)
    @test typeof(box.cutoff) == Float64
    @test eltype(box.input_unit_cell.matrix) == Float64
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1f0)
    @test typeof(box.cutoff) == Float64
    @test eltype(box.input_unit_cell.matrix) == Float64
    box = CellListMap.Box([10, 10, 10], 1)
    @test typeof(box.cutoff) == Float64
    @test eltype(box.input_unit_cell.matrix) == Float64
    box = CellListMap.Box([10.0f0, 10.0f0, 10.0f0], 1)
    @test typeof(box.cutoff) == Float32
    @test eltype(box.input_unit_cell.matrix) == Float32

    # update with tuples
    box = CellListMap.Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell = (2, 2, 2), cutoff = 0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell = (2, 2, 2), cutoff = 0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with SVector
    box = CellListMap.Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell = SVector(2, 2, 2), cutoff = 0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell = (2, 2, 2), cutoff = 0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with Limits
    x = rand(SVector{3, Float64}, 1000)
    box = CellListMap.Box(CellListMap.limits(x), 0.1)
    new_x = rand(SVector{3, Float64}, 1500)
    a = @ballocated CellListMap.update_box($box; unitcell = $(CellListMap.limits(new_x)), cutoff = 0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell = CellListMap.limits(new_x), cutoff = 0.2)
    @test new_box.cutoff == 0.2
    @test diag(new_box.input_unit_cell.matrix) ≈ CellListMap.limits(new_x).limits .+ 2.1 * 0.2

    # Update with SMatrix
    box = CellListMap.Box([1 0 0; 0 1 0; 0 0 1], 0.1)
    new_matrix = SMatrix{3, 3, Float64, 9}(2, 0, 0, 0, 2, 0, 0, 0, 2)
    a = @ballocated CellListMap.update_box($box; unitcell = $new_matrix, cutoff = 0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell = new_matrix, cutoff = 0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

end
