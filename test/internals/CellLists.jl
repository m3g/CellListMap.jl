@testitem "output of nbatches" begin
    using CellListMap, StaticArrays
    const PSP = CellListMap.ParticleSystemPositions
    x = rand(3, 100)
    y = rand(3, 10)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(PSP(x), box; nbatches = (2, 4))
    @test CellListMap.nbatches(cl) == (2, 4)
    @test CellListMap.nbatches(cl, :build) == 2
    @test CellListMap.nbatches(cl, :map) == 4
    cl = CellListMap.CellList(PSP(x), PSP(y), box; nbatches = (2, 4))
    @test CellListMap.nbatches(cl) == (2, 4)
    @test CellListMap.nbatches(cl, :build) == 2
    @test CellListMap.nbatches(cl, :map) == 4
    cl = CellListMap.CellList(PSP(x), PSP(y), box)
    # For CellListPair, build batches are independent, but map batches should match
    @test CellListMap.nbatches(cl.ref_list, :map) == CellListMap.nbatches(cl.target_list, :map)
    cl = CellListMap.CellList(PSP(x), box; nbatches = (2, 4), parallel = false)
    @test CellListMap.nbatches(cl) == (1, 1)
    # The automatic set of number of batches for this small system:
    if Threads.nthreads() == 10
        x = rand(SVector{3, Float64}, 10)
        sys = ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
        @test CellListMap.nbatches(sys) == (4, 4)
        x = rand(SVector{3, Float64}, 10000)
        resize!(sys.xpositions, length(x))
        sys.xpositions .= x
        pairwise!((pair, out) -> out += pair.d2, sys)
        @test CellListMap.nbatches(sys) == (10, 10)
    else
        @warn "Test not run because it is invalid for this number of threads, requires Threads.nthreads() = 10"
    end
end

@testitem "automatic nbatches update on UpdateCellList!" begin
    using CellListMap, StaticArrays
    # 3-arg UpdateCellList! (without preallocated aux) should update nbatches
    # when the number of particles changes
    x = rand(SVector{3, Float64}, 2)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(x, box)
    nb_initial = CellListMap.nbatches(cl)
    x = rand(SVector{3, Float64}, 10000)
    cl = CellListMap.UpdateCellList!(x, box, cl)
    expected = (
        CellListMap._nbatches_build_cell_lists(10000),
        CellListMap._nbatches_map_computation(10000),
    )
    @test CellListMap.nbatches(cl) == expected
    # nbatches should not change when particle count is unchanged
    x .= rand(SVector{3, Float64}, 10000)
    cl = CellListMap.UpdateCellList!(x, box, cl)
    @test CellListMap.nbatches(cl) == expected
end

@testitem "nbatches symbol variants, show, and CellListPair" begin
    using CellListMap, StaticArrays
    x = rand(SVector{3, Float64}, 100)
    y = rand(SVector{3, Float64}, 50)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(x, box; nbatches = (2, 4))
    # Full symbol names
    @test CellListMap.nbatches(cl, :map_computation) == 4
    @test CellListMap.nbatches(cl, :build_cell_lists) == 2
    # show method for NumberOfBatches
    @test sprint(show, MIME"text/plain"(), cl.nbatches) !== ""
    # Auto nbatches with parallel=true (unconditional)
    cl = CellListMap.CellList(x, box)
    nb = CellListMap.nbatches(cl)
    @test nb[1] >= 1 && nb[1] <= min(8, Threads.nthreads())
    @test nb[2] >= 1
    # nbatches on CellListPair
    cl_pair = CellListMap.CellList(x, y, box)
    @test CellListMap.nbatches(cl_pair) == CellListMap.nbatches(cl_pair.ref_list)
    @test CellListMap.nbatches(cl_pair, :build) == CellListMap.nbatches(cl_pair.ref_list, :build)
    @test CellListMap.nbatches(cl_pair, :map) == CellListMap.nbatches(cl_pair.ref_list, :map)
end

@testitem "set_idxs! error on nbatches mismatch" begin
    using CellListMap
    idxs = [1:10, 11:20]
    @test_throws ArgumentError CellListMap.set_idxs!(idxs, 100, 3)
end

@testitem "UpdateCellList! for CellListPair" begin
    using CellListMap, StaticArrays
    x = rand(SVector{3, Float64}, 100)
    y = rand(SVector{3, Float64}, 50)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(x, y, box)
    # 3-arg update, serial path
    x2 = rand(SVector{3, Float64}, 100)
    y2 = rand(SVector{3, Float64}, 50)
    cl = CellListMap.UpdateCellList!(x2, y2, box, cl; parallel = false)
    @test cl.ref_list.n_real_particles == 100  # x (first arg) has 100
    @test cl.target_list.n_real_particles == 50  # y (second arg) has 50
    # 3-arg update, parallel path
    x3 = rand(SVector{3, Float64}, 100)
    y3 = rand(SVector{3, Float64}, 50)
    cl = CellListMap.UpdateCellList!(x3, y3, box, cl; parallel = true)
    @test cl.ref_list.n_real_particles == 100  # x (first arg) has 100
    @test cl.target_list.n_real_particles == 50  # y (second arg) has 50
end

@testitem "celllists - validate coordinates" begin
    using StaticArrays
    x = rand(SVector{3, Float64}, 100)
    x[50] = SVector(1.0, NaN, 1.0)
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1)
    y = rand(SVector{3, Float64}, 100)
    @test_throws ArgumentError CellListMap.CellList(x, box)
    cl = CellListMap.CellList(y, box)
    @test_throws ArgumentError CellListMap.UpdateCellList!(x, box, cl)
    @test_throws ArgumentError CellListMap.CellList(x, y, box)
    @test_throws ArgumentError CellListMap.CellList(y, x, box)
    cl = CellListMap.CellList(y, y, box)
    @test_throws ArgumentError CellListMap.UpdateCellList!(x, y, box, cl)
    @test_throws ArgumentError CellListMap.UpdateCellList!(y, x, box, cl)
    x = rand(3, 100)
    x[2, 50] = NaN
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1)
    y = rand(3, 100)
    @test_throws ArgumentError CellListMap.CellList(x, box)
    cl = CellListMap.CellList(y, box)
    @test_throws ArgumentError CellListMap.UpdateCellList!(x, box, cl)
    @test_throws ArgumentError CellListMap.CellList(x, y, box)
    @test_throws ArgumentError CellListMap.CellList(y, x, box)
    cl = CellListMap.CellList(y, y, box)
    @test_throws ArgumentError CellListMap.UpdateCellList!(x, y, box, cl)
    @test_throws ArgumentError CellListMap.UpdateCellList!(y, x, box, cl)
end

@testitem "CellList dimension mismatch error" begin
    using CellListMap, StaticArrays
    # 3D positions with 2D box should throw DimensionMismatch
    x = rand(SVector{3, Float64}, 100)
    box = CellListMap.Box([1.0, 1.0], 0.1)
    @test_throws DimensionMismatch CellListMap.CellList(x, box)
    # 2D positions with 3D box should throw DimensionMismatch
    x = rand(SVector{2, Float64}, 100)
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1)
    @test_throws DimensionMismatch CellListMap.CellList(x, box)
    # Matrix input: 3 rows (3D) with 2D box
    x = rand(3, 100)
    box = CellListMap.Box([1.0, 1.0], 0.1)
    @test_throws DimensionMismatch CellListMap.CellList(x, box)
end

@testitem "real_particle_border_case" begin
    using CellListMap, StaticArrays
    # Test that particles at exact boundaries are handled correctly
    # The function adjusts cell indices at the boundary cells (lcell and nc - lcell + 1)

    # Create a box where we can control the boundary conditions
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1, lcell = 1)

    # Test lower boundary case: particle at cell index == lcell should be moved to lcell + 1
    cartesian_index = CartesianIndex(box.lcell, box.lcell + 1, box.lcell + 2)
    result = CellListMap.real_particle_border_case(cartesian_index, box)
    @test result[1] == box.lcell + 1  # was at lcell, moved to lcell + 1
    @test result[2] == box.lcell + 1  # unchanged
    @test result[3] == box.lcell + 2  # unchanged

    # Test upper boundary case: particle at cell index == nc - lcell + 1 should be moved to nc - lcell
    upper_boundary = box.nc[1] - box.lcell + 1
    cartesian_index = CartesianIndex(upper_boundary, box.lcell + 1, box.lcell + 2)
    result = CellListMap.real_particle_border_case(cartesian_index, box)
    @test result[1] == upper_boundary - 1  # was at upper boundary, moved down
    @test result[2] == box.lcell + 1       # unchanged
    @test result[3] == box.lcell + 2       # unchanged

    # Test both boundaries in different dimensions
    cartesian_index = CartesianIndex(box.lcell, box.nc[2] - box.lcell + 1, box.lcell + 5)
    result = CellListMap.real_particle_border_case(cartesian_index, box)
    @test result[1] == box.lcell + 1                # lower boundary adjusted
    @test result[2] == box.nc[2] - box.lcell        # upper boundary adjusted
    @test result[3] == box.lcell + 5                # unchanged (not at boundary)
end
