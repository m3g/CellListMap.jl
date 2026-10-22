@testitem "output of nbatches" begin
    using CellListMap, StaticArrays
    const PSP = CellListMap.ParticleSystemPositions
    x = rand(SVector{3,Float64}, 100)
    y = rand(SVector{3,Float64}, 10)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(PSP(x), box; nbatches=(2, 4))
    @test CellListMap.nbatches(cl) == (2, 4)
    @test CellListMap.nbatches(cl, :build) == 2
    @test CellListMap.nbatches(cl, :map) == 4
    cl = CellListMap.CellList(PSP(x), PSP(y), box; nbatches=(2, 4))
    @test CellListMap.nbatches(cl) == (2, 4)
    @test CellListMap.nbatches(cl, :build) == 2
    @test CellListMap.nbatches(cl, :map) == 4
    cl = CellListMap.CellList(PSP(x), PSP(y), box)
    # For CellListPair, build batches are independent, but map batches should match
    @test CellListMap.nbatches(cl.ref_list, :map) == CellListMap.nbatches(cl.target_list, :map)
    cl = CellListMap.CellList(PSP(x), box; nbatches=(2, 4), parallel=false)
    @test CellListMap.nbatches(cl) == (1, 1)
    # The automatic set of number of batches for this small system:
    if Threads.nthreads() == 10
        x = rand(SVector{3,Float64}, 10)
        sys = ParticleSystem(positions=x, cutoff=0.1, unitcell=[1, 1, 1], output=0.0)
        @test CellListMap.get_nbatches(sys) == (4, 4)
        x = rand(SVector{3,Float64}, 10000)
        resize!(sys.xpositions, length(x))
        sys.xpositions .= x
        pairwise!((pair, out) -> out += pair.d2, sys)
        @test CellListMap.get_nbatches(sys) == (10, 10)
    else
        @warn "Test not run because it is invalid for this number of threads, requires Threads.nthreads() = 10"
    end
end

@testitem "automatic nbatches update on UpdateCellList!" setup = [Testing] begin
    using CellListMap, StaticArrays
    # 3-arg UpdateCellList! (without preallocated aux) should update nbatches
    # when the number of particles changes
    x = rand(SVector{3,Float64}, 2)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(PSP(x), box)
    nb_initial = CellListMap.nbatches(cl)
    x = rand(SVector{3,Float64}, 10000)
    cl = CellListMap.UpdateCellList!(PSP(x), box, cl)
    expected = (
        CellListMap._nbatches_build_cell_lists(10000),
        CellListMap._nbatches_map_computation(10000),
    )
    @test CellListMap.nbatches(cl) == expected
    # nbatches should not change when particle count is unchanged
    x .= rand(SVector{3,Float64}, 10000)
    cl = CellListMap.UpdateCellList!(PSP(x), box, cl)
    @test CellListMap.nbatches(cl) == expected

    # Test updating to very few particles, in parallel runs (caused crash in a bug)
    sys = ParticleSystem(
        xpositions=rand(3, 10),
        ypositions=rand(3, 100),
        cutoff=0.2,
        output=0.0,
        unitcell=[1, 1, 1],
        parallel=true,
        nbatches=(4, 4),
    )
    @test CellListMap.get_nbatches(sys) == (4,4)
    pairwise!((pair, sumd) -> sumd += pair.d, sys)
    update!(sys; xpositions=rand(3,1))
    @test pairwise!((pair, sumd) -> sumd += pair.d, sys) > -1 
    update!(sys; xpositions=rand(3,100), ypositions=rand(3,1))
    @test pairwise!((pair, sumd) -> sumd += pair.d, sys) > -1 
    update!(sys; xpositions=rand(3,1), ypositions=rand(3,1))
    @test pairwise!((pair, sumd) -> sumd += pair.d, sys) > -1
    @test CellListMap.get_nbatches(sys) == (4,4)

end

@testitem "nbatches symbol variants, show, and CellListPair" setup = [Testing] begin
    using CellListMap, StaticArrays
    x = rand(SVector{3,Float64}, 100)
    y = rand(SVector{3,Float64}, 50)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(PSP(x), box; nbatches=(2, 4))
    # Full symbol names
    @test CellListMap.nbatches(cl, :map_computation) == 4
    @test CellListMap.nbatches(cl, :build_cell_lists) == 2
    # show method for NumberOfBatches
    @test sprint(show, MIME"text/plain"(), cl.nbatches) !== ""
    # Auto nbatches with parallel=true (unconditional)
    cl = CellListMap.CellList(PSP(x), box)
    nb = CellListMap.nbatches(cl)
    @test nb[1] >= 1 && nb[1] <= min(8, Threads.nthreads())
    @test nb[2] >= 1
    # nbatches on CellListPair
    cl_pair = CellListMap.CellList(PSP(x), PSP(y), box)
    @test CellListMap.nbatches(cl_pair) == CellListMap.nbatches(cl_pair.ref_list)
    @test CellListMap.nbatches(cl_pair, :build) == CellListMap.nbatches(cl_pair.ref_list, :build)
    @test CellListMap.nbatches(cl_pair, :map) == CellListMap.nbatches(cl_pair.ref_list, :map)
end

@testitem "UpdateCellList! for CellListPair" setup = [Testing] begin
    using CellListMap, StaticArrays
    x = rand(SVector{3,Float64}, 100)
    y = rand(SVector{3,Float64}, 50)
    box = CellListMap.Box([1, 1, 1], 0.1)
    cl = CellListMap.CellList(PSP(x), PSP(y), box)
    # 3-arg update, serial path
    x2 = rand(SVector{3,Float64}, 100)
    y2 = rand(SVector{3,Float64}, 50)
    cl = CellListMap.UpdateCellList!(PSP(x2), PSP(y2), box, cl; parallel=false)
    @test cl.ref_list.n_real_particles == 100  # x (first arg) has 100
    @test cl.target_list.n_real_particles == 50  # y (second arg) has 50
    # 3-arg update, parallel path
    x3 = rand(SVector{3,Float64}, 100)
    y3 = rand(SVector{3,Float64}, 50)
    cl = CellListMap.UpdateCellList!(PSP(x3), PSP(y3), box, cl; parallel=true)
    @test cl.ref_list.n_real_particles == 100  # x (first arg) has 100
    @test cl.target_list.n_real_particles == 50  # y (second arg) has 50
end

@testitem "celllists - validate coordinates" setup = [Testing] begin
    using StaticArrays
    x = rand(SVector{3,Float64}, 100)
    x[50] = SVector(1.0, NaN, 1.0)
    @test_throws "Invalid coordinates" ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    x[50] = SVector(1.0, 0.5, 1.0)
    sys = ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    sys.xpositions[1] = SVector(1.0, NaN, 1.0)
    @test_throws "Invalid coordinates" pairwise!(f1, sys)
    x[50] = SVector(1.0, 0.5, 1.0)
    y = rand(SVector{3,Float64}, 100)
    y[50] = SVector(1.0, NaN, 1.0)
    @test_throws "Invalid coordinates" ParticleSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    y[50] = SVector(1.0, 0.5, 1.0)
    sys = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    sys.ypositions[1] = SVector(1.0, NaN, 1.0)
    @test_throws "Invalid coordinates" pairwise!(f1, sys)
    x = rand(3, 100)
    x[2, 50] = NaN
    @test_throws "Invalid coordinates" ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
end

@testitem "celllists - validate coordinates in parallel builds" setup = [Testing] begin
    using StaticArrays
    # The default validation is done while locating the particles, thus the invalid
    # coordinate must be detected by any of the tasks of the parallel builds: the
    # parallel location of the serial counting sort (10_000 particles) and the
    # parallel counting sort (100_000 particles).
    for n in (10_000, 100_000), i in (1, n ÷ 2, n)
        x = rand(SVector{3,Float64}, n)
        x[i] = SVector(0.5, NaN, 0.5)
        @test_throws "Invalid coordinates" ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.05, output=0.0)
        x[i] = SVector(0.5, 0.5, 0.5)
        sys = ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.05, output=0.0)
        x[i] = SVector(0.5, NaN, 0.5)
        @test_throws "Invalid coordinates" update!(sys; xpositions=x, rebuild=true)
        # the list is rebuilt on the next call after the coordinates are fixed
        x[i] = SVector(0.5, 0.5, 0.5)
        update!(sys; xpositions=x)
        @test pairwise!((pair, out) -> out + 1, sys) == pairwise!((pair, out) -> out + 1,
            ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.05, output=0, parallel=false))
    end
    # Reference list of a two-set system with few particles relative to the number of cells
    y = rand(SVector{3,Float64}, 10_000)
    x = rand(SVector{3,Float64}, 10)
    sys = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.05, output=0.0)
    x[5] = SVector(0.5, NaN, 0.5)
    @test_throws "Invalid coordinates" update!(sys; xpositions=x, rebuild=true)
end

@testitem "CellList dimension mismatch error" begin
    using CellListMap, StaticArrays
    PSP = CellListMap.ParticleSystemPositions
    # 3D positions with 2D box should throw DimensionMismatch
    x = rand(SVector{3,Float64}, 100)
    box = CellListMap.Box([1.0, 1.0], 0.1)
    @test_throws DimensionMismatch CellListMap.CellList(PSP(x), box)
    # 2D positions with 3D box should throw DimensionMismatch
    x = rand(SVector{2,Float64}, 100)
    box = CellListMap.Box([1.0, 1.0, 1.0], 0.1)
    @test_throws DimensionMismatch CellListMap.CellList(PSP(x), box)
end

@testitem "parallel counting sort" begin
    using CellListMap, StaticArrays
    PSP = CellListMap.ParticleSystemPositions
    # Number of batches of the counting sort
    @test CellListMap._nbatches_counting_sort(8, 10_000, 100) == 1
    @test CellListMap._nbatches_counting_sort(8, 10^6, 10^5) == 8
    @test CellListMap._nbatches_counting_sort(8, 10^6, 2 * 10^6) == 4
    @test CellListMap._nbatches_counting_sort(8, 10^6, 3 * 10^6) == 1
    @test CellListMap._nbatches_counting_sort(1, 10^6, 10^5) == 1
    # The parallel build must be identical to the serial one (called directly, to
    # run independently of the number of threads)
    function same_build(x, box)
        cl = CellListMap.CellList(PSP(x), box; parallel = false)
        ref = deepcopy(cl)
        for nb in (2, 3, 8)
            CellListMap._build_parallel!(cl, PSP(x), box, nb)
            cl.particles == ref.particles || return false
            cl.shifts == ref.shifts || return false
            cl.cell_start == ref.cell_start || return false
        end
        return true
    end
    x = [10 * rand(SVector{3, Float64}) .- 2 for _ in 1:5000]
    @test same_build(x, CellListMap.Box([8, 8, 8], 1.0))
    @test same_build(x, CellListMap.Box([8, 8, 8], 1.0; lcell = 2))
    @test same_build(x, CellListMap.Box([8 1 0; 0 8 1; 0 0 8], 1.0))
    @test same_build(x, CellListMap.Box(CellListMap.limits(PSP(x)), 1.0))
    x = [10 * rand(SVector{2, Float64}) for _ in 1:5000]
    @test same_build(x, CellListMap.Box([10, 10], 0.5))
    @test same_build(x[1:3], CellListMap.Box([10, 10], 0.5))
    # Through UpdateCellList!, with enough particles to use the parallel sort
    x = [50 * rand(SVector{3, Float64}) for _ in 1:CellListMap.COUNTING_SORT_MIN_PARTICLES]
    box = CellListMap.Box([50, 50, 50], 2.0)
    cl_serial = CellListMap.CellList(PSP(x), box; parallel = false)
    cl_parallel = CellListMap.CellList(PSP(x), box; parallel = true)
    @test cl_parallel.particles == cl_serial.particles
    @test cl_parallel.shifts == cl_serial.shifts
    @test cl_parallel.cell_start == cl_serial.cell_start
    @test cl_parallel.nonempty_cells == cl_serial.nonempty_cells
end
