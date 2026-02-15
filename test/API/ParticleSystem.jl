@testitem "ParticleSystem properties" begin
    using CellListMap
    using StaticArrays

    if VERSION >= v"1.11"
        @test Base.ispublic(CellListMap, :AbstractParticleSystem)
        @test Base.ispublic(CellListMap, :ParticleSystem1)
        @test Base.ispublic(CellListMap, :ParticleSystem2)
    end

    sys = ParticleSystem(
        positions = rand(SVector{3, Float64}, 1000),
        cutoff = 0.1,
        unitcell = [1, 1, 1],
        output = 0.0,
        output_name = :test
    )
    @test length(sys.positions) == 1000
    @test sys.cutoff == 0.1
    @test sys.unitcell == @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test sys.output == 0
    @test sys.test == 0
    @test sys.parallel == true

    sys.parallel = false
    @test sys.parallel == false
    sys.cutoff = 0.2
    @test sys.cutoff == 0.2
    sys.positions[1] = SVector(0.0, 0.0, 0.0)
    @test sys.positions[1] == SVector(0.0, 0.0, 0.0)
    sys.unitcell = [1.2, 1.2, 1.2]
    @test sys.unitcell == @SMatrix [1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2]

    # test the construction with pathologically few particles
    for x in [
            SVector{3, Float64}[],
            Vector{Float64}[],
            Matrix{Float64}(undef, 3, 0),
            [rand(SVector{3, Float64})],
            [rand(3)],
            rand(3, 1),
        ]
        _sys = ParticleSystem(
            positions = x,
            cutoff = 0.1,
            unitcell = [1, 1, 1],
            output = 0.0,
            output_name = :test
        )
        @test CellListMap.pairwise!((pair, out) -> out += pair.d2, _sys) == 0.0
    end

    # unitcell type
    x = rand(SVector{3, Float64}, 100)
    @test CellListMap.unitcelltype(ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)) == CellListMap.OrthorhombicCell
    @test CellListMap.unitcelltype(ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1 0 0; 0 1 0; 0 0 1], output = 0.0)) == CellListMap.TriclinicCell
    @test CellListMap.unitcelltype(ParticleSystem(positions = x, cutoff = 0.1, output = 0.0)) == CellListMap.NonPeriodicCell

    # Argument errors
    @test_throws ArgumentError ParticleSystem(
        positions = rand(SVector{3, Float64}, 100),
        xpositions = rand(SVector{3, Float64}, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws ArgumentError ParticleSystem(
        ypositions = rand(SVector{3, Float64}, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        positions = rand(1, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(1, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(2, 100),
        ypositions = rand(1, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        positions = rand(2, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(2, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(2, 100),
        ypositions = rand(2, 100),
        cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        positions = rand(3, 100),
        cutoff = 0.1, unitcell = [1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(3, 100),
        cutoff = 0.1, unitcell = [1, 1], output = 0.0,
    )
    @test_throws "Incompatible dimensions" ParticleSystem(
        xpositions = rand(3, 100),
        ypositions = rand(3, 100),
        cutoff = 0.1, unitcell = [1, 1], output = 0.0,
    )

end

@testitem "reducer method basics and errors" begin
    using CellListMap
    @test CellListMap.reducer(1, 2) == 3
    @test CellListMap.copy_output(1) == 1
    @test CellListMap.reset_output(1) == 0
    struct A end
    @test_throws ArgumentError CellListMap.reset_output(A())
    @test_throws ArgumentError CellListMap.copy_output(A())
    @test_throws ArgumentError CellListMap.reducer(A(), A())
end

@testitem "update_unitcell!" setup = [AllocTest] begin
    using BenchmarkTools
    using LinearAlgebra: diag
    using StaticArrays
    using CellListMap
    x = rand(SVector{3, Float64}, 1000)
    sys1 = ParticleSystem(xpositions = x, unitcell = [1, 1, 1], cutoff = 0.1, output = 0.0)
    update_unitcell!(sys1, SVector(2, 2, 2))
    @test diag(sys1.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys1, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == Allocs(0)
    y = rand(SVector{3, Float64}, 1000)
    sys2 = ParticleSystem(xpositions = x, ypositions = y, unitcell = [1, 1, 1], cutoff = 0.1, output = 0.0)
    update_unitcell!(sys2, SVector(2, 2, 2))
    @test diag(sys2.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys2, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == Allocs(0)
    # Test throwing error on updating non-periodic unit cells
    sys = ParticleSystem(xpositions = x, cutoff = 0.1, output = 0.0)
    @test_throws ArgumentError update_unitcell!(sys, [1, 1, 1])
    sys = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, output = 0.0)
    @test_throws ArgumentError update_unitcell!(sys, [1, 1, 1])
end

@testitem "update_cutoff!" setup = [AllocTest, Testing] begin
    using BenchmarkTools
    using StaticArrays
    using CellListMap
    using PDBTools
    x = rand(SVector{3, Float64}, 1000)
    sys1 = ParticleSystem(xpositions = x, unitcell = [1, 1, 1], cutoff = 0.1, output = 0.0)
    update_cutoff!(sys1, 0.2)
    @test sys1.cutoff == 0.2
    a = @ballocated update_cutoff!($sys1, 0.1) evals = 1 samples = 1
    @test a == Allocs(0)
    y = rand(SVector{3, Float64}, 1000)
    sys2 = ParticleSystem(xpositions = x, ypositions = y, unitcell = [1, 1, 1], cutoff = 0.1, output = 0.0)
    update_cutoff!(sys2, 0.2)
    @test sys2.cutoff == 0.2
    a = @ballocated update_cutoff!($sys2, 0.1) evals = 1 samples = 1
    @test a == Allocs(0)

    # Update cutoff of non-periodic systems
    x = coor(read_pdb(CellListMap.argon_pdb_file))
    sys1 = ParticleSystem(xpositions = x, cutoff = 8.0, output = 0.0)
    @test CellListMap.unitcelltype(sys1) == CellListMap.NonPeriodicCell
    @test sys1.unitcell ≈ [35.63 0.0 0.0; 0.0 35.76 0.0; 0.0 0.0 35.79] atol = 1.0e-2
    update_cutoff!(sys1, 10.0)
    @test sys1.unitcell ≈ [39.83 0.0 0.0; 0.0 39.96 0.0; 0.0 0.0 39.99] atol = 1.0e-2
    a = @ballocated update_cutoff!($sys1, 8.0) evals = 1 samples = 1
    @test a == Allocs(0)
    sys2 = ParticleSystem(xpositions = x[1:50], ypositions = x[51:100], cutoff = 8.0, output = 0.0)
    @test CellListMap.unitcelltype(sys2) == CellListMap.NonPeriodicCell
    @test sys2.unitcell ≈ [35.63 0.0 0.0; 0.0 35.76 0.0; 0.0 0.0 35.79] atol = 1.0e-2
    update_cutoff!(sys2, 10.0)
    @test sys2.unitcell ≈ [39.83 0.0 0.0; 0.0 39.96 0.0; 0.0 0.0 39.99] atol = 1.0e-2
    a = @ballocated update_cutoff!($sys2, 8.0) evals = 1 samples = 1
    @test a == Allocs(0)
end

@testitem "get_computing_box" begin
    using StaticArrays
    using CellListMap
    x = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(xpositions = x, unitcell = [1, 1, 1], cutoff = 0.1, output = 0.0)
    @test CellListMap.get_computing_box(sys) == ([-0.1, -0.1, -0.1], [1.1, 1.1, 1.1])
end

@testitem "UpdateParticleSystem!" setup = [AllocTest] begin
    # this updates must be non-allocating in the serial case
    using BenchmarkTools
    using StaticArrays
    using CellListMap
    x = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(xpositions = x, unitcell = [1.0, 1.0, 1.0], cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)
    y = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(xpositions = x, ypositions = y, unitcell = [1.0, 1.0, 1.0], cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)

    # Test construction with more general abstract vectors
    x = @view(x[1:500])
    sys = ParticleSystem(xpositions = x, unitcell = [1.0, 1.0, 1.0], cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)
    y = @view(y[1:500])
    sys = ParticleSystem(xpositions = x, ypositions = y, unitcell = [1.0, 1.0, 1.0], cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)

    # Update with matrices
    x = rand(3, 500)
    sys = ParticleSystem(xpositions = x, unitcell = [1.0, 1.0, 1.0], cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)

    # Update non-periodic system
    x = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(xpositions = x, cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)
    y = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, output = 0.0, parallel = false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == Allocs(0)

end

@testitem "automatic nbatches update on ParticleSystem resize" begin
    using CellListMap, StaticArrays

    # ParticleSystem1: nbatches updates when particle count changes
    x = rand(SVector{3, Float64}, 2)
    sys = ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    x = rand(SVector{3, Float64}, 10000)
    resize!(sys.xpositions, length(x))
    sys.xpositions .= x
    pairwise!((pair, out) -> out += pair.d2, sys)
    expected = (
        CellListMap._nbatches_build_cell_lists(10000),
        CellListMap._nbatches_map_computation(10000),
    )
    @test CellListMap.nbatches(sys) == expected

    # ParticleSystem1: nbatches stable when particle count does not change
    x = rand(SVector{3, Float64}, 1000)
    sys = ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    nb_before = CellListMap.nbatches(sys)
    sys.xpositions .= rand(SVector{3, Float64}, 1000)
    pairwise!((pair, out) -> out += pair.d2, sys)
    @test CellListMap.nbatches(sys) == nb_before

    # ParticleSystem1: manually set nbatches are not overridden on resize
    x = rand(SVector{3, Float64}, 100)
    sys = ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0, nbatches = (2, 3))
    @test CellListMap.nbatches(sys) == (2, 3)
    x = rand(SVector{3, Float64}, 10000)
    resize!(sys.xpositions, length(x))
    sys.xpositions .= x
    pairwise!((pair, out) -> out += pair.d2, sys)
    @test CellListMap.nbatches(sys) == (2, 3)

    # ParticleSystem2: nbatches updates when particle count changes
    x = rand(SVector{3, Float64}, 2)
    y = rand(SVector{3, Float64}, 3)
    sys = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    x = rand(SVector{3, Float64}, 5000)
    y = rand(SVector{3, Float64}, 10000)
    resize!(sys.xpositions, length(x))
    sys.xpositions .= x
    resize!(sys.ypositions, length(y))
    sys.ypositions .= y
    pairwise!((pair, out) -> out += pair.d2, sys)
    # For CellListPair, nbatches(sys) returns ref_list (x) nbatches
    # Build: based on ref_list particle count, Map: based on product
    expected = (
        CellListMap._nbatches_build_cell_lists(5000),
        CellListMap._nbatches_map_computation(5000 * 10000),
    )
    @test CellListMap.nbatches(sys) == expected
end

@testitem "ParticleSystem - validate coordinates" begin
    using CellListMap
    using StaticArrays
    # 1-set system
    x = rand(SVector{3, Float64}, 100)
    x[50] = SVector(1.1, NaN, 1.1)
    @test_throws ArgumentError ParticleSystem(positions = x, cutoff = 0.1, output = 0.0)
    @test_throws ArgumentError ParticleSystem(positions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    y = rand(SVector{3, Float64}, 100)
    p = ParticleSystem(positions = copy(y), cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    p.positions .= x
    @test_throws ArgumentError pairwise!((pair, out) -> out += pair.d2, p)
    p = ParticleSystem(xpositions = copy(y), cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0, validate_coordinates = nothing)
    @test pairwise!((pair, out) -> out += pair.d2, p) > 0.0
    # 2-set system
    x = rand(SVector{3, Float64}, 100)
    x[50] = SVector(1.1, NaN, 1.1)
    y = rand(SVector{3, Float64}, 100)
    @test_throws ArgumentError ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, output = 0.0)
    @test_throws ArgumentError ParticleSystem(xpositions = y, ypositions = x, cutoff = 0.1, output = 0.0)
    @test_throws ArgumentError ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    @test_throws ArgumentError ParticleSystem(xpositions = y, ypositions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    p = ParticleSystem(xpositions = copy(y), ypositions = copy(y), cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    p.xpositions .= x
    @test_throws ArgumentError pairwise!((pair, out) -> out += pair.d2, p)
    p = ParticleSystem(xpositions = copy(y), ypositions = copy(y), cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    p.ypositions .= x
    @test_throws ArgumentError pairwise!((pair, out) -> out += pair.d2, p)
    p = ParticleSystem(xpositions = copy(y), ypositions = copy(y), cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0, validate_coordinates = nothing)
    @test pairwise!((pair, out) -> out += pair.d2, p) > 0.0
end

@testitem "reset output value" setup = [Testing] begin
    using CellListMap, StaticArrays
    x, box = xatomic(5000)
    uc = box.input_unit_cell.matrix
    sys = ParticleSystem(positions = x, unitcell = uc, cutoff = 12.0, output = 0.0)
    u = pairwise!((pair, u) -> u += inv(pair.d), sys)
    sys = ParticleSystem(positions = x, unitcell = uc, cutoff = 12.0, output = u)
    pairwise!((pair, u) -> u += inv(pair.d), sys; reset = false)
    @test sys.output ≈ 2 * u
    pairwise!((pair, u) -> u += inv(pair.d), sys)
    @test sys.output ≈ u
end
