@testitem "compare_with_namd.jl" begin

    using Test
    import Chemfiles
    using CellListMap
    using StaticArrays

    function lj_NE(d2, u)
        ε = 0.0441795
        σ = 2 * 1.64009
        d = sqrt(d2)
        u += ε * ((σ / d)^12 - 2 * (σ / d)^6)
    end

    function getcoor(file)
        traj = redirect_stdout(() -> Chemfiles.Trajectory(file), devnull)
        frame = Chemfiles.read_step(traj, 0)
        Chemfiles.close(traj)
        return reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
    end

    function test_newcl(file, unit_cell, correct, lcell)
        coordinates = getcoor(file)
        box = Box(unit_cell, 10.0, lcell=lcell)
        cl = CellList(coordinates, box)
        u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
        return u ≈ correct
    end

    dir = @__DIR__

    lcell = 1

    # Some orthorhombic cells

    unit_cell = [50.0, 50.0, 50.0]
    correct = 32230.01699504111
    @test test_newcl("$dir/o1.dcd", unit_cell, correct, lcell)

    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    @test test_newcl("$dir/o2.dcd", unit_cell, correct, lcell)

    # Orthorhombic but rotated

    unit_cell = [50.0 0.0 50.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 1724.3195067566828
    @test test_newcl("$dir/o3.dcd", unit_cell, correct, lcell)

    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1754.0802503953591
    @test test_newcl("$dir/o4.dcd", unit_cell, correct, lcell)

    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1765.1389457850137
    @test test_newcl("$dir/o5.dcd", unit_cell, correct, lcell)

    unit_cell = [80.0, 80.0, 80.0]
    correct = -158.04751357760088
    @test test_newcl("$dir/o6.dcd", unit_cell, correct, lcell)

    # Some triclinic cells

    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    @test test_newcl("$dir/t1.dcd", unit_cell, correct, lcell)

    unit_cell = [50.0 0.0 0.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 32096.48839031735
    @test test_newcl("$dir/t2.dcd", unit_cell, correct, lcell)

    #
    # Check cell list updating routine
    #

    coordinates = getcoor("$dir/o1.dcd")
    unit_cell = [50.0 0.0 0.0
        0.0 50.0 0.0
        0.0 0.0 50.0]
    correct = 32230.01699504111
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = CellList(coordinates, box)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

    coordinates = getcoor("$dir/o2.dcd")
    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = UpdateCellList!(coordinates, box, cl)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

    # Test preallocated AuxThreaded struct
    aux = CellListMap.AuxThreaded(cl)
    coordinates = getcoor("$dir/t1.dcd")
    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = UpdateCellList!(coordinates, box, cl, aux)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

    # Same thing with a different lcell

    lcell = 3

    # Some orthorhombic cells

    unit_cell = [50.0, 50.0, 50.0]
    correct = 32230.01699504111
    @test test_newcl("$dir/o1.dcd", unit_cell, correct, lcell)

    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    @test test_newcl("$dir/o2.dcd", unit_cell, correct, lcell)

    # Orthorhombic but rotated

    unit_cell = [50.0 0.0 50.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 1724.3195067566828
    @test test_newcl("$dir/o3.dcd", unit_cell, correct, lcell)

    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1754.0802503953591
    @test test_newcl("$dir/o4.dcd", unit_cell, correct, lcell)

    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1765.1389457850137
    @test test_newcl("$dir/o5.dcd", unit_cell, correct, lcell)

    unit_cell = [80.0, 80.0, 80.0]
    correct = -158.04751357760088
    @test test_newcl("$dir/o6.dcd", unit_cell, correct, lcell)

    # Some triclinic cells

    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    @test test_newcl("$dir/t1.dcd", unit_cell, correct, lcell)

    unit_cell = [50.0 0.0 0.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 32096.48839031735
    @test test_newcl("$dir/t2.dcd", unit_cell, correct, lcell)

    #
    # Check cell list updating routine
    #

    coordinates = getcoor("$dir/o1.dcd")
    unit_cell = [50.0 0.0 0.0
        0.0 50.0 0.0
        0.0 0.0 50.0]
    correct = 32230.01699504111
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = CellList(coordinates, box)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

    coordinates = getcoor("$dir/o2.dcd")
    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = UpdateCellList!(coordinates, box, cl)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

    # Test preallocated AuxThreaded struct
    aux = CellListMap.AuxThreaded(cl)
    coordinates = getcoor("$dir/t1.dcd")
    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    box = Box(unit_cell, 10.0, lcell=lcell)
    cl = UpdateCellList!(coordinates, box, cl, aux)
    u = map_pairwise!((x, y, i, j, d2, u) -> lj_NE(d2, u), 0.0, box, cl)
    @test u ≈ correct

end