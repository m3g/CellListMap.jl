@testitem "compare_with_gromacs.jl" begin

    using Test
    import Chemfiles
    using CellListMap
    using StaticArrays

    function lj_Argon_Gromacs(d2, u)
        c6 = 0.00622127e6
        c12 = 9.69576e6
        u += c12 / d2^6 - c6 / d2^3
    end

    function getcoor(file)
        traj = redirect_stdout(() -> Chemfiles.Trajectory(file), devnull)
        frame = Chemfiles.read_step(traj, 0)
        unitcell = Chemfiles.matrix(Chemfiles.UnitCell(frame))
        Chemfiles.close(traj)
        return unitcell, reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
    end

    function test_newcl(file, lcell)
        unitcell, coordinates = getcoor(file)
        unitcell = transpose(unitcell)
        box = Box(unitcell, 8.0, lcell=lcell)
        cl = CellList(coordinates, box)
        u = map_pairwise!((x, y, i, j, d2, u) -> lj_Argon_Gromacs(d2, u), 0.0, box, cl)
        return u
    end

    function check_unit_cell(m::AbstractMatrix)
        x, y, z = 1, 2, 3
        a = @view(m[:,1])
        b = @view(m[:,2])
        c = @view(m[:,3])
        if a[y] == 0 &&
           a[z] == 0 &&
           b[z] == 0 &&
           a[x] > 0  &&
           b[y] > 0  &&
           c[z] > 0  &&
           abs(b[x]) <= 0.5 * a[x] &&
           abs(c[x]) <= 0.5 * a[x] &&
           abs(c[y]) <= 0.5 * b[y] 
           return true
        else
           return false
        end
    end

    dir = @__DIR__

    lcell = 1

    # Minimal (2 atom) orthorhombic box
    correct = -0.358447
    @test test_newcl("$dir/argon_minimal/traj_comp.xtc", lcell) ≈ correct atol=1e-6

    # Some orthorhombic cells
    correct = -0.439334
    @test test_newcl("$dir/argon_original/traj_comp.xtc", lcell) ≈ correct atol=1e-6

    # dodecahedral box
#    correct = 40774.9
#    @test test_newcl("$dir/argon/traj_comp.xtc", lcell) ≈ correct atol=1e-6


end