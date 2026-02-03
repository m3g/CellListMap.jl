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
        # The gromacs unitcell is row-major as it comes from Chemfiles
        unitcell = Matrix(transpose(unitcell))
        Chemfiles.close(traj)
        return unitcell, reinterpret(reshape, SVector{3, Float64}, Chemfiles.positions(frame))
    end

    function test_newcl(file, lcell; cutoff = 8.0)
        unitcell, coordinates = getcoor(file)
        box = CellListMap.Box(unitcell, cutoff, lcell = lcell)
        cl = CellListMap.CellList(coordinates, box)
        u = pairwise!((pair, u) -> lj_Argon_Gromacs(pair.d2, u), 0.0, box, cl)
        return u
    end

    function grep_energy(file)
        found = false
        for line in eachline(file)
            if occursin("LJ (SR)", line)
                found = true
                continue
            end
            if found
                data = split(line)
                return parse(Float64, data[1])
            end
        end
        error("Could not find LJ (SR) in $file")
    end

    dir = "$(@__DIR__)/argon"

    minimal = grep_energy("$dir/minimal.log")
    cubic = grep_energy("$dir/cubic.log")
    dodecahedral = grep_energy("$dir/dodecahedron.log")
    octahedral = grep_energy("$dir/octahedron.log")
    for lcell in [1, 2, 3]
        @test test_newcl("$dir/minimal.xtc", lcell) ≈ minimal rtol = 1.0e-5
        @test test_newcl("$dir/cubic.xtc", lcell) ≈ cubic rtol = 1.0e-5
        @test test_newcl("$dir/dodecahedron.xtc", lcell) ≈ dodecahedral rtol = 1.0e-5
        @test test_newcl("$dir/octahedron.xtc", lcell) ≈ octahedral rtol = 1.0e-5
    end

end
