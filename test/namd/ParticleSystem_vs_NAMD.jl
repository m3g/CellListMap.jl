#
# Compare with NAMD results, calculations performed through the 
# ParticleSystem interface
#
@testitem "ParticleSystem_vs_NAMD.jl" begin

    using Test
    import Chemfiles
    using CellListMap
    using StaticArrays
    
    ε = 0.0441795
    σ = 2 * 1.64009
    
    function lj_NE(d2, u)
        d = sqrt(d2)
        u += ε * ((σ / d)^12 - 2 * (σ / d)^6)
    end
    
    function getcoor(file)
        traj = redirect_stdout(() -> Chemfiles.Trajectory(file), devnull)
        frame = Chemfiles.read_step(traj, 0)
        Chemfiles.close(traj)
        return reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
    end
    
    function copy_to_svector(positions)
        if positions isa AbstractVector{<:AbstractVector}
            posvec = [ SVector(ntuple(i -> v[i], Val(length(v)))) for v in positions ]
        elseif positions isa AbstractMatrix
            posvec = [ SVector(ntuple(i -> v[i], Val(length(v)))) for v in eachcol(positions) ]
        end
        return posvec
    end
    
    function with_ParticleSystem(file, unit_cell, correct, lcell)
        coordinates = copy_to_svector(getcoor(file))
        system = ParticleSystem(
            xpositions=coordinates,
            unitcell=unit_cell,
            cutoff=10.0,
            output=0.0,
            lcell=lcell,
        )
        u = map_pairwise!(
            (pair, u) -> lj_NE(pair.d2, u),
            system,
        )
        return u ≈ correct
    end
    
    dir = @__DIR__
    
    lcell = 1
    
    # Some orthorhombic cells
    
    unit_cell = [50.0, 50.0, 50.0]
    correct = 32230.01699504111
    @test with_ParticleSystem("$dir/o1.dcd", unit_cell, correct, lcell)
    
    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    @test with_ParticleSystem("$dir/o2.dcd", unit_cell, correct, lcell)
    
    # Orthorhombic but rotated
    
    unit_cell = [50.0 0.0 50.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 1724.3195067566828
    @test with_ParticleSystem("$dir/o3.dcd", unit_cell, correct, lcell)
    
    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1754.0802503953591
    @test with_ParticleSystem("$dir/o4.dcd", unit_cell, correct, lcell)
    
    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1765.1389457850137
    @test with_ParticleSystem("$dir/o5.dcd", unit_cell, correct, lcell)
    
    unit_cell = [80.0, 80.0, 80.0]
    correct = -158.04751357760088
    @test with_ParticleSystem("$dir/o6.dcd", unit_cell, correct, lcell)
    
    # Some triclinic cells
    
    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    @test with_ParticleSystem("$dir/t1.dcd", unit_cell, correct, lcell)
    
    unit_cell = [50.0 0.0 0.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 32096.48839031735
    @test with_ParticleSystem("$dir/t2.dcd", unit_cell, correct, lcell)
    
    #
    # Check cell list updating routine
    #
    
    coordinates = getcoor("$dir/o1.dcd")
    unit_cell = [50.0 0.0 0.0
        0.0 50.0 0.0
        0.0 0.0 50.0]
    correct = 32230.01699504111
    system = ParticleSystem(
        xpositions=copy_to_svector(coordinates),
        unitcell=unit_cell,
        cutoff=10.0,
        output=0.0,
    )
    u = map_pairwise!(
        (pair, u) -> lj_NE(pair.d2, u),
        system,
    )
    @test u ≈ correct
    
    
    coordinates = getcoor("$dir/o2.dcd")
    unit_cell = [80.0 0.0 0.0; 0.0 70.0 0.0; 0.0 0.0 50.0] 
    correct = 1093.7225407797744
    
    # Update coordinates and unit cell of ParticleSystem
    system.xpositions .= copy_to_svector(coordinates)
    update_unitcell!(system, unit_cell)
    u = map_pairwise!((pair, u) -> lj_NE(pair.d2, u), system)
    @test u ≈ correct
    
    # Same thing with a different lcell
    
    lcell = 3
    
    # Some orthorhombic cells
    
    unit_cell = [50.0, 50.0, 50.0]
    correct = 32230.01699504111
    @test with_ParticleSystem("$dir/o1.dcd", unit_cell, correct, lcell)
    
    unit_cell = [80.0, 70.0, 50.0]
    correct = 1093.7225407797744
    @test with_ParticleSystem("$dir/o2.dcd", unit_cell, correct, lcell)
    
    # Orthorhombic but rotated
    
    unit_cell = [50.0 0.0 50.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 1724.3195067566828
    @test with_ParticleSystem("$dir/o3.dcd", unit_cell, correct, lcell)
    
    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1754.0802503953591
    @test with_ParticleSystem("$dir/o4.dcd", unit_cell, correct, lcell)
    
    unit_cell = transpose([70.7107 0.0 0.0
        35.3553 61.2372 0.0
        35.3553 20.4124 57.735])
    correct = 1765.1389457850137
    @test with_ParticleSystem("$dir/o5.dcd", unit_cell, correct, lcell)
    
    unit_cell = [80.0, 80.0, 80.0]
    correct = -158.04751357760088
    @test with_ParticleSystem("$dir/o6.dcd", unit_cell, correct, lcell)
    
    # Some triclinic cells
    
    unit_cell = [80.0 0.0 30.0
        30.0 80.0 0.0
        0.0 40.0 80.0]
    correct = -116.53213607052128
    @test with_ParticleSystem("$dir/t1.dcd", unit_cell, correct, lcell)
    
    unit_cell = [50.0 0.0 0.0
        50.0 50.0 0.0
        0.0 50.0 50.0]
    correct = 32096.48839031735
    @test with_ParticleSystem("$dir/t2.dcd", unit_cell, correct, lcell)
    
    #
    # Check cell list updating routine
    #
    
    coordinates = getcoor("$dir/o1.dcd")
    unit_cell = [50.0 0.0 0.0
        0.0 50.0 0.0
        0.0 0.0 50.0]
    correct = 32230.01699504111
    system = ParticleSystem(
        xpositions=copy_to_svector(coordinates),
        unitcell=unit_cell,
        cutoff=10.0,
        output=0.0,
        lcell=lcell,
    )
    u = map_pairwise!((pair, u) -> lj_NE(pair.d2, u), system)
    @test u ≈ correct
    
    coordinates = getcoor("$dir/o2.dcd")
    unit_cell = [80.0 0.0 0.0; 0.0 70.0 0.0; 0.0 0.0 50.0]
    correct = 1093.7225407797744
    system.xpositions .= copy_to_svector(coordinates)
    update_unitcell!(system, unit_cell)
    u = map_pairwise!((pair, u) -> lj_NE(pair.d2, u), system)
    @test u ≈ correct

end