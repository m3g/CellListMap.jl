#
# Testing show methods
#
@testitem "show methods" begin
    using CellListMap
    using StaticArrays
    using ShowMethodTesting

    x = rand(SVector{3, Float64}, 100)
    y = rand(SVector{3, Float64}, 10)
    box = CellListMap.Box([10, 10, 10], 1)

    @test parse_show(box; repl = Dict("CellListMap." => "")) ≈
        """
            Box{CellListMap.OrthorhombicCell, 3}
            unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
            cutoff = 1.0
            number of computing cells on each dimension = [13, 13, 13]
            computing cell sizes = [1.0, 1.0, 1.0] (lcell: 1)
            Total number of cells = 2197    
        """

    cl = CellListMap.CellList(CellListMap.ParticleSystemPositions(x), box; nbatches = (2, 4))
    @test parse_show(cl; repl = Dict("CellListMap." => "")) ≈ """
            CellList{3, Float64}
            100 real particles.
            1 cells with real particles.
            800 particles in computing box, including images.
        """

    @test parse_show(CellListMap.nbatches(cl)) ≈ "(2, 4)"
    @test parse_show(CellListMap.AuxThreaded(cl); repl = Dict("CellListMap." => "")) ≈ 
        """
        CellListMap.AuxThreaded{3, Float64}
        Auxiliary arrays for nbatches = 2
        """

    PSP = CellListMap.ParticleSystemPositions
    cl = CellListMap.CellList(PSP(x), PSP(y), box; nbatches = (2, 4))
    @test parse_show(cl; repl = Dict("CellListMap." => "")) ≈  
        """
            CellListMap.CellListPair{3, Float64}
            1 cells with real particles in reference set.
            1 cells with real particles in target set.
        """

    @test parse_show(CellListMap.nbatches(cl); repl = Dict("CellListMap." => "")) ≈ "(2, 4)"
    @test parse_show(CellListMap.AuxThreaded(cl); repl = Dict("CellListMap." => "")) ≈ 
        """
        CellListMap.AuxThreadedPair{3, Float64}
        Auxiliary arrays for nbatches = 2
        """

    s = ParticleSystem(xpositions = x, cutoff = 0.1, unitcell = [10, 10, 10], output = 0.0, nbatches = (2, 2))
    @test isapprox(parse_show(s; repl = Dict("CellListMap." => "")),
        """ 
            ParticleSystem1{output_name} of dimension 3, composed of:
                Box{CellListMap.OrthorhombicCell, 3}
                  unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
                  cutoff = 0.1
                  number of computing cells on each dimension = [103, 103, 103]
                  computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
                  Total number of cells = 1092727
                CellList{3, Float64}
                  100 real particles.
                  95 cells with real particles.
                  139 particles in computing box, including images.
                Parallelization auxiliary data set for 2 batch(es).
                Type of output variable (output_name): Float64
        """; int_match = (x1, x2) -> isapprox(x1, x2; rtol = 1))

    s = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [10, 10, 10], output = 0.0, nbatches = (2, 2))
    @test isapprox(parse_show(s; repl = Dict("CellListMap." => "")),
        """
            ParticleSystem2{output_name} of dimension 3, composed of:
            Box{CellListMap.OrthorhombicCell, 3}
            unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
            cutoff = 0.1
            number of computing cells on each dimension = [103, 103, 103]
            computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
            Total number of cells = 1092727
            CellListMap.CellListPair{3, Float64}
            10 cells with real particles in reference set.
            97 cells with real particles in target set.
            Parallelization auxiliary data set for 2 batch(es).
            Type of output variable (output_name): Float64
        """; int_match = (x1, x2) -> isapprox(x1, x2; rtol = 1))

    @test parse_show(InPlaceNeighborList(x = x, cutoff = 0.1, unitcell = [1, 1, 1]); repl = Dict("CellListMap." => "")) ≈  
        """
            InPlaceNeighborList with types: 
            CellList{3, Float64}
            Box{CellListMap.OrthorhombicCell, 3, Float64, Float64, 9, Float64}
            Current list buffer size: 0
        """

    @test parse_show(InPlaceNeighborList(x = x, y = x, cutoff = 0.1, unitcell = [1, 1, 1]); repl = Dict("CellListMap." => "")) ≈
        """
            InPlaceNeighborList with types: 
            CellListMap.CellListPair{3, Float64}
            Box{CellListMap.OrthorhombicCell, 3, Float64, Float64, 9, Float64}
            Current list buffer size: 0
        """

end
