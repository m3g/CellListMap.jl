#
# Testing show methods
#
@testitem "show methods" begin
    using CellListMap, StaticArrays
    using CellListMap: Box, CellList, nbatches

    function test_show(
        x, s::String; 
        f64 = (x1,x2) -> isapprox(x1,x2,rtol=1e-3),
        i64 = (x1,x2) -> x1 == x2, 
        rep = Dict{String,String}(),
    )
        match(f,x1,x2) = begin
            if !f(x1,x2)
                println("show method test failed with $x1 == $x2")
                return false
            end
            return true
        end
        buff = IOBuffer()
        show(buff, MIME"text/plain"(), x)
        ss = String(take!(buff))
        s = replace(s, rep...)
        ss = replace(ss, rep...)
        sfields = split(s)
        ssfields = split(ss)
        all_match = true
        for (f1, f2) in zip(sfields, ssfields)
            !all_match && break
            value = tryparse(Int, f1)
            if !isnothing(value)
                all_match = match(i64, value, tryparse(Int, f2))
                continue
            end
            value = tryparse(Float64, f1)
            if !isnothing(value)
                all_match = match(f64, value, tryparse(Float64,f2))
                continue
            end
            all_match = match(isequal, f1, f2)
        end
        return all_match
    end

    x = rand(SVector{3,Float64}, 100)
    y = rand(SVector{3,Float64}, 10)
    box = Box([10,10,10], 1)

    @test test_show(box, """
        Box{OrthorhombicCell, 3}
        unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
        cutoff = 1.0
        number of computing cells on each dimension = [13, 13, 13]
        computing cell sizes = [1.0, 1.0, 1.0] (lcell: 1)
        Total number of cells = 2197    
    """;
        rep = Dict("CellListMap." => "")
    )

    cl = CellList(x, box; nbatches=(2,4))    
    @test test_show(cl,"""
        CellList{3, Float64}
        100 real particles.
        1 cells with real particles.
        800 particles in computing box, including images.
    """;
        rep = Dict("CellListMap." => "")
    )

    @test test_show(nbatches(cl), "(2, 4)")
    @test test_show(CellListMap.AuxThreaded(cl),"""
    CellListMap.AuxThreaded{3, Float64}
    Auxiliary arrays for nbatches = 2
    """;
        rep = Dict("CellListMap." => "")
    )

    cl = CellList(x,y,box; nbatches=(2,4))
    @test test_show(cl,"""
        CellListMap.CellListPair{3, Float64}
        1 cells with real particles of the smallest set.
        1 cells with real particles of the largest set.
    """;
        rep = Dict("CellListMap." => "")
    )

    @test test_show(nbatches(cl), "(2, 4)")
    @test test_show(CellListMap.AuxThreaded(cl),"""
    CellListMap.AuxThreadedPair{3, Float64}
    Auxiliary arrays for nbatches = 2
    """;
        rep = Dict("CellListMap." => "")
    )

    s = ParticleSystem(xpositions=x,cutoff=0.1,unitcell=[10,10,10],output=0.0,nbatches=(2,2))
    @test test_show(s, """ 
        ParticleSystem1{output} of dimension 3, composed of:
            Box{OrthorhombicCell, 3}
              unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
              cutoff = 0.1
              number of computing cells on each dimension = [103, 103, 103]
              computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
              Total number of cells = 1092727
            CellList{3, Float64}
              100 real particles.
              97 cells with real particles.
              139 particles in computing box, including images.
            Parallelization auxiliary data set for 2 batch(es).
            Type of output variable (output): Float64
    """;
        i64 = (x1,x2) -> isapprox(x1, x2, rtol=1),
        rep = Dict("CellListMap." => "")
    )

    s = ParticleSystem(xpositions=x,ypositions=y,cutoff=0.1,unitcell=[10,10,10],output=0.0,nbatches=(2,2))
    @test test_show(s, """
        ParticleSystem2{output} of dimension 3, composed of:
        Box{OrthorhombicCell, 3}
        unit cell matrix = [ 10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0 ]
        cutoff = 0.1
        number of computing cells on each dimension = [103, 103, 103]
        computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
        Total number of cells = 1092727
        CellListMap.CellListPair{3, Float64}
        10 cells with real particles of the smallest set.
        97 cells with real particles of the largest set.
        Parallelization auxiliary data set for 2 batch(es).
        Type of output variable (output): Float64
    """;
        i64 = (x1,x2) -> isapprox(x1, x2, rtol=1),
        rep = Dict("CellListMap." => "")
    )

    @test test_show(InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1]),"""
        InPlaceNeighborList with types: 
        CellList{3, Float64}
        Box{OrthorhombicCell, 3, Float64, Float64, 9, Float64}
        Current list buffer size: 0
    """;
        rep = Dict("CellListMap." => "")
    )

    @test test_show(InPlaceNeighborList(x=x, y=x, cutoff=0.1, unitcell=[1,1,1]),"""
        InPlaceNeighborList with types: 
        CellListMap.CellListPair{3, Float64}
        Box{OrthorhombicCell, 3, Float64, Float64, 9, Float64}
        Current list buffer size: 0
    """;
        rep = Dict("CellListMap." => "")
    )

end