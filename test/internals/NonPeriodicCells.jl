@testitem "NonPeriodicCells"  begin
    using CellListMap
    using StaticArrays
    using ShowMethodTesting

    x = rand(SVector{3,Float64}, 90)
    y = rand(SVector{3,Float64}, 130)

    sys = ParticleSystem(positions=x, cutoff=0.1, output=0.0, nbatches=(4,4))
    @test parse_show(sys._aux; repl = Dict("CellListMap." => "")) ≈ """ 
        CellListMap.AuxNonPeriodic{3, Float64}
        Auxiliary arrays for nbatches = 4
    """

    sys2 = ParticleSystem(xpositions=x, ypositions=y, cutoff=0.1, output=0.0, nbatches=(4,4))
    @test parse_show(sys2._aux; repl = Dict("CellListMap." => "")) ≈ """ 
        CellListMap.AuxNonPeriodicPair{3, Float64}
        Auxiliary arrays for nbatches = 4
    """

end