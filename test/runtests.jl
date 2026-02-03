using TestItemRunner
@run_package_tests

# Aqua
@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(CellListMap)
end

# load test modules
include("./modules/AllocTest.jl")
include("./modules/Testing.jl")
include("./modules/TestingNeighborLists.jl")

# Define testitems
include("$(@__DIR__)/API/BasicForParticleSystem.jl")
include("$(@__DIR__)/API/neighborlists.jl")
include("$(@__DIR__)/API/ParticleSystem.jl")

include("$(@__DIR__)/internals/auxiliary_functions.jl")
include("$(@__DIR__)/internals/Box.jl")
include("$(@__DIR__)/internals/CellLists.jl")
include("$(@__DIR__)/internals/CellOperations.jl")
include("$(@__DIR__)/internals/tests.jl")

include("$(@__DIR__)/applications/gromacs/compare_with_gromacs.jl")
include("$(@__DIR__)/applications/namd/compare_with_namd.jl")
include("$(@__DIR__)/applications/namd/ParticleSystem_vs_NAMD.jl")
include("$(@__DIR__)/API/test_show.jl")
