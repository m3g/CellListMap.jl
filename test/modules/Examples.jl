@testmodule Examples begin

using CellListMap

include("$(@__DIR__)/../examples/average_displacement.jl")
include("$(@__DIR__)/../examples/distance_histogram.jl")
include("$(@__DIR__)/../examples/gravitational_force.jl")
include("$(@__DIR__)/../examples/gravitational_potential.jl")
include("$(@__DIR__)/../examples/nearest_neighbor.jl")
include("$(@__DIR__)/../examples/nearest_neighbor_nopbc.jl")
include("$(@__DIR__)/../examples/neighborlist.jl")
include("$(@__DIR__)/../examples/pairwise_velocities.jl")
include("$(@__DIR__)/../examples/packmol.jl")
include("$(@__DIR__)/../examples/generic_types.jl")

end