using CellListMap
using StaticArrays
import Random

#
# In this test we compute the complete neighbor list of particles, meaning all the pairs
# that are within the cutoff distance. This function is implemented in the
#
# CellListMap.neighborlist
#
# routine. Here we illustrate how is is this implementation, to facilitate customization.
#

struct PairList
    list::Vector{Tuple{Int,Int,Float64}}
end

# Function to be evaluated for each pair: push pair if d<cutoff
function push_pair!(pair, pair_list::PairList)
    (; i, j, d) = pair
    push!(pair_list.list, (i, j, d))
    return pair_list
end

CellListMap.copy_output(x::PairList) = PairList(copy(x.list))
function CellListMap.reset_output!(x::PairList) 
    empty!(x.list)
    return x
end

# Avoid multiple resize! calls defining a custom reduce_output! directly
function CellListMap.reduce_output!(x_reduced::PairList, x_threaded::Vector{PairList})
    n = sum(length(x.list) for x in x_threaded)
    resize!(x_reduced.list, n)
    offset = 0
    for x in x_threaded
        nx = length(x.list)
        x_reduced.list[offset+1:offset+nx] .= x.list
        offset += nx
    end
    return x_reduced
end

function neighborlist(; N = 100_000, parallel = true, x = nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250, 250, 250]
    cutoff = 10.0

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N ]
    end

    sys = ParticleSystem(
        positions=x,
        cutoff=cutoff,
        unitcell=sides,
        output=PairList(Tuple{Int,Int,Float64}[]),
        parallel=parallel,
    )

    # Run pairwise computation
    CellListMap.pairwise!(push_pair!, sys)

    return sys.output.list
end
