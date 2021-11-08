using CellListMap
using StaticArrays
import Random

#
# In this test we compute the complete neighbour list of particles, meaning all the pairs
# that are within the cutoff distance. This function is implemented in the 
# 
# CellListMap.neighbourlist
#
# routine. Here we illustrate how is is this implementation, to facilitate customization.
#
function neighbourlist(;N=100_000,parallel=true,x=nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250,250,250]
    cutoff = 10.
    box = Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing 
        x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    end

    # Initialize auxiliary linked lists
    cl = CellList(x, box, parallel=parallel)

    # Function to be evalulated for each pair: push pair if d<cutoff
    function push_pair!(i, j, d2, pairs, cutoff) 
        d = sqrt(d2)
        if d < cutoff
            push!(pairs, (i, j, d))
        end
        return pairs
    end

    # Reduction function
    function reduce_pairs(pairs, pairs_threaded)
        pairs = pairs_threaded[1]
        for i in 2:Threads.nthreads()
            append!(pairs, pairs_threaded[i])
        end
        return pairs
    end

    # Initialize output
    pairs = Tuple{Int,Int,Float64}[]

    # Run pairwise computation
    pairs = map_pairwise!(
        (x, y, i, j, d2, pairs) -> push_pair!(i, j, d2, pairs, cutoff),
        pairs,box,cl,
        reduce=reduce_pairs,
        parallel=parallel
    )

    return pairs
end
