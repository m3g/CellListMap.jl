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
function neighborlist(; N = 100_000, parallel = true, x = nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250, 250, 250]
    cutoff = 10.0
    box = CellListMap.Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N ]
    end

    # Initialize auxiliary linked lists
    cl = CellListMap.CellList(PSP(x), box, parallel = parallel)

    # Function to be evaluated for each pair: push pair if d<cutoff
    function push_pair!(i, j, d2, pairs, cutoff)
        d = sqrt(d2)
        if d < cutoff
            push!(pairs, (i, j, d))
        end
        return pairs
    end

    # Reduction function
    function reduce_pairs(pairs, pairs_threaded)
        for i in eachindex(pairs_threaded)
            append!(pairs, pairs_threaded[i])
        end
        return pairs
    end

    # Initialize empty output
    pairs = Tuple{Int, Int, Float64}[]

    # Run pairwise computation
    CellListMap._pairwise!(
        (pair, pairs) -> push_pair!(pair.i, pair.j, pair.d2, pairs, cutoff),
        pairs, box, cl,
        reduce = reduce_pairs,
        parallel = parallel
    )

    return pairs
end
