using CellListMap
using StaticArrays
import Random

#
# In this test we compute the minimum distance between two independent sets of particles
#
function nearest_neighbor(;N1=1_500,N2=1_500_000,parallel=true,x=nothing,y=nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250,250,250]
    cutoff = 10.
    box = CellListMap.Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing 
        x = [ sides .* rand(SVector{3,Float64}) for i in 1:N1 ]
    end
    if y === nothing 
        y = [ sides .* rand(SVector{3,Float64}) for i in 1:N2 ]
    end

    # Initialize auxiliary linked lists
    cl = CellListMap.CellList(x, y, box, parallel=parallel)

    # Function that keeps the minimum distance
    f(i, j, d2, mind) = d2 < mind[3] ? (i, j, d2) : mind

    # We have to define our own reduce function here (for the parallel version)
    function reduce_mind(output, output_threaded)
        mind = output_threaded[1]
        for i in firstindex(output_threaded)+1:lastindex(output_threaded)
            if output_threaded[i][3] < mind[3]
                mind = output_threaded[i]
            end
        end
        return mind
    end 

    # Initialize
    mind = (0, 0, +Inf)
  
    # Run pairwise computation
    mind = pairwise!(
        (pair, mind) -> f(pair.i, pair.j, pair.d2, mind),
        mind,box,cl;reduce=reduce_mind, parallel=parallel
    )
  
    # Take the square root of the minimum distance to return
    return (mind[1], mind[2], sqrt(mind[3]))

end
