using CellListMap
using StaticArrays
import Random

# Structure that stores a minimum distance pair
struct Mind
    i::Int
    j::Int
    d::Float64
end

# Auxiliary functions for parallelization
CellListMap.copy_output(m::Mind) = m
CellListMap.reset_output!(::Mind) = Mind(0, 0, +Inf) 
CellListMap.reducer!(m1::Mind, m2::Mind) =  m1.d <= m2.d ? m1 : m2

#
# In this test we compute the minimum distance between two independent sets of particles,
# without periodic conditions. We need to define a special reducing function to the
# parallel computation.
#
function nearest_neighbor_nopbc(; N1 = 1_500, N2 = 1_500_000, parallel = true, x = nothing, y = nothing)

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ rand(SVector{3, Float64}) for i in 1:N1 ]
    else
        N1 = length(x)
    end
    if y === nothing
        y = [ rand(SVector{3, Float64}) for i in 1:N2 ]
    else
        N2 = length(y)
    end

    # Obtain one upper bound for dmin by computing one distance for each element
    # of the smallest vector
    cutoff = +Inf
    for v in x
        iy = rand(1:N2)
        cutoff = min(CellListMap.norm(v - y[iy]), cutoff)
    end

    sys = ParticleSystem(
        xpositions=x,
        ypositions=y,
        cutoff=cutoff,
        output=Mind(0,0,+Inf),
        parallel=parallel,
    )

    # Function that keeps the minimum distance
    f(pair, mind) = pair.d < mind.d ? Mind(pair.i, pair.j, pair.d) : mind

    # Run pairwise computation
    mind = pairwise!(f, sys)

    # Take the square root of the minimum distance to return
    return (mind.i, mind.j, mind.d)
end
