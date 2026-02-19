using CellListMap
using StaticArrays
import Random
const PSP = ParticleSystemPositions

#
# In this test we compute the minimum distance between two independent sets of particles
#

struct Mind
    i::Int
    j::Int
    d::Float64
end
CellListMap.reducer(x::Mind, y::Mind) = x.d <= y.d ? x : y
CellListMap.copy_output(x::Mind) = x
CellListMap.reset_output(::Mind) = Mind(0,0,+Inf)
 
function nearest_neighbor(; N1 = 1_500, N2 = 1_500_000, parallel = true, x = nothing, y = nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250, 250, 250]
    cutoff = 10.0

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N1 ]
    end
    if y === nothing
        y = [ sides .* rand(SVector{3, Float64}) for i in 1:N2 ]
    end

    sys = ParticleSystem(
        xpositions=x,
        ypositions=y,
        cutoff=cutoff,
        unitcell=sides,
        output=Mind(0, 0, +Inf),
        parallel=parallel,
    )

    # Function that keeps the minimum distance
    f(pair, mind) = pair.d < mind.d ? Mind(pair.i, pair.j, pair.d) : mind

    # Run pairwise computation
    mind = pairwise!(f, sys)

    # Take the square root of the minimum distance to return
    return (mind.i, mind.j, mind.d)

end
