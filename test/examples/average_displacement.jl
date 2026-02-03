using CellListMap
using StaticArrays
import Random

#
# In this test we compute the average displacement of the x coordinates of the particles
#              
function average_displacement(; N=100_000, parallel=true, x=nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250, 250, 250]
    cutoff = 10
    box = CellListMap.Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [sides .* rand(SVector{3,Float64}) for i in 1:N]
    end

    # Function to be evaluated for each pair: sum of displacements on x
    f(pair, avg_dx) = avg_dx + abs(pair.x[1] - pair.y[1])

    sys = ParticleSystem(
        positions=x,
        unitcell=sides,
        cutoff=cutoff,
        output=0.0,
        parallel=parallel,
    )

    avg_dx = (N / (N * (N - 1) / 2)) * pairwise!(f, sys)

    return avg_dx
end
