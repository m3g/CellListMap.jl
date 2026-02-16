using CellListMap
using StaticArrays
import Random

#
# In this example we compute the "gravitational potential", assigning to each particle
# has a different mass. In this case, the closure is used to pass the masses to the
# function that computes the potential.
#
function gravitational_potential(; N = 100_000, parallel = true, x = nothing)

    # Particle positions
    Random.seed!(321)
    sides = @SVector [250, 250, 250]
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N ]
    end
    sys = ParticleSystem(
        xpositions=x,
        unitcell=sides,
        cutoff=10.0,
        parallel=parallel,
        output=0.0,
    )

    # masses
    mass = [ 5 * x[i][1] for i in 1:N ]

    # Function to be evaluated for each pair: build distance histogram
    function potential(pair, u, mass)
        (; i, j, d) = pair
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end

    # Run pairwise computation
    u = pairwise!(
        (pair, u) -> potential(pair, u, mass),
        sys
    )

    return u

end
