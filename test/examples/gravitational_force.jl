using CellListMap
using StaticArrays
import Random

#
# In this test we compute the "gravitational force", assigning to each particle
# a different mass. In this case, the closure is used to pass the masses and
# the force vector to the function that computes the potential.
#
function gravitational_force(; N = 100_000, parallel = true, x = nothing)

    # Number of particles, sides and cutoff
    sides = @SVector [250, 250, 250]
    cutoff = 10.0

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N ]
    end

    # masses
    mass = [ 5 * x[i][1] for i in 1:N ]

    # Preallocate and initialize forces
    forces = [ zeros(SVector{3, Float64}) for i in 1:N ]

    sys = ParticleSystem(
        positions = x,
        unitcell = sides,
        cutoff = cutoff,
        parallel = parallel,
        output = forces,
    )

    # Run pairwise computation
    forces = pairwise!(
        (pair, forces) -> begin
            G = 9.8 * mass[pair.i] * mass[pair.j] / pair.d2
            df = (G / pair.d) * (pair.x - pair.y)
            forces[pair.i] = forces[pair.i] - df
            forces[pair.j] = forces[pair.j] + df
            return forces
        end,
        sys
    )

    return forces
end
