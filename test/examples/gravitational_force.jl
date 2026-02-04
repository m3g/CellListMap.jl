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
    box = CellListMap.Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [ sides .* rand(SVector{3, Float64}) for i in 1:N ]
    end

    # masses
    mass = [ 5 * x[i][1] for i in 1:N ]

    # Initialize auxiliary linked lists
    cl = CellListMap.CellList(x, box, parallel = parallel)

    # Function to be evaluated for each pair
    function calc_forces!(x, y, i, j, d2, mass, forces)
        G = 9.8 * mass[i] * mass[j] / d2
        d = sqrt(d2)
        df = (G / d) * (x - y)
        forces[i] = forces[i] - df
        forces[j] = forces[j] + df
        return forces
    end

    # Preallocate and initialize forces
    forces = [ zeros(SVector{3, Float64}) for i in 1:N ]

    # Run pairwise computation
    CellListMap._pairwise!(
        (pair, forces) -> calc_forces!(pair.x, pair.y, pair.i, pair.j, pair.d2, mass, forces),
        forces, box, cl,
        parallel = parallel
    )

    return forces
end
