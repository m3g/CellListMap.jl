using CellListMap
using StaticArrays
import Random

#
# In this example we compute the "gravitational potential", assigning to each particle
# has a different mass. In this case, the closure is used to pass the masses to the
# function that computes the potential.
#
function gravitational_potential(;N=100_000,parallel=true,x=nothing)

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

    # masses
    mass = [ 5 * x[i][1] for i in 1:N ]

    # Function to be evalulated for each pair: build distance histogram
    function potential(i, j, d2, u, mass) 
        d = sqrt(d2)
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end

    # Run pairwise computation
    u = map_pairwise(
        (x, y, i, j, d2, u) -> potential(i, j, d2, u, mass),
        0.0,box,cl,
        parallel=parallel
    )

    return u

end
