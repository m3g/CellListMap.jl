using CellListMap
using StaticArrays
import Random

#
# In this example we compute the histogram of distances, expected to follow the
# function f(f) = ρ(4/3)π(r[i+1]^3 - r[i]^3) with ρ being the density of the system.
#
function distance_histogram(;N=100_000,parallel=true,x=nothing)

    # Number of particles, sides and cutoff
    sides = SVector{3,Float64}(250, 250, 250)
    cutoff = 10.
    box = Box(sides, cutoff)

    # Particle positions
    Random.seed!(321)
    if x === nothing 
        x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    end

    # Initialize auxiliary linked lists
    cl = CellList(x, box, parallel=parallel)

    # Function to be evalulated for each pair: build distance histogram
    function build_histogram!(d2, hist) 
        d = sqrt(d2)
        ibin = floor(Int, d) + 1
        hist[ibin] += 1
        return hist
    end

    # Preallocate and initialize histogram
    hist = zeros(Int, 10)

    # Run pairwise computation
    hist = (N / (N * (N - 1) / 2)) * map_pairwise!(
        (x, y, i, j, d2, hist) -> build_histogram!(d2, hist),
        hist,box,cl,
        parallel=parallel
    )
    return hist

end
