using CellListMap
using StaticArrays
import Random

#
# In this example we compute the histogram of distances, expected to follow the
# function f(f) = ρ(4/3)π(r[i+1]^3 - r[i]^3) with ρ being the density of the system.
#
function distance_histogram(; N=100_000, parallel=true, x=nothing)

    # Number of particles, sides and cutoff
    sides = SVector{3,Float64}(250, 250, 250)
    cutoff = 10.

    # Particle positions
    Random.seed!(321)
    if x === nothing
        x = [sides .* rand(SVector{3,Float64}) for i in 1:N]
    end

    # Function to be evaluated for each pair: build distance histogram
    function build_histogram!(pair, hist)
        ibin = floor(Int, pair.d) + 1
        hist[ibin] += 1
        return hist
    end

    # output histogram
    hist = zeros(Int, 10)

    sys = ParticleSystem(
        positions=x,
        cutoff=cutoff,
        unitcell=sides,
        parallel=parallel,
        output=hist,
        output_name=:hist,
    )

    # Run pairwise computation
    pairwise!(build_histogram!, sys)

    # Normalization (convert to float)
    hist = (N / (N * (N - 1) / 2)) * hist

    return hist
end
