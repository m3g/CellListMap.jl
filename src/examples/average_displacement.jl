using CellListMap
using StaticArrays
import Random

#
# In this test we compute the average displacement of the x coordinates of the particles
#              
function average_displacement(;N=100_000,parallel=true,x=nothing)
  
    # Number of particles, sides and cutoff
    sides = @SVector [250,250,250]
    cutoff = 10
    box = Box(sides, cutoff)
  
    # Particle positions
    Random.seed!(321)
    if x === nothing 
        x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    end
  
    # Initialize auxiliary linked lists
    cl = CellList(x, box, parallel=parallel)
  
    # Function to be evaluated for each pair: sum of displacements on x
    f(x, y, avg_dx) = avg_dx + abs(x[1] - y[1])
  
    avg_dx = (N / (N * (N - 1) / 2)) * map_pairwise(
        (x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx),
        0.,box,cl,
        parallel=parallel
    )
  
    return avg_dx
end
