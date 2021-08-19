# Examples

The full code of the examples described here is available at the [examples.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples.jl) file. 

### Mean difference of coordinates 

Computing the mean difference in `x` position between random particles. The closure is used to remove the indexes and the distance of the particles from the parameters of the input function, as they are not needed in this case.

```julia
using CellListMap

# System properties
N = 100_000
sides = [250,250,250]
cutoff = 10

# Particle positions
x = [ sides .* rand(3) for i in 1:N ]

# Initialize linked lists and box structures
box = Box(sides,cutoff)
cl = CellList(x,box)

# Function to be evaluated from positions 
f(x,y,sum_dx) = sum_dx + abs(x[1] - y[1])
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation (0.0 is the initial value)
avg_dx = normalization * map_pairwise!(
  (x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, box, cl 
)
```

The example above can be run with `CellListMap.test1()`. 

## Histogram of distances

Computing the histogram of the distances between particles (considering the same particles as in the above example). Again,
we use a closure to remove the positions and indexes of the particles from the function arguments, because they are not needed. The distance, on the other side, is needed in this example:

```julia
# Function that accumulates the histogram of distances
function build_histogram!(d2,hist)
  d = sqrt(d2)
  ibin = floor(Int,d) + 1
  hist[ibin] += 1
  return hist
end;

# Initialize (and preallocate) the histogram
hist = zeros(Int,10);
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation
hist = normalization * map_pairwise!(
  (x,y,i,j,d2,hist) -> build_histogram!(d2,hist),
  hist,box,cl
)

```

The example above can be run with `CellListMap.test2()`. 

## Gravitational potential

In this test we compute the "gravitational potential", assigning to each particle a different mass. In this case, the closure is used to pass the masses to the function that computes the potential.

```julia
# masses
const mass = rand(N)

# Function to be evaluated for each pair 
function potential(i,j,d2,mass,u)
  d = sqrt(d2)
  u = u - 9.8*mass[i]*mass[j]/d
  return u
end

# Run pairwise computation
u = map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,mass,u),0.0,box,cl)
```

The example above can be run with `CellListMap.test3()`. 

## Gravitational force

In the following example, we update a force vector of for all particles.

```julia
# masses
const mass = rand(N)

# Function to be evalulated for each pair: build distance histogram
function calc_forces!(x,y,i,j,d2,mass,forces)
  G = 9.8*mass[i]*mass[j]/d2
  d = sqrt(d2)
  df = (G/d)*(x - y)
  forces[i] = forces[i] - df
  forces[j] = forces[j] + df
  return forces
end

# Initialize and preallocate forces
forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

# Run pairwise computation
forces = map_pairwise!(
  (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
  forces,box,cl
)

```

The example above can be run with `CellListMap.test4()`. 

## Nearest neighbour

Here we compute the indexes of the particles that satisfy the minimum distance between two sets of points, using the linked lists. The distance and the indexes are stored in a tuple, and a reducing method has to be defined for that tuple to run the calculation.  The function does not need the coordinates of the points, only their distance and indexes.

```julia
# Number of particles, sides and cutoff
N1=1_500
N2=1_500_000
sides = [250,250,250]
cutoff = 10.
box = Box(sides,cutoff)

# Particle positions
x = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N1 ]
y = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N2 ]

# Initialize auxiliary linked lists (largest set!)
cl = CellList(x,y,box)

# Function that keeps the minimum distance
f(i,j,d2,mind) = d2 < mind[3] ? (i,j,d2) : mind

# We have to define our own reduce function here
function reduce_mind(output,output_threaded)
  mind = output_threaded[1]
  for i in 2:Threads.nthreads()
    if output_threaded[i][3] < mind[3]
      mind = output_threaded[i]
    end
  end
  return mind
end

# Initial value
mind = ( 0, 0, +Inf )

# Run pairwise computation
mind = map_pairwise!( 
  (x,y,i,j,d2,mind) -> f(i,j,d2,mind),
  mind,box,cl;reduce=reduce_mind
)
```

The example above can be run with `CellListMap.test5()`. The example `CellListMap.test6()` of [examples.jl](https://github.com/m3g/CellListMap.jl/blob/8661ae692abf3f44094f1fc41076c464300729b6/src/examples.jl#L219) describes a similar problem but *without* periodic boundary conditions. Depending on the distribution of points and size it is a faster method than usual ball-tree methods. 

## Neighbour list

Obs: The package provides a `neighbourlist` function that implements this calculation, and can be used with:

```julia-repl
julia> box = Box(sides,cutoff)

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   1500000 particles in the reference vector.
   1440 cells with real particles of target vector.

julia> CellListMap.neighbourlist(box,cl)
602882-element Vector{Tuple{Int64, Int64, Float64}}:
 (89, 1, 7.540726921707267)
 (1113, 9, 8.868103242326466)
 (660, 25, 6.384673918633141)
 (482, 57, 8.259716098064608)
 (189, 65, 6.251864921270615)
 ⋮
 (733, 1499960, 5.84040412313063)
 (1362, 1499968, 9.790301512741848)
 (1164, 1499968, 8.076434679245747)
 (632, 1499992, 5.3515755610344105)

```

The returning array contains tuples with the index of the particle in the first vector, the index of the particule in the second vector, and their distance.

The implementation of this function follows the principles below. 
 The empty `pairs` output array will be split in one vector for each thread, and reduced with a custom reduction function. 

```julia
# Function to be evalulated for each pair: push pair
function push_pair!(i,j,d2,pairs)
  d = sqrt(d2)
  push!(pairs,(i,j,d))
  return pairs
end

# Reduction function
function reduce_pairs(pairs,pairs_threaded)
  pairs = pairs_threaded[1]
  for i in 2:Threads.nthreads()
    append!(pairs,pairs_threaded[i])
  end
  return pairs
end

# Initialize output
pairs = Tuple{Int,Int,Float64}[]

# Run pairwise computation
pairs = map_pairwise!(
  (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs),
  pairs,box,cl,
  reduce=reduce_pairs
)
```

The full example can be run with `CellListMap.test7()`. 