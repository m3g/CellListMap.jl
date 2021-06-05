# CellListMap.jl

Maps a function to be computed pairwise using cell lists. For computing short-ranged particle interactions or other properties, considering orthorhombic (for the moment) periodic boundary conditions.

It allows the computation of any quantity from the pairs that are within the desired cutoff, for example an average distance or an histogram of distances, forces, potentials, minimum distances, etc., as the examples below illustrate. This is done by passing the function to be evaluated as a parameter of the `map_pairwise` function. 

## Installation

```julia
julia> ] add https://github.com/m3g/CellListMap.jl
```

## Examples:

```
map_pairwise(f::Function,output,x::AbstractVector,box::Box{N,T},lc::LinkedLists) where {N,T}
```

This function will run over every pair of particles which are closer than `box.cutoff` and compute
the Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a function `f` of the coordinates
of the two particles will be computed. 

This function `f` receives six arguments as input: 
```
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the   first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it. With that, it is possible to compute an average property of the     distance of the particles or, for example, build a histogram. The squared distance `d2` is computed   internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used   for the desired computation. 

### Mean difference of `x` coordinates 

Computing the mean difference in `x` position between random particles, remembering the the number of pairs of `n` particles is `n(n-1)/2`. The closure is used to remove the indexes and the distance of the atoms from the parameters of the input function, as they are not needed in this case.

```julia
using CellListMap

# System properties
n = 100_000
sides = [250,250,250]
cutoff = 10

# Initialize linked lists and box structures
lc = LinkedLists(n)
box = Box(sides,cutoff)

# Particle positions
x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:n ]

# Initialize cells (must be updated if positions change)
initcells!(x,box,lc)

# Function to be evaluated from positions 
f(x,y,sum_dx) = sum_dx + x[1] - y[1] 
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation
avg_dx = normalization * map_parwise((x,y,i,j,d2,sum_dx) -> (x,y,sum_dx), 0.0, x, box, lc)

```

The example above can be run with `CellListMap.test1()`. 

### Histogram of distances

Computing the histogram of the distances between particles (considering the same particles as in the above example). Again,
we use a closure to remove the indexes of the atoms from the computation, because they are not needed. The distance, on the other side, is needed in this example:

```julia

# Function that accumulates the histogram of distances
function build_histogram!(x,y,d2,hist)
  d = sqrt(d2)
  ibin = floor(Int,d) + 1
  hist[ibin] += 1
  return hist
end;

# Initialize (and preallocate) the histogram
hist = zeros(Int,10);
normalization = N / (N*(N-1)/2) # number of particles / number of pairs

# Run calculation
hist = normalization * map_pairwise((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),hist,x,box,lc)

```

The example above can be run with `CellListMap.test2()`. 

### Gravitational potential

In this test we compute the "gravitational potential", pretending that each particle
has a different mass. In this case, the closure is used to pass the masses to the
function that computes the potential.

```julia
# masses
mass = rand(N)

# Function to be evalulated for each pair: build distance histogram
function potential(x,y,i,j,d2,u,mass)
  d = sqrt(d2)
  u = u - 9.8*mass[i]*mass[j]/d
  return u
end

# Run pairwise computation
u = map_pairwise((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc)
```

The example above can be run with `CellListMap.test3()`. 

### Gravitational force

In the following example, we update a force vector of for all particles.

```julia

# masses
mass = rand(N)

# forces
forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

# Function to be evalulated for each pair: build distance histogram
function calc_forces!(x,y,i,j,d2,mass,forces)
  G = 9.8*mass[i]*mass[j]/d2
  d = sqrt(d2)
  df = (G/d)*(x - y)
  forces[i] = forces[i] - df
  forces[j] = forces[j] + df
  return forces
end

# Run pairwise computation
forces = map_pairwise((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),forces,x,box,lc)

```

The example above can be run with `CellListMap.test4()`. 

### Compute minimum distance between two sets of particles

Here we compute the minimum distance between two sets of points, using the linked
lists. The distance and the indexes are stored in a tuple, and a reducing method
has to be defined for that tuple to run the caculation. 

```julia
# Number of particles, sides and cutoff
N1=1_500,
N2=1_500_000
sides = [250,250,250]
cutoff = 10.
box = Box(sides,cutoff)

# Initialize auxiliary linked lists (largest set!)
lc = LinkedLists(N2)

# Particle positions
x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N1 ]
y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N2 ]

# Initializing linked cells with these positions (largest set!)
initcells!(y,box,lc)

# Function that keeps the minimum distance
f(x,y,i,j,d2,mind) = d2 < mind[3] ? mind = (i,j,d2) : mind

# We have to define our own reduce function here
function reduce(output_threaded::Vector{Tuple{Int,Int,Float64}})
  mind = output_threaded[1]
  for i in 2:nthreads()
    if output_threaded[i][3] < mind[3]
      mind = output_threaded[i]
    end
  end
  return (mind[1],mind[2],sqrt(mind[3]))
end

# Result
mind = ( 0, 0, +Inf )

# Run pairwise computation
mind = map_pairwise(f,mind,x,y,box,lc;reduce=reduce)
```

The example above can be run with `CellList.test5()`. 


