# CellListMap.jl

This package is for computing short-ranged particle interactions or any other property that is dependent on the distances between pairs of particles, within a cutoff. It maps a function to be computed pairwise using cell lists using, for the moment, orthorhombic periodic boundary conditions. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example an average distance or an histogram of distances, forces, potentials, minimum distances, etc., as the examples below illustrate. This is done by passing the function to be evaluated as a parameter of the `map_pairwise!` function. 

## Contents

1. [Installation](#installation)
2. [Overview](#overview)
3. [Examples](#examples)
      1. [Mean difference of coordinates](#mean-difference-of-coordinates) 
      2. [Histogram of distances](#histogram-of-distances) 
      3. [Gravitational potential](#gravitational-potential) 
      4. [Gravitational force](#gravitational-force) 
      5. [Nearest neighbor](#nearest-neighbor) 
      6. [Neighbor list](#neighbor-list) 
4. [Parallelization splitting and reduction](#parallelization-splitting-and-reduction)
      1. [Preallocating auxiliary arrays](#preallocating-auxiliary-arrays) 
      2. [Custom reduction functions](#custom-reduction-functions) 

## Installation

```julia
julia> ] add https://github.com/m3g/CellListMap.jl
```

## Overview

The main function is `map_parwise!`: 

```julia
map_pairwise!(f::Function,output,x::AbstractVector,box::Box{N,T},lc::LinkedLists) where {N,T}
```

This function will run over every pair of particles which are closer than `box.cutoff` and compute
the (squared) Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the (squared) cutoff, a function `f` of the coordinates
of the two particles will be computed. 

The function `f` receives six arguments as input: 
```julia
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it. With that, it is possible to compute an average property of the distance of the particles or, for example, build a histogram. The squared distance `d2` is computed   internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used for the desired computation. Thus, the function `f` that is passed to `map_pairwise!` must be always of the form:
```julia
function f(x,y,i,j,d2,output)
  # update output
  return output
end
```
and the user can define more or less parameters using closures, as shown in the examples.

Parallel calculations are the default if more than one thread is available. Use `parallel=false` as an optional argument to `map_pairwise!` to run the serial version instead.

## Examples

The full code of the examples described here is available at the [examples.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples.jl) file. 

### Mean difference of coordinates 

Computing the mean difference in `x` position between random particles. The closure is used to remove the indexes and the distance of the atoms from the parameters of the input function, as they are not needed in this case.

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
initlists!(x,box,lc)

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
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation
hist = normalization * map_pairwise!((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),hist,x,box,lc)

```

The example above can be run with `CellListMap.test2()`. 

### Gravitational potential

In this test we compute the "gravitational potential", assigning to each particle a different mass. In this case, the closure is used to pass the masses to the function that computes the potential.

```julia
# masses
const mass = rand(N)

# Function to be evalulated for each pair 
function potential(x,y,i,j,d2,u,mass)
  d = sqrt(d2)
  u = u - 9.8*mass[i]*mass[j]/d
  return u
end

# Run pairwise computation
u = map_pairwise!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc)
```

The example above can be run with `CellListMap.test3()`. 

### Gravitational force

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
forces = map_pairwise!((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),forces,x,box,lc)

```

The example above can be run with `CellListMap.test4()`. 

### Nearest neighbor

Here we compute the indexes of the atoms that satisfy the minimum distance between two sets of points, using the linked lists. The distance and the indexes are stored in a tuple, and a reducing method has to be defined for that tuple to run the calculation.  The function does not need the coordinates of the points, only their distance and indexes.

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
initlists!(y,box,lc)

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
  return (mind[1],mind[2],mind[3])
end

# Initial value
mind = ( 0, 0, +Inf )

# Run pairwise computation
mind = map_pairwise!( 
  (x,y,i,j,d2,mind) -> f(i,j,d2,mind),
  mind,x,y,box,lc;reduce=reduce_mind
)
```

The example above can be run with `CellListMap.test5()`. The example `CellListMap.test6()` of [examples.jl](https://github.com/m3g/CellListMap.jl/blob/8661ae692abf3f44094f1fc41076c464300729b6/src/examples.jl#L219) describes a similar problem but *without* periodic boundary conditions. Depending on the distribution of points it is a faster method than usual ball-tree methods. 

### Neighbor list

In this example we compute the complete neighbor list, of all pairs of atoms which are closer than the desired cutoff. The implementation returns a vector of tuples, in which each tuple contains the indexes of the atoms and the corresponding distance. The empty `pairs` output array will be split in one vector for each thread, and reduced with a custom reduction function. 

```julia
# Function to be evalulated for each pair: push pair if d<cutoff
function push_pair!(i,j,d2,pairs,cutoff)
  d = sqrt(d2)
  if d < cutoff
    push!(pairs,(i,j,d))
  end
  return pairs
end

# Reduction function
function reduce_pairs(pairs,pairs_threaded)
  pairs = pairs_threaded[1]
  for i in 2:nthreads()
    append!(pairs,pairs_threaded[i])
  end
  return pairs
end

# Initialize output
pairs = Tuple{Int,Int,Float64}[]

# Run pairwise computation
pairs = map_pairwise!(
  (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs,cutoff),
  pairs,x,box,lc,
  reduce=reduce_pairs,
  parallel=parallel
)
```

The full example can be run with `CellListMap.test7()`. 

## Parallelization splitting and reduction

The parallel execution requires the splitting of the computation among threads, obviously. Thus, the output variable must be split and then reduced to avoid concurrency. To control these steps, set manually the `output_threaded` and `reduce` optional input parameters of the `map_pairwise!` function. 

By default, we define (here `using Base.Threads`):
```julia
output_threaded = [ deepcopy(output) for i in 1:nthreads() ]
```
and, for scalars and vectors, the reduction is just the sum of the output per thread:
```julia
reduce(output::Number,output_threaded) = sum(output_threaded)
function reduce(output::Vector,output_threaded) 
  for i in 1:nthreads()
     @. output += output_threaded[i] 
  end
  return output
end
```

### Preallocating auxiliary arrays

Note, however, that `output_threaded` is defined on the call to `map_pairwise!`. Thus, if the calculation must be run multiple times (for example, for several steps of a trajectory), it is probably a good idea to preallocate the threaded output, particularly if it is a large array. For example, the arrays of forces should be created only once, and reset to zero after each use:
```julia
forces = zeros(SVector{3,Float64},N)
forces_threaded = [ deepcopy(forces) for i in 1:nthreads() ]
for i in 1:nsteps
  map_pairwise!(f, forces, x, box, lc, output_threaded=forces_threaded)
  # work with the final forces vector
  ...
  # Reset forces_threaded
  for i in 1:nthreads()
    @. forces_threaded[i] = zero(SVector{3,Float64}) 
  end
end
```
In this case, the `forces` vector will be updated by the default reduction method.

### Custom reduction functions

In some cases, as in the [Nearest neighbor](#nearest-neighbor) example, the output is a tuple and reduction consists in keeping the output from each thread having the minimum value for the distance. Thus, the reduction operation is not a simple sum over the elements of each threaded output. We can, therefore, overwrite the default reduction method, by passing the reduction function as the `reduce` parameter of `map_pairwise!`:
```julia
mind = map_pairwise!( 
  (x,y,i,j,d2,mind) -> f(i,j,d2,mind), mind,x,y,box,lc;
  reduce=reduce_mind
)
```
where here the `reduce` function is set to be the custom function that keeps the tuple associated to the minimum distance obtained between threads:
```julia
function reduce_mind(output,output_threaded)
  mind = output_threaded[1]
  for i in 2:Threads.nthreads()
    if output_threaded[i][3] < mind[3]
      mind = output_threaded[i]
    end
  end
  return (mind[1],mind[2],mind[3])
end
```
This function *must* return the updated `output` variable, being it mutable or not, to be compatible with the interface.  







