# CellListMap.jl

This package is for computing interactions or any other property that is dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. It maps a function to be computed pairwise using cell lists, using periodic boundary conditions of any type. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example an average distance or an histogram of distances, forces, potentials, minimum distances, etc., as the examples below illustrate. This is done by passing the function to be evaluated as a parameter of the `map_pairwise!` function. 

## Contents

1. [Installation](#installation)
2. [Overview](#overview)
3. [Examples](#examples)
      1. [Mean difference of coordinates](#mean-difference-of-coordinates) 
      2. [Histogram of distances](#histogram-of-distances) 
      3. [Gravitational potential](#gravitational-potential) 
      4. [Gravitational force](#gravitational-force) 
      5. [Nearest neighbour](#nearest-neighbour) 
      6. [Neighbour list](#neighbour-list) 
      7. [Periodic boundary conditions](#periodic-boundary-conditions)
4. [Parallelization splitting and reduction](#parallelization-splitting-and-reduction)
      1. [Custom reduction functions](#custom-reduction-functions) 
5. [Preallocating auxiliary arrays: threaded output and cell lists](#preallocating-auxiliary-arrays-threaded-output-and-cell-lists)
      1. [Preallocating the cell lists](#preallocating-the-cell-lists)
      2. [Preallocating threaded output auxiliary arrays](#preallocating-threaded-output-auxiliary-arrays) 
6. [Performance tunning and additional options](#performance-tunning-and-additional-options)
      1. [Optimizing the cell grid](#optimizing-the-cell-grid)
      2. [Output progress](#output-progress)
7. [Some benchmarks](#some-benchmarks)
8. [Citation](#citation)

## Installation

```julia
julia> ] add CellListMap
```

## Overview

The main function is `map_parwise!`: 

If the analysis is performed on the pairs of a single vector `x` (`n*(n-1)/2` pairs), the function can be called with:
```julia
map_pairwise!(f::Function,output,box::Box,cl::CellList)
```
while if two distinct sets of points are provided (`n*m` pairs), it is called with:
```julia
map_pairwise!(f::Function,output,box::Box,cl::CellListPair)
```
where the `cl` variable of type `CellList` or `CellListPair` contains the cell lists built from the coordinates of the system, and `box` contains the system box properties.

These functions will run over every pair of particles which are closer than `box.cutoff` and compute the (squared) Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a user defined function `f` of the coordinates of the two particles will be computed. 

The function `f` receives six arguments as input: 
```julia
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it.  The squared distance `d2` is computed   internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used for the desired computation. Thus, the function `f` that is passed to `map_pairwise!` must be always of the form:
```julia
function f(x,y,i,j,d2,output)
  # update output
  return output
end
```
and the user can define more or less parameters or additional data required to compute the function using closures, as shown in the examples.

Parallel calculations are the default if more than one thread is available. Use `parallel=false` as an optional argument to `map_pairwise!` to run the serial version instead.

## Examples

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

### Histogram of distances

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

### Gravitational potential

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
forces = map_pairwise!(
  (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
  forces,box,cl
)

```

The example above can be run with `CellListMap.test4()`. 

### Nearest neighbour

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

### Neighbour list

In this example we compute the complete neighbour list, of all pairs of particles which are closer than the desired cutoff. The implementation returns a vector of tuples, in which each tuple contains the indexes of the particles and the corresponding distance. The empty `pairs` output array will be split in one vector for each thread, and reduced with a custom reduction function. 

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
  for i in 2:Threads.nthreads()
    append!(pairs,pairs_threaded[i])
  end
  return pairs
end

# Initialize output
pairs = Tuple{Int,Int,Float64}[]

# Run pairwise computation
pairs = map_pairwise!(
  (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs,cutoff),
  pairs,box,cl,
  reduce=reduce_pairs
)
```

The full example can be run with `CellListMap.test7()`. 

## Periodic boundary conditions

Triclinic periodic boundary conditions of any kind can be used. However, the input has some limitations for the moment. The lattice vectors must have strictly positive coordinates, and the smallest distance within the cell cannot be smaller than twice the size of the cutoff. An error will be produced if the cell does not satisfy these conditions. 

Let us illustrate building a two-dimensional cell, for easier visualization. A matrix of column-wise lattice vectors is provided in the construction of the box, and that is all. 

Here, the lattice vectors are `[1,0]` and `[0.5,1]` (and we illustrate with `cutoff=0.1`): 

```julia
julia> box = Box([ 1.0  0.5
                     0  1.0 ], 0.1);

julia> x = 10*rand(SVector{2,Float64},1000);
```
We have created random coordinates for `1000` particles, that are not necessarily wrapped according to the periodic boundary conditions. We can see the coordinates in the minimum image cell with:
```julia
julia> using Plots

julia> CellListMap.draw_computing_cell(x,box)
```

<img src=./src/assets/lattice.png>

The construction of the cell list is, as always, done with:

```julia
julia> cl = CellList(x,box)
CellList{2, Float64}
  109 cells with real particles.
  2041 particles in computing box, including images.

```

Upon construction of the cell lists, the particles are replicated to fill a rectangular box (or orthorhombic box, in three-dimensions), with boundaries that exceed the actual system size. This improves the performance of the pairwise computations by avoding the necessity of wrapping coordinates on the main loop (this is an implementation detail only). 

In summary, to use arbitrary periodic boundary conditions, just initialize the box with the matrix of lattice vectors. In three dimensions, for example, one could use:

```julia
julia> box = Box([ 50.  0. 00. 
                    0. 30. 30.          
                    0. 00. 50. ],  2.)

julia> x = 100*rand(SVector{3,Float64},10000);

julia> p = [ CellListMap.wrap_to_first(x,box) for x in x ];

julia> scatter(Tuple.(p),aspect_ratio=1,framestyle=:box,label=:none)
```
to work with an arbitrary 3D lattice, Which in this case looks like:

<img src=./src/assets/3Dlattice.png>

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

### Custom reduction functions

In some cases, as in the [Nearest neighbour](#nearest-neighbour) example, the output is a tuple and reduction consists in keeping the output from each thread having the minimum value for the distance. Thus, the reduction operation is not a simple sum over the elements of each threaded output. We can, therefore, overwrite the default reduction method, by passing the reduction function as the `reduce` parameter of `map_pairwise!`:
```julia
mind = map_pairwise!( 
  (x,y,i,j,d2,mind) -> f(i,j,d2,mind), mind,box,cl;
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
  return mind
end
```
This function *must* return the updated `output` variable, being it mutable or not, to be compatible with the interface.  

## Preallocating auxiliary arrays: threaded output and cell lists

### Preallocating the cell lists

The arrays containing the cell lists can be initialized only once, and then updated. This is useful for iterative runs. Note that, since the list size depends on the box size and cutoff, if the box properties changes some arrays might be increased (never shrinked) on this update. 

```julia
# Initialize cell lists with initial coordinates
cl = CellList(x,box)
for i in 1:nsteps
  x = ... # new coordinates
  box = Box(sides,cutoff) # perhaps the box has changed
  cl = UpdateCellList!(x,box) 
end
```

The procedure is identical if using two sets of coordinates, in which case, one would do:

```julia
# Initialize cell lists with initial coordinates
cl = CellList(x,y,box)
for i in 1:nsteps
  x = ... # new coordinates
  box = Box(sides,cutoff) # perhaps the box has changed
  cl = UpdateCellList!(x,y,box)
end
```

### Preallocating threaded output auxiliary arrays

On parallel runs, note that `output_threaded` is, by default, initialized on the call to `map_pairwise!`. Thus, if the calculation must be run multiple times (for example, for several steps of a trajectory), it is probably a good idea to preallocate the threaded output, particularly if it is a large array. For example, the arrays of forces should be created only once, and reset to zero after each use:
```julia
forces = zeros(SVector{3,Float64},N)
forces_threaded = [ deepcopy(forces) for i in 1:nthreads() ]
for i in 1:nsteps
  map_pairwise!(f, forces, box, cl, output_threaded=forces_threaded)
  # work with the final forces vector
  ...
  # Reset forces_threaded
  for i in 1:Threads.nthreads()
    @. forces_threaded[i] = zero(SVector{3,Float64}) 
  end
end
```
In this case, the `forces` vector will be updated by the default reduction method.

## Performance tunning and additional options

### Optimizing the cell grid

Until this is automatized (hopefuly soon), the partition of the space into cells is dependent on a parameter `lcell` which can be passed to `Box`. For example:
```julia
box = Box(x,box,lcell=2)
cl = CellList(x,box)
map_pairwise!(...)
```
This parameter determines how fine is the mesh of cells. There is a trade-off between the number of cells and the number of particles per cell. For low-density systems, greater meshes are better, because each cell will have only a few particles and the computations loop over a samller number of cells. For dense systems, it is better to run over more cells with less particles per cell. It is a good idea to test different values of `lcell` to check which is the optimal choice for your system. Usually the best value is between `lcell=1` and `lcell=6`, but for large and dense systems a larger value may be optimal. For molecular systems with normal densities `lcell=1` is likely the optimal choice. The peformance can be tested using the progress meter, as explained below.  

### Output progress 

For long-running computations, the user might want to see the progress. A progress meter can be turned on with the `show_progress` option. For example:
```julia
map_pairwise!(f,output,box,cl,show_progress=true)
```
whill print something like:
```julia-repl
Progress:  43%|█████████████████                    | ETA: 0:18:25
```

Thus, besides being useful for following the progress of a long run, it is useful to test different values of `lcell` to tune the peformance of the code, by looking at the estimated time to finish (ETA) and killing the execution after a sample run. The default and recommended option for production runs is to use `show_progress=false`, because tracking the progress introduces a small overhead into the computation. 

## Some benchmarks

The goal here is to provide a good implementation of cell lists. We compare it with the implementation of the nice cython/python [halotools](https://github.com/astropy/halotools) package, in the computation of an histogram of mean pairwise velocities. 

<img src=https://github.com/lmiq/PairVelocities/blob/main/data/cd_v0.5.3.png>

<img src=https://github.com/lmiq/PairVelocities/blob/main/data/cv_v0.5.3.png>

The full test is available [at this](https://github.com/lmiq/PairVelocities) repository. And we kindly thank [Carolina Cuesta](https://github.com/florpi) for providing the example. These benchmarks were run on an Intel i7 8th gen laptop, with 4 cores (8 threads). 

## Citation

If you use this software and need to cite it, please use the following reference:

Martínez, Leandro. (2021, June 11). CellListMap.jl: Flexible implementation of cell lists to map the calculations of short-ranged particle-pair dependent functions, such as forces, energies, neighbour lists, etc. Zenodo. http://doi.org/10.5281/zenodo.4927541






