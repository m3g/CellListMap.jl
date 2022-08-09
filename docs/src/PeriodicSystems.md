# PeriodicSystems interface

The `PeriodicSystems` interface facilitates the use of `CellListMap` for the majority of cases. To use it, load the `PeriodicSystems` module directly, with:

```julia
using CellListMap.PeriodicSystems
```

## Basic usage

The `PeriodicSystem` constructor receives the properties of the system and sets up automatically the most commonly used data structures necessary. For example, let us build a system of random particles in a cubic box:

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> system = PeriodicSystem(
           positions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0,
           output_name = :energy
       )
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [12, 12, 12]
      computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
      Total number of cells = 1728
    CellListMap.CellList{3, Float64}
      1000 real particles.
      646 cells with real particles.
      1648 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 12
    Type of output variable (energy): Float64
```

Now, directly, let us compute a putative energy of the particles, assuming a simple formula which depends on the inverse of the distance between pairs:

```julia-repl
julia> map_pairwise!((x,y,i,j,d2,output) -> output += 1 / (1 + sqrt(d2)), system)
1914.7714461887822
```

The `system.energy` field accesses the resulting value of the computation:
```julia-repl
julia> system.energy
1914.7714461887822
```

## Nearest neighbor

Here we compute the indexes of the particles that satisfy the minimum distance between two sets of points, using the linked lists. The distance and the indexes are stored in a tuple, and a reducing method has to be defined for that tuple to run the calculation.  The function does not need the coordinates of the points, only their distance and indexes.

```julia
# Number of particles, sides and cutoff
N1=1_500
N2=1_500_000
sides = [250,250,250]

# Particle positions
x = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N1 ]
y = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N2 ]

const MinimumDistance = Tuple{Int,Int,Float64}
PeriodicSystems.copy_output(x::MinimumDistance) = x
PeriodicSystems.reset_output(x::MinimumDistance) = (0, 0, +Inf)

system = PeriodicSystem(
    xpositions = x,
    ypositions = y,
    unitcell = sides,
    cutoff = 10.0,
    output = (0, 0, +Inf),
    output_name = :minimum_distance,
)

# Function that keeps the minimum distance
f(i,j,d2,mind) = d2 < mind[3] ? (i,j,d2) : mind

# We have to define our own reduce function here
function PeriodicSystems.reduce_output(
    output::T,
    output_threaded::Vector{T}
) where T <: MinimumDistance
    mind = (0, 0, +Inf)
    for i in eachindex(output_threaded)
        if output_threaded[i][3] < mind[3]
            mind = output_threaded[i]
        end
    end
    return mind
end

# Run pairwise computation
mind = map_pairwise((x,y,i,j,d2,mind) -> f(i,j,d2,mind), system)
```

The example above can be run with `CellListMap.Examples.nearest_neighbor()` and is available in the
[nearest_neighbor.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/nearest_neighbor.jl) file.

The example `CellListMap.Examples.nearest_neighbor_nopbc()` of [nearest\_neighbor\_nopbc.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/nearest_neighbor_nopbc.jl) describes a similar problem but *without* periodic boundary conditions. Depending on the distribution of points and size it is a faster method than usual ball-tree methods. 

## Neighbor lists

### The `CellListMap.neighborlist` function
The package provides a `neighborlist` function that implements this calculation. Without periodic boundary conditions, just do:

```julia-repl
julia> x = [ rand(3) for _ in 1:10_000 ];

julia> CellListMap.neighborlist(x,0.05)
24778-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 62, 0.028481068525796384)
 ⋮
 (9954, 1749, 0.04887502372299809)
 (9974, 124, 0.040110356034451795)
```
or `CellListMap.neighborlist(x,y,r)` for computing the lists of pairs of two sets closer than `r`.

The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

If periodic boundary conditions are used, the `Box` and `CellList` must be constructed in advance:
```julia-repl
julia> x = [ rand(3) for _ in 1:10_000 ]; 

julia> box = Box([1,1,1],0.1);

julia> cl = CellList(x,box);

julia> CellListMap.neighborlist(box,cl)

julia> CellListMap.neighborlist(box,cl)
209506-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 121, 0.05553035041478053)
 (1, 1589, 0.051415489701932444)
 ⋮
 (7469, 7946, 0.09760096646331885)
```

### Implementation

The implementation of the above function follows the principles below. 
 The empty `pairs` output array will be split in one vector for each thread, and reduced with a custom reduction function. 

```julia
# Function to be evaluated for each pair: push pair
function push_pair!(i,j,d2,pairs)
    d = sqrt(d2)
    push!(pairs,(i,j,d))
    return pairs
end

# Reduction function
function reduce_pairs(pairs,pairs_threaded)
    for i in 1:length(pairs_threaded)
        append!(pairs,pairs_threaded[i])
    end
    return pairs
end

# Initialize output
pairs = Tuple{Int,Int,Float64}[]

# Run pairwise computation
map_pairwise!(
    (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs),
    pairs,box,cl,
    reduce=reduce_pairs
)
```

The full example can be run with `CellListMap.Examples.neighborlist()`, available in the file 
[neighborlist.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/neighborlist.jl).
