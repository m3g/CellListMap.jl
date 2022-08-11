# PeriodicSystems interface

The `PeriodicSystems` interface facilitates the use of `CellListMap` for the majority of cases. To use it, load the `PeriodicSystems` module directly, with:

```julia
using CellListMap.PeriodicSystems
```

## Common output data types

The `output` of the `CellListMap` computation may be of any kind. Most commonly, it is an energy, a set of forces, or other data type that can be represented either as a number, an array of numbers, or an array of vectors (`SVectors` in particular), such as arrays of forces.  

Additionally, the combination rules for output data are assumed to be just the sum (the energy is the sum of the energy of the particles, or the forces are added by summation). 

For these types of `output` data the usage of `CellListMap.PeriodicSystems` is the simplest, and does not require the implementation of any data-type dependent function. 

### Computing the energy of a set of particles

For example, let us build a system of random particles in a cubic box, and compute an "energy", which in this case is simply the sum of `1/d` over all pair of particles, within a cutoff.

The `PeriodicSystem` constructor receives the properties of the system and sets up automatically the most commonly used data structures necessary. 

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> system = PeriodicSystem(
           positions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0,
           output_name = :energy
       );
```

Now, directly, let us compute a putative energy of the particles, assuming a simple formula which depends on the inverse of the distance between pairs:

```julia-repl
julia> map_pairwise!((x,y,i,j,d2,energy) -> energy += 1 / sqrt(d2), system)
30679.386366872823
```

The `system.energy` field accesses the resulting value of the computation:
```julia-repl
julia> system.energy
1914.7714461887822
```

### Computing forces between particles

Following the example above, let us compute the forces between the particles. We have to define the function that computes the force between a pair of particles and updates the array of forces:

```julia
function update_forces!(x,y,i,j,d2,forces)
    d = sqrt(d2)
    df = (1/d2)*(1/d)*(y - x)
    forces[i] += df
    forces[j] -= df
    return forces
end
```

Importantly, the function *must* return the `forces` array to follow the API. 

Now, let us setup the system with the new type of output variable, which will be now an array of forces with the same type as the positions:



```julia-repl
julia> positions = rand(SVector{3,Float64},1000);

julia> system = PeriodicSystem(
           positions = positions,
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = similar(positions),
           output_name = :forces
       );
```

Let us note that the `forces` where reset upon the construction of the system:
```julia-repl
julia> system.forces
1000-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 ⋮
 [0.0, 0.0, 0.0]
```

A call to `map_pairwise!` with the appropriate function definition will update the forces:
```julia-repl
julia> map_pairwise!((x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces), system)
1000-element Vector{SVector{3, Float64}}:
 [-151.19529230407284, 159.33819000196905, -261.3055111242796]
 [-173.02442398784672, -178.782819965489, 4.570607952876692]
 ⋮
 [-82.96794866090711, -635.9779270880592, -279.84420678948067]
 [-722.5400961501635, 182.65287417718935, 380.0394926753039]
```

## Updating coordinates, unit cell, and cutoff

If the `map_pairwise!` function will compute energy and/or forces in a iterative procedure (a simulation, for instance), we need to update the coordinates, and perhaps the unit cell and the cutoff.

### Updating coordinates

The coordinates can be updated (mutated, or the array of coordinates can change in size by pushing or deleting particles), simply by directly acessing the `positions` field of the system:

```julia-repl
julia> system.positions[1]
3-element SVector{3, Float64} with indices SOneTo(3):
 0.6391290709055079
 0.43679325975360894
 0.8231829019768698

julia> system.positions[1] = zeros(SVector{3,Float64})
3-element SVector{3, Float64} with indices SOneTo(3):
 0.0
 0.0
 0.0

julia> push!(system.positions, rand(SVector{3,Float64}))
1001-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.5491373098208292, 0.23899915605319244, 0.49058287555218516]
 [0.4790813256330335, 0.9922655645411307, 0.8489638318699169]
 ⋮
 [0.4700394061063937, 0.5440026379397457, 0.7411235688716618]
 [0.2973028974000733, 0.9566251992966597, 0.7427323891563248]
```

The `output` arrays may have to be resized accordingly, depending on
the calculation being performed. For example, that is the case in
the caculation of `forces`, above. In this case, we need to resize
the output arrays with the function `resize_output!`:

```julia-repl
julia> resize_output!(system, length(system.positons));

julia> map_pairwise!((x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces), system)
1001-element Vector{SVector{3, Float64}}:
 [756.2076075886971, -335.1637545330828, 541.8627090466914]
 [-173.02442398784672, -178.782819965489, 4.570607952876692]
 ⋮
 [-722.5400961501635, 182.65287417718935, 380.0394926753039]
 [20.27985502389337, -193.77607810950286, -155.28968519541544]
```

In this case, if the `output` is not resized, a `BoundsError:` is
be obtained, because updates of forces at unavailable positions will
be attempted. 

### Updating the unit cell

The unit cell can be updated to new dimensions at any moment, with the `update_unitcell!` function:

```julia-repl

```

### Updating the cutoff

The cutoff can also be updated, using the `update_cutoff!` function:

```julia-repl
julia> update_cutoff!(system, 0.2)
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.2
      number of computing cells on each dimension = [7, 7, 7]
      computing cell sizes = [0.2, 0.2, 0.2] (lcell: 1)
      Total number of cells = 343
    CellListMap.CellList{3, Float64}
      1000 real particles.
      125 cells with real particles.
      2792 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 8
    Type of output variable (forces): Vector{SVector{3, Float64}}

julia> map_pairwise!((x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces), system)
1000-element Vector{SVector{3, Float64}}:
 [306.9612911344924, -618.7375562535321, -607.1449767066479]
 [224.0803003775478, -241.05319348787023, 67.53780411933884]
 ⋮
 [2114.4873184508524, -3186.265279868732, -6777.748445712408]
 [-25.306486853608945, 119.69319481834582, 104.1501577339471]
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
