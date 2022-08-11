# PeriodicSystems interface

The `PeriodicSystems` interface facilitates the use of `CellListMap` for the majority of cases. To use it, load the `PeriodicSystems` module directly, with:

```julia
using CellListMap.PeriodicSystems
```

## The mapped function

The function to be mapped for every pair of particles within the cutoff follows the same interface as the standard interface. It must be of the form
```julia
function f(x, y, i, j, d2, output)
    # update output variable
    return output
end
```
where `x` and `y` are the positions of the particles, already wrapped relative to each other according to the periodic boundary conditions (a minimum-image set of positions), `i` and `j` are the indexes of the particles in the arrays of coordinates, `d2` is the squared distance between the particles, and `output` is the variable to be computed. 

For example, computing the energy, as the sum of the inverse of the distance between particles, can be done with a function like:
```julia
function energy(d2,u)
    u += 1 / sqrt(d2)
    return u
end
```
and the additional parameters required by the interface can be eliminated by the use of an anonymous function, directly on the call to the `map_pairwise!  function:
```julia
u = map_pairwise((x,y,i,j,d2,u) -> energy(d2,u), system)
```
(what `system` is will be explained in the examples below).

Alternatively, the function might require additional parameters, such as the masses of the particles. In this case, we can use a closure to provide such data:
```julia
function energy(i,j,d2,u,masses)
    u += masses[i]*masses[j] / sqrt(d2)
    return u
end
const masses = # ... some masses
u = map_pairwise((x,y,i,j,d2,u) -> energy(d2,u,masses), system)
```

## Potential energy example

!!! note
    The `output` of the `CellListMap` computation may be of any kind. Most commonly, it is an energy, a set of forces, or other data type that can be represented either as a number, an array of numbers, or an array of vectors (`SVectors` in particular), such as arrays of forces.  

    Additionally, the properties are frequently additive (the energy is the sum of the energy of the particles, or the forces are added by summation). 

    For these types of `output` data the usage of `CellListMap.PeriodicSystems` is the simplest, and does not require the implementation of any data-type dependent function. 

For example, let us build a system of random particles in a cubic box, and compute an "energy", which in this case is simply the sum of `1/d` over all pair of particles, within a cutoff.

The `PeriodicSystem` constructor receives the properties of the system and sets up automatically the most commonly used data structures necessary. 

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> system = PeriodicSystem(
           positions = rand(SVector{3,Float64},1000), 
           unitcell=[1.0,1.0,1.0], 
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
because the `output_name` field was provided. If it is not provided, you can access the output value from the `system.output` field.

!!! note
    - Systems can be 2 or 3-dimensional. 
    - The `unitcell` parameter may be either a vector, as in the example, or a unit cell matrix, for general boundary conditions.
    - `Unitful` quantities can be provided, given appropriate types for all input parameters. 

## Computing forces

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
 [-722.5400961501635, 182.65287417718935, 380.0394926753039]
```

## Computing both energy and forces

!!! note
    In this example we define a general type of `output` variable, for which custom copy, reset, and reduction functions
    must be defined. It can be followed for the computation of other general properties from the particle positions.

The computation of energies and forces in a single call is an interesting example for the definition of a custom `output` type and the required interface functions. 
Let us first define an output variable containing both quantities:
```julia-repl
julia> mutable struct EnergyAndForces
           energy::Float64
           forces::Vector{SVector{3,Float64}}
       end
```

Now we need to define what it means to copy, reset, and reduce this new type of output. We overload
the default corresponding functions, for our new output type:
```julia-repl
julia> import CellListMap.PeriodicSystems: copy_output, reset_output!, reduce_output!
```

The copy method creates a new instance of the `EnergyAndForces` type, with copied data:
```julia-repl
julia> copy_output(x::EnergyAndForces) = EnergyAndForces(copy(x.energy), copy(x.forces))
copy_output (generic function with 4 methods)
```

The reset method will zero both the energy and all forces:
```julia-repl
julia> function reset_output!(output::EnergyAndForces)
           output.energy = 0.0
           for i in eachindex(output.forces)
               output.forces[i] = SVector(0.0, 0.0, 0.0)
           end
           return output
       end
reset_output! (generic function with 4 methods)
```

The reduction function defines what it means to combine two output variables obtained on
independent threads. In this case, we sum the energies and forces. Different reduction functions
might be necessary for other custom types (for example if computing minimum distances).
```julia-repl
julia> function reduce_output!(output::T, output_threaded::Vector{T}) where {T <: EnergyAndForces}
           output = reset_output!(output)
           for i in eachindex(output_threaded)
               output.energy += output_threaded[i].energy
               output.forces .+= output_threaded[i].forces
           end
           return output
       end
reduce_output! (generic function with 2 methods)
```

!!! warning 
    All these functions **must** return the modified `output` variable, to adhere to the interface.

Now we can proceed as before, defining a function that updates the output variable appropriately:
```julia-repl
julia> function energy_and_forces!(x,y,i,j,d2,output::EnergyAndForces)
           d = sqrt(d2)
           output.energy += 1/d
           df = (1/d2)*(1/d)*(y - x)
           output.forces[i] += df
           output.forces[j] -= df
           return output
       end
energy_and_forces! (generic function with 1 method)
```

To finally define the system and compute the properties:

```julia-repl
julia> positions = rand(SVector{3,Float64},1000);

julia> system = PeriodicSystem(
           positions = positions,
           unitcell=[1.0,1.0,1.0], 
           cutoff = 0.1, 
           output = EnergyAndForces(0.0, similar(positions)),
           output_name = :energy_and_forces
       );

julia> system.energy_and_forces.energy
31696.94766439311

julia> system.energy_and_forces.forces
1000-element Vector{SVector{3, Float64}}:
 [-338.1909601911842, 7.7663690656924445, 202.25889647151405]
 [33.67299655756128, 282.7581453168999, -79.09639223837306]
 ⋮
 [38.83014327604529, -204.45236278342745, 249.307871211616]
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

julia> push!(system.positions, SVector(1.0, 1.0, 1.0))
1001-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.5491373098208292, 0.23899915605319244, 0.49058287555218516]
 ⋮
 [0.4700394061063937, 0.5440026379397457, 0.7411235688716618]
 [1.0, 1.0, 1.0]
```

!!! warning
    The `output` variable may have to be resized accordingly, depending on
    the calculation being performed. 

If the `output` array has to be resized, that has to be done
with the  `resize_output!` function, which will keep the consistency
of the auxiliary multi-threading buffers. This is, for instance, the case 
in the example of computation of forces, as the `forces` array must be of the
same length as the array of positions:

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
julia> update_unitcell!(system, [1.2, 1.2, 1.2])
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.2, 0.0, 0.0; 0.0, 1.2, 0.0; 0.0, 0.0, 1.2 ]
      cutoff = 0.1
      number of computing cells on each dimension = [13, 13, 13]
      computing cell sizes = [0.11, 0.11, 0.11] (lcell: 1)
      Total number of cells = 2197
    CellListMap.CellList{3, Float64}
      1000 real particles.
      623 cells with real particles.
      1719 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 12
    Type of output variable (forces): Vector{SVector{3, Float64}}
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

## Computations for two-particle sets

If the computation involves two particle sets, a similar interface is available. 
The only difference is that the coordinates of the two sets must be provided to
the `PeriodicSystem` constructor as the `xpositions` and `ypositions` arrays.

We will illustrate this interface by computing the minimum distance between two
sets of particles, which allows us to showcase further the definition of custom
type interfaces:

```julia-repl
julia> xpositions = rand(SVector{3,Float64},1000);

julia> ypositions = rand(SVector{3,Float64},1000);

julia> system = PeriodicSystem(
           xpositions = xpositions,
           ypositions = ypositions, 
           unitcell=[1.0,1.0,1.0], 
           cutoff = 0.1, 
           output = +Inf,
           output_name = :minimum_distance,
        )
```

First, we define a variable type that will carry the indexes and 
the distance of the closest pair of particles:
```julia-repl
julia> struct MinimumDistance
           i::Int
           j::Int
           d::Float64
       end
```

The function that, given two particles, retains the minimum distance, is:
```julia-repl
julia> function minimum_distance(i, j, d2, md)
           d = sqrt(d2)
           if d < md.d
               md = MinimumDistance(i, j, d)
           end
           return md
       end
minimum_distance (generic function with 1 method)
```

We overload copy, reset, and reduce functions, accordingly:
```julia-repl
julia> import CellListMap.PeriodicSystems: copy_output, reset_output!, reduce_output!

julia> copy_output(md::MinimumDistance) = md
copy_output (generic function with 5 methods)

julia> reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
reset_output! (generic function with 5 methods)

julia> function reduce_output!(md::T, md_threaded::Vector{T}) where {T <: MinimumDistance}
           md = reset_output!(md)
           for i in eachindex(md_threaded)
               if md_threaded[i].d < md.d
                   md = md_threaded[i]
               end
           end
           return md
       end
reduce_output! (generic function with 4 methods)
```

Note that since `MinimumDistance` is immutable, copying it is the same as returning the value. 
Also, resetting the minimum distance consists of setting its `d` field to `+Inf`. And, finally,
reducing the threaded distances consists of keeping the pair with the shortest distance. 

Now, we are ready to obtain the minimum distance between the sets: 
```julia-repl
julia> map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), system)
MinimumDistance(276, 617, 0.006009804808785543)
```


