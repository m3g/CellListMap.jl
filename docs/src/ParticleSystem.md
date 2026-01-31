# ParticleSystem interface

The `ParticleSystem` interface facilitates the use of `CellListMap` for the majority of cases. 

!!! note
    - This interface requires `CellListMap.jl` version `0.8.30` or greater.
    - The complete codes of the examples are at the end of this page, with examples of:
        - [Simple energy computation](@ref)
        - [Force computation](@ref)
        - [Energy and forces](@ref)
        - [Two sets of particles](@ref)
        - [Particle simulation](@ref)

!!! compat
    The `ParticleSystem` interface is available since version `0.9.0` of CellListMap.jl.
    It replaces the `PeriodicSystems` interface available in previous versions.

## The mapped function

The purpose of CellListMap is to compute a pairwise-dependent function for all pairs of particles that are closer to
each other than a defined cutoff. This pairwise function must be implemented by the user and adhere to the following 
interface:  

```julia
function f(x, y, i, j, d2, output)
    # update output variable
    return output
end
```
where `x` and `y` are the positions of the particles, already wrapped relative to each other according to the periodic boundary conditions (a minimum-image set of positions), `i` and `j` are the indexes of the particles in the arrays of coordinates, `d2` is the squared distance between the particles, and `output` is the variable to be computed. 

!!! info
    #### Details of the mapped function interface

    The input parameters `x`, `y`, `i`, `j`, and `d2` must not be modified by the user. They are the
    the input data that the user may use to update the `output` variable.  

    |  Input Parameter |  Type   |  Meaning   |
    |:----------------:|:-------:|:-----------|
    | `x`              |  `SVector`   |  The coordinates of particle `i` of the pair.  |
    | `y`              |  `SVector`   |  The coordinates of particle `j` of the pair (minimum-image relative to `x`).  |
    | `i`              |  `Int`       |  Index of first particle in the original array of coordinates. |
    | `j`              |  `Int`       |  Index of second particle in the original array of coordinates. |
    | `d2`             |  `<:Real`    |  Squared distance between the particles. |
    | `output`         |  user defined   |  the value to be updated  |

    **Notes:** `x` and `y` may be 2D or 3D vectors, depending on the dimension of the system. The type of 
    the coordinates of `x`, `y`, and of `d2` are dependent on the input arrays and cutoff, and can be `Float64`,
    `Float32`, unitful quantities, etc. 

    |  Return value |  Type   |  Meaning   |
    |:----------------:|:-------:|:-----------|
    | `output`         |  user defined   |  the updated value of output.  |

    The `output` variable **must** be returned by the function, being it mutable or immutable. 


### Basic examples

For example, computing the energy, as the sum of the inverse of the distance between particles, can be done with a function like:
```julia
function energy(d2,u)
    u += 1 / sqrt(d2)
    return u
end
```
and the additional parameters required by the interface can be eliminated by the use of an anonymous function, directly on the call to the `map_pairwise` function:
```julia
u = map_pairwise(
    (x,y,i,j,d2,u) -> energy(d2,u), 
    system
)
```
(what `system` is will be explained in the examples below). Note that the `energy` function does not use the `x`, `y`, `i`, and `j` input parameters, such 
that the anonymous function managing the interface could also be written as `(_, _, _, _, d2, u) -> energy(d2, u)`, making explicit the dummy character of
these variables in the example.

Alternatively, the function might require additional parameters, such as the masses of the particles. In this case, we can use a closure to provide such data:
```julia
function energy(i,j,d2,u,masses)
    u += masses[i]*masses[j] / sqrt(d2)
    return u
end
const masses = # ... some masses
u = map_pairwise((x,y,i,j,d2,u) -> energy(d2,u,masses), system)
```

Here we reinforce the fact that the `energy` functions defined above compute the contribution to the energy of the interaction of *a single* pair
of particles. This function will be called for every pair of particles within the cutoff, automatically, in the `map_pairwise` call. 

!!! note
    The `output` of the `CellListMap` computation may be of any kind. Most commonly, it is an energy, a set of forces, or other data type that can be represented either as a number, an array of numbers, or an array of vectors (`SVectors` in particular), such as an arrays of forces.  

    Additionally, the properties are frequently additive (the energy is the sum of the energy of the particles, or the forces are added by summation). 

    For these types of `output` data the usage does not require the implementation of any data-type dependent function. 

## The ParticleSystem constructor

## Potential energy example

For example, here we read the coordinates of Argon atoms from a PDB file. The coordinates are given as 
vector of `SVector`s. We then compute an "energy", which in this case is simply the sum of `1/d` over all pair of particles, within a cutoff.

The `ParticleSystem` constructor receives the properties of the system and sets up automatically the most commonly used data structures necessary. 

```julia-repl
julia> using CellListMap, PDBTools

julia> argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

julia> system = ParticleSystem(
           xpositions=argon_coordinates,
           unitcell=[21.0,21.0,21.0], 
           cutoff = 8.0, 
           output = 0.0,
           output_name = :energy
       );
```

!!! note
    - Systems can be 2 or 3-dimensional. 
    - The `unitcell` parameter may be:
        - a vector, in which case the system periodic boundaries are Orthorhombic, this is faster.
        - a matrix, in which case the system periodic boundaries are Triclinic (general). The lattice
          vectors correspond to the *columns* of the matrix.
        - `nothing` (by default), in which case no periodic boundary conditions will be used.
    - `Unitful` quantities can be provided, given appropriate types for all input parameters. 

Now, let us compute the energy of the particles, assuming a simple formula which depends on the inverse of the distance between pairs:

```julia-repl
julia> function energy(x, y, i, j, d2, energy)
           energy += 1 / sqrt(d2)
           return energy
       end

julia> map_pairwise!(energy, system)
207.37593043370865
```
Note that the first four parameters of `energy` are not used here but are needed to adhere to the interface. The function 
input could be written as `(_, _, _, _, d2, energy)` to make that explicit. 

Because `output_name` was set to `:energy`, the `system.energy` field accesses the resulting value of the computation:
```julia-repl
julia> system.energy
207.37593043370865
```
If the `output_name` field is not provided, the output value from the `system.output` field.

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
julia> using CellListMap, PDBTools

julia> argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

julia> system = ParticleSystem(
           xpositions=argon_coordinates,
           unitcell=[21.0, 21.0, 21.0], 
           cutoff = 8.0, 
           output = similar(argon_coordinates),
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
100-element Vector{SVector{3, Float64}}:
 [0.026493833307357332, 0.18454277989323772, -0.012253902366284965]
 [0.07782602581235695, 0.2791082233740261, 0.21926615329195248]
 ⋮
 [0.11307234751448932, 0.006353545239676281, -0.05955687310348302]
 [-0.03101200918307673, 0.03543655648545697, 0.031849121630976335]
```

## Computing both energy and forces

In this example we define a general type of `output` variable, for which custom copy, reset, and reduction functions
must be defined. It can be followed for the computation of other general properties from the particle positions.

!!! note
    Interface to be implemented:

    |   Method              |   Return    |  What it does    |
    |:----------------------|:-----------:|:---------------- |
    | `copy_output(x::T)`   | new instance of type `T` | Copies an element of the output type `T`. |
    | `reset_output!(x::T)` | mutated `x` | Resets (usually zero) the value of x to the initial value it must assume before mapping.  If `x` is immutable, the function can return a new instance of `T`. |
    | `reducer(x::T,y::T)` | `mutated x` | Reduces `x` and `y` into `x` (for example `x = x + y`). If `x` is immutable, returns a new instance of type `T`.    

    **Remark:** if the output is an array of an immutable type `T`, the methods above can be defined for single *instances* of `T`, which is simpler
    than for the arrays.

```julia
using CellListMap, StaticArrays, PDBTools
```

The computation of energies and forces in a single call is an interesting example for the definition of a custom `output` type and the required interface functions. 
Let us first define an output variable containing both quantities:
```julia
mutable struct EnergyAndForces
    energy::Float64
    forces::Vector{SVector{3,Float64}}
end
```

Now we need to define what it means to copy, reset, and reduce this new type of output. We overload
the default corresponding functions, for our new output type:

The copy method creates a new instance of the `EnergyAndForces` type, with copied data:
```julia
function CellListMap.copy_output(x::EnergyAndForces) 
    return EnergyAndForces(copy(x.energy), copy(x.forces))
end
```

The reset method will zero both the energy and all forces:
```julia
function CellListMap.reset_output!(output::EnergyAndForces)
    output.energy = 0.0
    for i in eachindex(output.forces)
        output.forces[i] = SVector(0.0, 0.0, 0.0)
    end
    return output
end
```

The reducer function defines what it means to combine two output variables obtained on
independent threads. In this case, we sum the energies and forces. Different reduction functions
might be necessary for other custom types (for example if computing minimum distances).
```julia
function CellListMap.reducer(x::EnergyAndForces, y::EnergyAndForces)
    e_tot = x.energy + y.energy
    x.forces .+= y.forces
    return EnergyAndForces(e_tot, x.forces)
end
```
Note that in the above example, we reuse the `x.forces` array in the return instance
of `EnergyAndForces`. You must always reduce from right to left, and reuse the
possible buffers of the first argument of the reducer (in this case, `x`).

!!! warning 
    - All these functions **must** return the modified `output` variable, to adhere to the interface.
    - The proper definition of a reduction function is crucial for correctness. Please verify
      your results if using the default reducer function, which sums the elements.

Now we can proceed as before, defining a function that updates the output variable appropriately:
```julia
function energy_and_forces!(x,y,i,j,d2,output::EnergyAndForces)
    d = sqrt(d2)
    output.energy += 1/d
    df = (1/d2)*(1/d)*(y - x)
    output.forces[i] += df
    output.forces[j] -= df
    return output
end
```

To finally define the system and compute the properties:

```julia
argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

system = ParticleSystem(
    xpositions = argon_coordinates,
    unitcell = [21.0, 21.0, 21.0], 
    cutoff = 8.0, 
    output = EnergyAndForces(0.0, similar(argon_coordinates)),
    output_name = :energy_and_forces
);

map_pairwise((x,y,i,j,d2,output) -> energy_and_forces!(x,y,i,j,d2,output), system);
```

The output can be seen with the aliases of the `system.output` variable:
```julia-repl
julia> system.energy_and_forces.energy
207.37593043370862

julia> system.energy_and_forces.forces
100-element Vector{SVector{3, Float64}}:
 [0.02649383330735732, 0.18454277989323772, -0.012253902366284958]
 [0.07782602581235692, 0.27910822337402613, 0.21926615329195248]
 ⋮
 [0.11307234751448932, 0.006353545239676281, -0.05955687310348303]
 [-0.031012009183076745, 0.03543655648545698, 0.03184912163097636]
```

## Updating coordinates, unit cell, and cutoff

If the `map_pairwise!` function will compute energy and/or forces in a iterative procedure (a simulation, for instance), we need to update the coordinates, and perhaps the unit cell and the cutoff.

- [Updating coordinates](@ref)
- [Updating the unit cell](@ref)
- [Updating the cutoff](@ref)

### Updating coordinates

The coordinates can be updated (mutated, or the array of coordinates can change in size by pushing or deleting particles), simply by directly accessing the `xpositions` field of the system. The `xpositions` array is a `Vector` of `SVector` (from `StaticArrays`), with coordinates copied from the input array provided. Thus, the coordinates in the `ParticleSystem` structure must be updated independently of updates in the original array of coordinates. 

Let us exemplify the interface with the computation of forces:

```julia-repl
julia> using CellListMap, StaticArrays

julia> positions = rand(SVector{3,Float64}, 1000);

julia> system = ParticleSystem(
           xpositions = positions,
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = similar(positions),
           output_name = :forces
       );

julia> system.xpositions[1]
3-element SVector{3, Float64} with indices SOneTo(3):
 0.6391290709055079
 0.43679325975360894
 0.8231829019768698

julia> system.xpositions[1] = zeros(SVector{3,Float64})
3-element SVector{3, Float64} with indices SOneTo(3):
 0.0
 0.0
 0.0

julia> push!(system.xpositions, SVector(0.5, 0.5, 0.5))
1001-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.5491373098208292, 0.23899915605319244, 0.49058287555218516]
 ⋮
 [0.4700394061063937, 0.5440026379397457, 0.7411235688716618]
 [0.5, 0.5, 0.5]
```

!!! warning
    The `output` variable may have to be resized accordingly, depending on
    the calculation being performed. Use the `resize_output!` function 
    (do **not** use `Base.resize!` on your output array directly).

If the `output` array has to be resized, that has to be done
with the  `resize_output!` function, which will keep the consistency
of the auxiliary multi-threading buffers. This is, for instance, the case 
in the example of computation of forces, as the `forces` array must be of the
same length as the array of positions:

```julia-repl
julia> resize_output!(system, length(system.xpositions));

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
julia> using CellListMap, StaticArrays

julia> system = ParticleSystem(;
           positions=rand(SVector{3,Float64}, 1000),
           unitcell=[1.0, 1.0, 1.0],
           cutoff=0.1,
           output = 0.0,
        );

julia> update_unitcell!(system, SVector(1.2, 1.2, 1.2))
ParticleSystem1 of dimension 3, composed of:
    Box{OrthorhombicCell, 3}
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

!!! note
    - The unit cell can be set initially using a vector or a unit cell matrix. If a vector is provided
      the system is considered Orthorhombic, if a matrix is provided, a Triclinic system is built. 
      Unit cells updates must preserve the system type. 
    - The unit cell of non-periodic systems (initialized with `nothing`) cannot be updated manually.

    - It is recommended (but not mandatory) to use static arrays (or Tuples) to update the unitcell, 
      as in this case the update will be non-allocating. 


### Updating the cutoff

The cutoff can also be updated, using the `update_cutoff!` function:

```julia-repl
julia> using CellListMap, StaticArrays

julia> system = ParticleSystem(;
           positions=rand(SVector{3,Float64}, 1000),
           unitcell=[1.0, 1.0, 1.0],
           cutoff=0.1,
           output = 0.0,
        );

julia> update_cutoff!(system, 0.2)
ParticleSystem1{output} of dimension 3, composed of:
    Box{OrthorhombicCell, 3}
      unit cell matrix = [ 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0 ]
      cutoff = 0.2
      number of computing cells on each dimension = [8, 8, 8]
      computing cell sizes = [0.2, 0.2, 0.2] (lcell: 1)
      Total number of cells = 512
    CellList{3, Float64}
      1000 real particles.
      636 cells with real particles.
      1738 particles in computing box, including images.
    Parallelization auxiliary data set for 8 batch(es).
    Type of output variable (output): Float64

```

## Computations for two sets of particles

If the computation involves two sets of particle, a similar interface is available. 
The only difference is that the coordinates of the two sets must be provided to
the `ParticleSystem` constructor as the `xpositions` and `ypositions` arrays.

We will illustrate this interface by computing the minimum distance between two
sets of particles, which allows us to showcase further the definition of custom
type interfaces:

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
julia> import CellListMap: copy_output, reset_output!, reducer!

julia> copy_output(md::MinimumDistance) = md
copy_output (generic function with 5 methods)

julia> reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
reset_output! (generic function with 5 methods)

julia> reducer!(md1::MinimumDistance, md2::MinimumDistance) = md1.d < md2.d ? md1 : md2
reducer! (generic function with 2 methods)
```
Note that since `MinimumDistance` is immutable, copying it is the same as returning the value. 
Also, resetting the minimum distance consists of setting its `d` field to `+Inf`. And, finally,
reducing the threaded distances consists of keeping the pair with the shortest distance. 

Next, we build the system

```julia-repl
julia> xpositions = rand(SVector{3,Float64},1000);

julia> ypositions = rand(SVector{3,Float64},1000);

julia> system = ParticleSystem(
           xpositions = xpositions,
           ypositions = ypositions, 
           unitcell=[1.0,1.0,1.0], 
           cutoff = 0.1, 
           output = MinimumDistance(0,0,+Inf),
           output_name = :minimum_distance,
        )
```

And finally we can obtain the minimum distance between the sets: 

```julia-repl
julia> map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), system)
MinimumDistance(276, 617, 0.006009804808785543)
```

## Additional options

- [Turn parallelization on and off](@ref)
- [Displaying a progress bar](@ref)
- [Fine control of the parallelization](@ref)
- [Avoid cell list updating](@ref)
- [Control CellList cell size](@ref)
- [Coordinates as matrices](@ref)

### Turn parallelization on and off

The use of parallel computations can be tunned on and of by the `system.parallel` boolean flag.
For example, using 6 cores (12 threads) for the calculation of the minimum-distance example: 

```julia-repl
julia> f(system) = map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), system)
f (generic function with 1 method)

julia> Threads.nthreads()
8

julia> system.parallel = true
true

julia> @btime f($system)
  268.265 μs (144 allocations: 16.91 KiB)
MinimumDistance(783, 497, 0.007213710914619913)

julia> system.parallel = false
false

julia> @btime f($system)
  720.304 μs (0 allocations: 0 bytes)
MinimumDistance(783, 497, 0.007213710914619913)
```

### Displaying a progress bar

Displaying a progress bar: for very long runs, the user might want to see the progress
of the computation. Use the `show_progress` keyword parameter of the `map_pairwise!` 
function for that.

For example, we execute the computation above, but with much more
particles:

```julia-repl
julia> xpositions = rand(SVector{3,Float64},10^6);

julia> ypositions = rand(SVector{3,Float64},10^6);

julia> system = ParticleSystem(
                  xpositions = xpositions,
                  ypositions = ypositions, 
                  unitcell=[1.0,1.0,1.0], 
                  cutoff = 0.1, 
                  output = MinimumDistance(0,0,+Inf),
                  output_name = :minimum_distance,
               );

julia> map_pairwise(
           (x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), system; 
           show_progress = true
       )
Progress:  24%|██████████▏                               |  ETA: 0:00:29
```

By activating the `show_progress` flag, a nice progress bar is shown. 

### Fine control of the parallelization

The number of batches launched in parallel runs can be tunned by the 
`nbatches` keyword parameter of the `ParticleSystem` constructor. 
By default, the number of batches is defined as heuristic function 
dependent on the number of particles, and possibly returns optimal
values in most cases. For a detailed discussion about this parameter, 
see [Number of batches](@ref Number-of-batches).

For example, to set the number of batches for cell list calculation
to 4 and the number of batches for mapping to 8, we can do:

```julia-repl
julia> system = ParticleSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0,
           output_name = :energy,
           nbatches=(4,8), # use this keyword 
       );
```

Most times it is expected that the default parameters are optimal. But particularly for
inhomogeneous systems increasing the number of batches of the mapping phase (second
parameter of the tuple) may improve the performance by reducing the idle time of
threads.

When the number of batches is left at the default (i.e., `nbatches=(0,0)` or omitted),
it is automatically recomputed whenever `UpdateParticleSystem!` detects that the number
of particles has changed. This allows adding or removing particles from the system
without having to manually adjust the parallelization parameters:

```julia-repl
julia> system = ParticleSystem(
           xpositions = rand(SVector{3,Float64}, 1000),
           unitcell = [1,1,1],
           cutoff = 0.1,
           output = 0.0,
           output_name = :energy,
       );

julia> nbatches(system) # default batches for 1000 particles
(2, 4)

julia> system.xpositions = rand(SVector{3,Float64}, 100000); # resize

julia> UpdateParticleSystem!(system); # nbatches recomputed automatically

julia> nbatches(system) # updated for 100000 particles
(8, 32)
```

If the number of batches is explicitly set to non-zero values, they will be kept fixed
and will not change when the number of particles changes.

### Avoid cell list updating

To compute different properties without recomputing cell lists, use `update_lists=false` in 
the call of `map_pairwise` methods, for example,
```julia
using CellListMap, StaticArrays
system = ParticleSystem(xpositions=rand(SVector{3,Float64},1000), output=0.0, cutoff=0.1, unitcell=[1,1,1])
# First call, will compute the cell lists
map_pairwise((x,y,i,j,d2,u) -> u += d2, system)
# Second run: do not update the cell lists but compute a different property
map_pairwise((x,y,i,j,d2,u) -> u += sqrt(d2), system; update_lists = false)
```
in which case we are computing the sum of distances from the same cell lists used to compute the energy in the previous example
(requires version 0.8.9). 

!!! warning
    This option will skip the updating of the cell lists, thus be careful to **not** use this
    option if the coordinates, cutoff, unitcell, or any other property of the system changed. 

### Control CellList cell size

The cell sizes of the construction of the cell lists can be controlled with the keyword `lcell`
of the `ParticleSystem` constructor. For example:
```julia-repl
julia> system = ParticleSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0,
           output_name = :energy,
           lcell=2,
       );
```
Most times using `lcell=1` (default) or `lcell=2` will provide the optimal performance. For very
dense systems, or systems for which the number of particles within the cutoff is very large,
larger values of `lcell` may improve the performance. To be tested by the user.

!!! note
    The number of cells in which the particles will be classified is, for each dimension `lcell*length/cutoff`. 
    Thus if the `length` of the box is too large relative to the `cutoff`, many cells will be created, and this
    imposes a perhaps large memory requirement. Usually, it is a good practice to limit the number of cells to
    be not greater than the number of particles, and for that the cutoff may have to be increased, if there is
    a memory bottleneck. A reasonable choice is to use `cutoff = max(real_cutoff, length/n^(1/D))` where `n` is the 
    number of particles and `D` is the dimension (2 or 3). With that the number of cells will be close to `n` in the worst case.  

### Coordinates as matrices

Coordinates can also be provided as matrices of size `(D,N)` where `D` is the dimension (2 or 3) and `N` is the number of particles. For example:

```jldoctest; filter = r"\d+" => "" 
julia> using CellListMap

julia> system = ParticleSystem(
           xpositions=rand(2,100),
           ypositions=rand(2,200),
           cutoff=0.1,
           unitcell=[1,1],
           output=0.0,
       )
ParticleSystem2{output} of dimension 2, composed of:
    Box{OrthorhombicCell, 2}
      unit cell matrix = [ 1.0 0.0; 0.0 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [13, 13]
      computing cell sizes = [0.1, 0.1] (lcell: 1)
      Total number of cells = 169
    CellListMap.CellListPair{2, Float64}
       63 cells with real particles of the smallest set.
       85 cells with real particles of the largest set.
    Parallelization auxiliary data set for 1 batch(es).
    Type of output variable (output): Float64
```
!!! warning
    This interface less flexible than when the coordinates are input as vectors of vectors, because
    *the number of particles* cannot be changed, because matrices cannot be resized. Otherwise, matrices can
    be used as input.

## Complete example codes

- [Simple energy computation](@ref)
- [Force computation](@ref)
- [Energy and forces](@ref)
- [Two sets of particles](@ref)
- [Particle simulation](@ref)

### Simple energy computation

In this example, a simple potential energy defined as the sum of the 
inverse of the distance between the particles is computed.

```julia
using CellListMap
using StaticArrays
system = ParticleSystem(
    xpositions = rand(SVector{3,Float64},1000), 
    unitcell=[1.0,1.0,1.0], 
    cutoff = 0.1, 
    output = 0.0,
    output_name = :energy
)
map_pairwise!((x,y,i,j,d2,energy) -> energy += 1 / sqrt(d2), system)
```

### Force computation

Here we compute the force vector associated to the potential energy
function of the previous example.

```julia
using CellListMap
using StaticArrays
positions = rand(SVector{3,Float64},1000) 
system = ParticleSystem(
    xpositions = positions, 
    unitcell=[1.0,1.0,1.0], 
    cutoff = 0.1, 
    output = similar(positions),
    output_name = :forces
)
function update_forces!(x,y,i,j,d2,forces)
    d = sqrt(d2)
    df = (1/d2)*(1/d)*(y - x)
    forces[i] += df
    forces[j] -= df
    return forces
end
map_pairwise!((x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces), system)
```

### Energy and forces

In this example, the potential energy and the forces are computed in a single
run, and a custom data structure is defined to store both values.

```julia
using CellListMap
using StaticArrays
# Define custom type
mutable struct EnergyAndForces
    energy::Float64
    forces::Vector{SVector{3,Float64}}
end
# Custom copy, reset and reducer functions
import CellListMap: copy_output, reset_output!, reducer
copy_output(x::EnergyAndForces) = EnergyAndForces(copy(x.energy), copy(x.forces))
function reset_output!(output::EnergyAndForces)
    output.energy = 0.0
    for i in eachindex(output.forces)
        output.forces[i] = SVector(0.0, 0.0, 0.0)
    end
    return output
end
function reducer(x::EnergyAndForces, y::EnergyAndForces)
    e_tot = x.energy + y.energy
    x.forces .+= y.forces
    return EnergyAndForces(e_tot, x.forces)
end
# Function that updates energy and forces for each pair
function energy_and_forces!(x,y,i,j,d2,output::EnergyAndForces)
    d = sqrt(d2)
    output.energy += 1/d
    df = (1/d2)*(1/d)*(y - x)
    output.forces[i] += df
    output.forces[j] -= df
    return output
end
# Initialize system
positions = rand(SVector{3,Float64},1000);
system = ParticleSystem(
    xpositions = positions,
    unitcell=[1.0,1.0,1.0], 
    cutoff = 0.1, 
    output = EnergyAndForces(0.0, similar(positions)),
    output_name = :energy_and_forces
)
# Compute energy and forces
map_pairwise((x,y,i,j,d2,output) -> energy_and_forces!(x,y,i,j,d2,output), system)
```

### Two sets of particles

In this example we illustrate the interface for the computation of properties
of two sets of particles, by computing the minimum distance between the two sets.

```julia
using CellListMap
using StaticArrays
# Custom structure to store the minimum distance pair
struct MinimumDistance
    i::Int
    j::Int
    d::Float64
end
# Function that updates the minimum distance found
function minimum_distance(i, j, d2, md)
    d = sqrt(d2)
    if d < md.d
        md = MinimumDistance(i, j, d)
    end
    return md
end
# Define appropriate methods for copy, reset and reduce 
import CellListMap: copy_output, reset_output!, reducer!
copy_output(md::MinimumDistance) = md
reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
reducer!(md1::MinimumDistance, md2::MinimumDistance) = md1.d < md2.d ? md1 : md2
# Build system 
xpositions = rand(SVector{3,Float64},100);
ypositions = rand(SVector{3,Float64},1000);
system = ParticleSystem(
       xpositions = xpositions,
       ypositions = ypositions, 
       unitcell=[1.0,1.0,1.0], 
       cutoff = 0.1, 
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
# Function following the required interface of the mapped function
get_md(_, _, i, j, d2, md) = minimum_distance(i, j, d2, md)
# Compute the minimum distance
map_pairwise(get_md, system)
```

In the above example, the function is used such that cell lists are constructed for both
sets. There are situations where this is not optimal, in particular:

1. When one of the sets if very small. In this case, constructing a cell list for the largest
   set becomes the bottleneck. Therefore, it is better to construct a cell list for the smallest
   set and loop over the particles of the largest set.
2. When one of the set is fixed and the second set is variable. In this case, it is better to
   construct the cell list for the fixed set only and loop over the variables of the variable set.

For dealing with these possibilities, an additional two-set interface is available, where one maps
the computation over an array of particles relative to a previously computed cell list. Complementing
the example above, we could compute the same minimum distance using:

```julia
# Construct the cell list system only for one of the sets: ypositions
ysystem = ParticleSystem(
       positions = ypositions,
       unitcell=[1.0,1.0,1.0], 
       cutoff = 0.1, 
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
# obtain the minimum distance between xpositions and the cell list in system
# Note the additional `xpositions` parameter in the call to map_pairwise.
map_pairwise(get_md, xpositions, ysystem)
```

Additionally, if the `xpositions` are updated, we can obtain compute the function relative to `ysystem` without
having to update the cell lists: 

```julia-repl
julia> xpositions = rand(SVector{3,Float64},100);

julia> map_pairwise(get_md, xpositions, ysystem)
MinimumDistance(67, 580, 0.008423693268450603)
```

while with the two-set cell list system one would need to update the cell lists for this new computation.

!!! compat
    The single-set cross-interaction was introduced in v0.10.0. It uses the method previously implemented
    for all cross-interactions. 

### Benchmarking of cross-interaction alternatives

With the following functions we will benchmark the performance of the two alternatives for computing
cross-set interactions, **including** the time required to build the cell lists (the initialization
of the `ParticleSystem`) objects:

```julia
# First alternative: compute cell lists for the two sets
function two_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       xpositions = xpositions,
       ypositions = ypositions, 
       unitcell=[1.0,1.0,1.0], 
       cutoff = 0.1, 
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
    return map_pairwise(get_md, system) 
end
# Second alternative: compute cell lists for one set
function one_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       positions = ypositions,
       unitcell=[1.0,1.0,1.0], 
       cutoff = 0.1, 
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
    return map_pairwise(get_md, xpositions, system) 
end
```

If one of the sets is small, the one-set alternative is clearly faster, if we 
construct the cell lists for the smaller set:

```julia-repl
julia> using BenchmarkTools 

julia> xpositions = rand(SVector{3,Float64}, 10^6);

julia> ypositions = rand(SVector{3,Float64}, 100);

julia> @btime one_set_celllist($xpositions, $ypositions) samples=1 evals=1
  25.165 ms (1531 allocations: 575.72 KiB)
MinimumDistance(65937, 63, 0.00044803040276614203)

julia> @btime two_set_celllist($xpositions, $ypositions) samples=1 evals=1
  207.129 ms (154794 allocations: 478.00 MiB)
MinimumDistance(65937, 63, 0.00044803040276614203)
```

For much larger system, though, the computation of the cell lists become less relevant and the first alternative 
might be the most favorable, even including the cell lists updates:

```julia-repl
julia> @btime one_set_celllist($xpositions, $ypositions) samples=1 evals=1
  12.196 s (153327 allocations: 478.02 MiB)
MinimumDistance(627930, 889247, 7.59096139675071e-5)

julia> @btime two_set_celllist($xpositions, $ypositions) samples=1 evals=1
  2.887 s (306416 allocations: 952.00 MiB)
MinimumDistance(627930, 889247, 7.59096139675071e-5)
```

This performance advantage of the two-set cell lists arises because more interactions can be skipped
by [precomputing properties of the cells involved](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20563). 
On the other side, when the lists are available 
for only one set, the loop over all the particles of the second set is mandatory. Since this loop
is fast, it is favorable over the construction of the cell lists for smaller sets.  

### Particle simulation

In this example, a complete particle simulation is illustrated, with a simple potential.  This example can illustrate how particle positions and forces can be updated. Run this
simulation with:

```julia-repl
julia> system = init_system(N=200); # number of particles

julia> trajectory = simulate(system);

julia> animate(trajectory)
```

One important characteristic of this example is that the `system` is built outside the function that performs the simulation. This is done because the construction of the system is type-unstable (it is dimension, geometry and output-type dependent). Adding a function barrier avoids type-instabilities to propagate to the simulation causing possible performance problems. 

```julia
using StaticArrays
using CellListMap
import CellListMap.wrap_relative_to
# Function that updates the forces, for potential of the form:
# if d < cutoff k*(d^2-cutoff^2)^2 else 0.0 with k = 10^6
function update_forces!(x, y, i, j, d2, forces, cutoff)
    r = y - x
    dudr = 10^6 * 4 * r * (d2 - cutoff^2)
    forces[i] += dudr
    forces[j] -= dudr
    return forces
end
# Function that initializes the system: it is preferable to initialize
# the system outside the function that performs the simulation, because
# the system (data)type is defined on initialization. Initializing it outside
# the simulation function avoids possible type-instabilities. 
function init_system(;N::Int=200)
    Vec2D = SVector{2,Float64}
    positions = rand(Vec2D, N)
    unitcell = [1.0, 1.0]
    cutoff = 0.1
    system = ParticleSystem(
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=similar(positions),
        output_name=:forces,
    )
    return system
end
function simulate(system=init_system(); nsteps::Int=100, isave=1)
    # initial velocities
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    dt = 1e-3
    trajectory = typeof(system.positions)[]
    for step in 1:nsteps
        # compute forces at this step
        map_pairwise!(
            (x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces,system.cutoff),
            system
        )
        # Update positions and velocities
        for i in eachindex(system.positions, system.forces)
            f = system.forces[i]
            x = system.positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)
            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x
            velocities[i] = v
        end
        # Save step for printing
        if step % isave == 0
            push!(trajectory, copy(system.positions))
        end
    end
    return trajectory
end

using Plots
function animate(trajectory)
    anim = @animate for step in trajectory
        scatter(
            Tuple.(step),
            label=nothing,
            lims=(-0.5, 0.5),
            aspect_ratio=1,
            framestyle=:box,
        )
    end
    gif(anim, "simulation.gif", fps=10)
end
```

## Docstrings

```@meta
CollapsedDocStrings = true
```

```@autodocs
Modules = [CellListMap]
Pages = ["ParticleSystem.jl"]
Order = [:function, :type]
```

