# Updating the system

If the `pairwise!` function will compute energy and/or forces in an iterative procedure (a simulation, for instance), we need to update the coordinates, and perhaps the unit cell and the cutoff.

## Updating coordinates

The coordinates can be updated (mutated, or the array of coordinates can change in size by pushing or deleting particles), simply by directly accessing the `xpositions` (and/or `ypositions`) field of the system. These position arrays are of type [`ParticleSystemPositions`](@ref ParticleSystemPositions), which wraps a `Vector{SVector{N,T}}` and tracks mutations automatically. When positions are modified through the supported interface (see [The `ParticleSystemPositions` type](@ref ParticleSystemPositions)), the cell lists are recomputed automatically on the next call to `pairwise!`. Thus, the coordinates in the `ParticleSystem` structure must be updated independently of updates in the original array of coordinates.

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

    In the case of compound outputs (custom output structures) like that
    of the [Computing both energy and forces](@ref) example, calling `resize_output!` will
    return an error and require the user to define `Base.resize!` for the
    custom type.

If the `output` array has to be resized, that has to be done
with the  `resize_output!` function, which will keep the consistency
of the auxiliary multi-threading buffers. This is, for instance, the case
in the example of computation of forces, as the `forces` array must be of the
same length as the array of positions:

```julia-repl
julia> resize_output!(system, length(system.xpositions));

julia> pairwise!(update_forces!, system)
1001-element Vector{SVector{3, Float64}}:
 [756.2076075886971, -335.1637545330828, 541.8627090466914]
 [-173.02442398784672, -178.782819965489, 4.570607952876692]
 ⋮
 [-722.5400961501635, 182.65287417718935, 380.0394926753039]
 [20.27985502389337, -193.77607810950286, -155.28968519541544]
```

In this case, if the `output` is not resized, a `BoundsError` is
be obtained, because updates of forces at unavailable positions will
be attempted.

## Updating the unit cell

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


## Updating the cutoff

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
ParticleSystem1{default_output_name} of dimension 3, composed of:
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
    Type of output variable (default_output_name): Float64

```
