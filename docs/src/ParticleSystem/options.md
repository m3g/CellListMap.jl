# Additional options

This section covers additional options for fine-tuning the behavior and performance
of the `ParticleSystem` interface.

## Turn parallelization on and off

The use of parallel computations can be tunned on and of by the `system.parallel` boolean flag.
For example, using 6 cores (12 threads) for the calculation of the minimum-distance example:

```julia-repl
julia> f(system) = map_pairwise(minimum_distance, system)
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

## Displaying a progress bar

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

julia> map_pairwise(minimum_distance, system; show_progress = true)
Progress:  24%|██████████▏                               |  ETA: 0:00:29
```

By activating the `show_progress` flag, a nice progress bar is shown.

## Fine control of the parallelization

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

## Avoid cell list updating

To compute different properties without recomputing cell lists, use `update_lists=false` in
the call of `map_pairwise` methods, for example,
```julia
using CellListMap, StaticArrays
system = ParticleSystem(xpositions=rand(SVector{3,Float64},1000), output=0.0, cutoff=0.1, unitcell=[1,1,1])
# First call, will compute the cell lists
map_pairwise!((pair, u) -> u += pair.d2, system)
# Second run: do not update the cell lists but compute a different property
map_pairwise!((pair, u) -> u += pair.d, system; update_lists = false)
```
in which case we are computing the sum of distances from the same cell lists used to compute the energy in the previous example
(requires version 0.8.9).

!!! warning
    This option will skip the updating of the cell lists, thus be careful to **not** use this
    option if the coordinates, cutoff, unitcell, or any other property of the system changed.

## Control CellList cell size

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

## Coordinates as matrices

Coordinates can also be provided as matrices of size `(D,N)` where `D` is the dimension (2 or 3) and `N` is the number of particles. For example:

```jldoctest; filter = r"(\d+|CellListMap\.)" => ""
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
    CellListPair{2, Float64}
       70 cells with real particles of the smallest set.
       85 cells with real particles of the largest set.
    Parallelization auxiliary data set for 4 batch(es).
    Type of output variable (output): Float64
```
!!! warning
    This interface less flexible than when the coordinates are input as vectors of vectors, because
    *the number of particles* cannot be changed, because matrices cannot be resized. Otherwise, matrices can
    be used as input.
