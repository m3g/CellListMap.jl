# Neighbor lists

Neighbor lists can be computed, returning all pairs of particles that are found within the cutoff,
and the corresponding distances.  

- [Non-periodic systems](@ref)
- [Periodic systems](@ref)
- [In-place computation of neighbor lists](@ref)
- [Options](@ref)

!!! note
    When computing neighbor lists with cell-lists, it is possible for pairs of particles that are at a distance equal to the cutoff to either be included or excluded due to numerical rounding. As a result, these neighbor lists should only be utilized for calculating properties that vanish at the cutoff distance.
    
## Non-periodic systems

Without periodic boundary conditions, just provide the coordinates and the cutoff:

```julia-repl
julia> using CellListMap

julia> x = [ rand(2) for _ in 1:10_000 ];

julia> neighborlist(x,0.05)
376457-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 363, 0.04855594810064624)
 (1, 513, 0.03356381123125866)
 (1, 1209, 0.005159666709130686)
 ⋮
 (6575, 7378, 0.03791567990447959)
 (7378, 3450, 0.01748757015908321)
```

If the neighbor lists between two sets of points are required, use the following notation, 
in this case using coordinates as arrays of static arrays:
```julia-repl
julia> using StaticArrays

julia> x = rand(SVector{3,Float64},10^4);

julia> y = rand(SVector{3,Float64},10^3);

julia> list = neighborlist(x,y,0.1)
37309-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 971, 0.09867846773727411)
 (1, 567, 0.06630101425431004)
 (1, 3, 0.04103170149300593)
 ⋮
 (10000, 156, 0.08549899843141298)
 (10000, 444, 0.0737386384422871)
```

where, similarly, the third parameter is the cutoff.
The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

## Periodic systems

If periodic boundary conditions are used, the `unitcell` can be provided explicitly as keyword parameters:

```julia-repl
julia> x = [ rand(2) for _ in 1:10_000 ]; 

julia> neighborlist(x, 0.05; unitcell=[1,1])
392100-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 5, 0.03445098850037766)
 (1, 393, 0.039448810592487206)
 (1, 1632, 0.02276457565643465)
 ⋮
 (9501, 9781, 0.03351665194098955)
 (9501, 5429, 0.04199258248973222)
```

In the example above, an `Orthorhombic` cell was assumed, and thus a vector of sides was provided. For general
periodic boundary conditions, a unit cell matrix can be provided, for example:

```julia-repl
julia> neighborlist(x, 0.05; unitcell=[1.0 0.5; 0.5 1.0])
580693-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 457, 0.03935441952786555)
 (1, 1467, 0.033407692174569875)
 (1, 1767, 0.04490555313598093)
 ⋮
 (3652, 8475, 0.04721628783510375)
 (6260, 8475, 0.04946130971686825)
```

!!! note
    Positions and unit cells can be 2 or 3-dimensional.

## In-place computation of neighbor lists

If neighbor lists are computed within a interactive scenario, it is interesting preallocate all the necessary
data and just update the lists at every iteration. This can be achieved by constructing the `InPlaceNeighborList` 
object in advance. The performance gain of performing the operations in place might vary and may not be 
important for single runs, as the allocations do not dominate the computing time. 

We will first illustrate the interface for a non-parallel run:
```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 10^4);

julia> system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1], parallel=false)
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{OrthorhombicCell, 3, Float64, Float64, 9}
Current list buffer size: 0
```

Note that the buffer size has size `0`. The first time the neighbor lists are computed, the list will
be allocated. We will use the `neighborlist!` (with the bang) function, because it will effectively 
mutate the `system`, by allocating all necessary data:

```julia-repl
julia> @time list = neighborlist!(system)
  0.017765 seconds (12 allocations: 7.445 MiB)
209190-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1375, 0.09425551992016712)
 (1, 3076, 0.045320021406080775)
 (1, 3666, 0.07780146666634076)
 ⋮
 (9962, 6983, 0.07355578793348823)
 (9962, 7457, 0.07597724209140656)
```

Now, if we modify the coordinates, we can update the system and recompute the neighbor lists:
```julia-repl
julia> x_new = rand(SVector{3,Float64}, 10^4);

julia> @time update!(system, x_new)
  0.003562 seconds
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{OrthorhombicCell, 3, Float64, Float64, 9}
Current list buffer size: 209190

julia> @time list = neighborlist!(system);
  0.012338 seconds
```

!!! note
    - Here we illustrate the behavior of the functions in their second calls, to remove the 
      effects of compilation on the allocation results.
    - The `cutoff` and `unitcell`  can be modified by providing additional keyword parameters
      to the `update!` function (for example `update!(system, x; cutoff=0.1)`).
    - Allocations can occur if the cutoff, unit cell, or number of particles change such
      that greater buffers are required. The number of allocations tend to diminish as 
      the buffers become large enough to accommodate the possible variations of the computation.

For parallel runs, the allocations are minimal, but some small auxiliary data is required for the
launching of multiple threads. We illustrate here the convergence of the allocations to the 
minimum required for multi-threaded calculations:

```julia-repl
julia> system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1], parallel=true);

julia> @time list = neighborlist!(system);
  0.007762 seconds (230 allocations: 18.142 MiB)

julia> x_new = rand(SVector{3,Float64},10^4);

julia> @time update!(system, x_new)
  0.005283 seconds (20.30 k allocations: 6.200 MiB)
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{OrthorhombicCell, 3, Float64, Float64, 9}
Current list buffer size: 209190

julia> @time neighborlist!(system);
  0.008190 seconds (166 allocations: 6.461 MiB)

julia> x_new = rand(SVector{3,Float64},10^4);

julia> @time update!(system, x_new);
  0.002723 seconds (221 allocations: 208.922 KiB)

julia> @time neighborlist!(system);
  0.006227 seconds (165 allocations: 2.863 MiB)

julia> x_new = rand(SVector{3,Float64},10^4);

julia> @time update!(system, x_new);
  0.002396 seconds (275 allocations: 144.078 KiB)

julia> @time neighborlist!(system);
  0.004996 seconds (161 allocations: 15.141 KiB)
```

## Options

Additional optional parameters can be used in a `neighborlist` call:

| Keyword |  Values types | Default | About |
|:-------:|:-------:|:-------:|:------|
| `parallel` | `Bool`  | `true` | turns on and off parallelization |
| `show_progress` | `Bool` | `false` |  turns on and off progress bar | 
| `nbatches` | `Tuple{Int,Int}` | `(0,0)` |  Number of batches used in parallelization (see [here](@ref Number-of-batches)) | 
| `autoswap` | `Bool` | `true` |  (advanced) automatically choose set to construct the cell lists |




