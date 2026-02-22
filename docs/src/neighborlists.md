# Neighbor lists

Neighbor lists can be computed, returning all pairs of particles that are found within the cutoff,
and the corresponding distances.  

- [Non-periodic systems](@ref)
- [Periodic systems](@ref)
- [In-place computation of neighbor lists](@ref)
- [Options](@ref)

!!! note
    - When computing neighbor lists with cell-lists, it is possible for pairs of particles that are at a distance equal to the cutoff to either be included or excluded due to numerical rounding. As a result, these neighbor lists should only be utilized for calculating properties that vanish at the cutoff distance.
    - Positions and unit cells can be 2 or 3-dimensional.
    - The order of the pairs in the output list is not guaranteed and may change, in particular,
    for parallel executions.

    
## Non-periodic systems

Without periodic boundary conditions, just provide the coordinates and the cutoff:

```@example nb
using CellListMap
x = [ rand(2) for _ in 1:10_000 ]
neighborlist(positions=x, cutoff=0.01)
```

If the neighbor lists between two sets of points are required, use the following notation, where `xpositions` and `ypositions` define the input coordinates. In this case using coordinates as arrays of static arrays:

```@example nb
using StaticArrays
x = rand(SVector{3,Float64}, 10^4)
y = rand(SVector{3,Float64}, 10^3)
nblist = neighborlist(xpositions=x, ypositions=y, cutoff=0.01)
```

The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

## Periodic systems

If periodic boundary conditions are used, the `unitcell` can be provided explicitly as keyword parameters:

```@example nb
x = [ rand(2) for _ in 1:10_000 ] 
neighborlist(positions=x, cutoff=0.01, unitcell=[1,1])
```

In the example above, an `Orthorhombic` cell was assumed automatically from the fact that a vector of sides was provided. For general periodic boundary conditions, a unit cell matrix must be given:

```@example nb
neighborlist(positions=x, cutoff=0.01, unitcell=[1.0 0.5; 0.5 1.0])
```

## In-place computation of neighbor lists

If neighbor lists are computed within a interactive scenario, it is interesting preallocate all the necessary data and just update the lists at every iteration. This can be achieved by constructing the `InPlaceNeighborList` object in advance. The performance gain of performing the operations in place might vary and may not be important for single runs, as the allocations do not dominate the computing time. 

We will first illustrate the interface for a non-parallel run:
```@example nb
x = rand(SVector{3,Float64}, 10^4)
system = InPlaceNeighborList(
    positions=x, 
    cutoff=0.1, 
    unitcell=[1,1,1], 
    parallel=false
)
```

Note that the buffer length has size `0`. The first time the neighbor lists are computed, the list will be allocated. We will use the `neighborlist!` (with the bang) function, because it will effectively mutate the `system`, by allocating all necessary data:

```@example nb
@time list = neighborlist!(system)
```

Now, if we modify the coordinates, we can update the system and recompute the neighbor lists:

```@example nb
x_new = rand(SVector{3,Float64}, 10^4)
update!(system; positions=x_new)
@time list = neighborlist!(system)
```

where minimal allocations occur (in fact, if the length of the neighborlist does not increase, there will be no allocations).

## Updating the `InPlaceNeighborList` object

The `update!` function provides an interface to modify the coordinates, unitcell, and cutoff, and perform a new neighbor list search. In the example above we have updated the positions, `positions`, of the system. The full keyword list for updating the system are:

| Keyword | Type | Default | Description |
|:--------|:-----|:-------:|:------------|
| `positions` / `xpositions` | `AbstractVector` or `AbstractMatrix` | `nothing` | New coordinates for the (first) set of particles. The internal buffer is resized automatically if the number of particles changes. `positions` is an alias for `xpositions`. |
| `ypositions` | `AbstractVector`, `AbstractMatrix`, or `Nothing` | `nothing` | New coordinates for the second set of particles (two-set systems only). |
| `cutoff` | `Number` | `nothing` | New cutoff distance. |
| `unitcell` | vector or matrix | `nothing` | New unit cell. Cannot be changed for non-periodic systems. |
| `parallel` | `Bool` | `nothing` | Enable or disable multi-threading. |

!!! note
    - Allocations can occur if the cutoff, unit cell, or number of particles change such
      that greater buffers are required. The number of allocations tend to diminish as
      the buffers become large enough to accommodate the possible variations of the computation.

For parallel runs, the allocations are minimal, but some small auxiliary data is required for the launching of multiple threads. We illustrate here the convergence of the allocations to the minimum required for multi-threaded calculations:

```@example nb
using Chairmarks
update!(system; parallel=true)
@b neighborlist!(system)
```

## Options

### Input parameters: `neighborlist` and `InPlaceNeighborList`

| Keyword | Type | Default | Description |
|:--------|:-----|:-------:|:------------|
| `positions` / `xpositions` | `AbstractVector` or `AbstractMatrix` | — (required) | Coordinates of the (first) set of particles. `positions` is an alias for `xpositions` for single-set computations. |
| `ypositions` | `AbstractVector`, `AbstractMatrix`, or `Nothing` | `nothing` | Coordinates of the second set of particles. If provided, cross-pair neighbor lists between the two sets are computed. |
| `cutoff` | `Number` | — (required) | Cutoff distance. Pairs within this distance are included in the list. |
| `unitcell` | vector or matrix | `nothing` | Unit cell sides (orthorhombic, as a vector) or unit cell matrix (triclinic). If `nothing`, the system is treated as non-periodic. |
| `parallel` | `Bool` | `true` | Enable multi-threading for the computation. |
| `show_progress` | `Bool` | `false` | Display a progress bar during computation. Only available for `InPlaceNeighborList` / `neighborlist!`. |
| `nbatches` | `Tuple{Int,Int}` | `(0,0)` | Number of batches used in parallelization (see [Number of batches](@ref Number-of-batches)). |

### Output

The neighbor list is a `Vector{Tuple{Int,Int,T}}`, where each element `(i, j, d)` contains
the indices of the two particles and their distance `d`. The floating-point type `T` is
promoted from the `cutoff` and `unitcell` types.

`neighborlist` returns this vector directly. `neighborlist!` returns the same vector, stored
in `system.nb.list`, and updates it in place on repeated calls.

