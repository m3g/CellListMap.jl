# CellListMap.jl

[CellListMap.jl](https://github.com/m3g/CellListMap.jl) implements an efficient cell list scheme for the computation of interactions, neighbor lists, or any other property dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. 

The package provides an interface to compute a generic function for each pair of particles closer 
than a cutoff, using general periodic boundary conditions. Parallel and serial implementations can be used. 

## Overview

`CellListMap` is a package that implements a fast scheme for computing properties of systems of particles in 2 or 3 dimensions, within a cutoff. In brief, it is designed to replace double loops running over the pairs of particles of a system. Naively, a loop over all pair of particles is written as:
```julia
for i in 1:N
    for j in i+1:N
        # compute distance, possibly considering periodic boundary conditions
        d = distance(particle[i],particle[j]) 
        if d <= cutoff 
            # compute property dependent on d
        end
    end
end
```
where `N` is the number of particles. 

Alternatively, if the interaction is between two disjoint sets of particles, the naive loop is
```julia
for i in 1:N 
    for j in 1:M
        # compute distance, possibly considering periodic boundary conditions
        d = distance(first_set[i], second_set[j])
        if d <= cutoff
            # compute property dependent on d
        end
    end
end
```
where `N` and `M` are the numbers of particles of each set. If the cutoff is significantly smaller than the dimension of the system,
these loops are very expensive, and it is possible to avoid looping over particles that are farther from each other than the cutoff.

CellListMap implements an efficient and parallel cell-list method, with optimizations, to substitute such double-loops while taking into account
periodic boundary conditions. Cell lists are an alternative to distance trees and are particularly effective when the distribution
of the particles is roughly homogeneous. For highly heterogeneous systems distance trees like those implemented in
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) might be more performant.

## Quick start: computing neighbor lists

The simplest use case of `CellListMap` is to compute the list of neighboring particles within a cutoff distance.
This can be implemented as a particular usage case of `CellListMap`, but given its generality it is provided as an independent interface: the `neighborlist` function:

An example of the application of the `neighborlist` function follows:

```julia-repl
julia> using CellListMap

julia> x = rand(3, 10_000); # 10,000 particles in 3D

julia> neighborlist(x, 0.05) # cutoff of 0.05
209452-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1055, 0.04770602450750023)
 (1, 1261, 0.048087995599677004)
 ⋮
 (9989, 8269, 0.0442616226498007)
 (9989, 6452, 0.042189574206507035)
```

Each element of the list is a tuple `(i, j, d)` containing the indices of the particles and their distance.

`CellListMap` supports general periodic boundary conditions, by providing the unitcell as vector (for orthorhombic systems) or a unitcell matrix (for general triclinic systems):

```julia-repl
julia> neighborlist(x, 0.05; unitcell=[1,1,1]) # periodic box of side 1
305006-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1055, 0.04770602450750023)
 ⋮
```

See the [Neighbor lists](@ref) section for more details, including in-place computations for iterative workflows.

!!! note
    The downside of explicitly computing neighbor lists is that it implies storing the list of neighbors. This can be memory consuming and slow in situations where a property can be computed by iterating over the neighbors without materializing the neighbor lists. `CellListMap` is designed to provide this alternative.

## General pairwise computations

For more general pairwise computations (energies, forces, etc.) without the materialization of the neighbor lists, the `ParticleSystem` interface provides a flexible way to define custom functions that are applied to all pairs of particles within the cutoff, and reducing an output value. See the [ParticleSystem interface](@ref) section for details.

A concise example is the computation of the sum of the inverse of the distance between particles:

```julia-repl
julia> using CellListMap

julia> x = rand(3, 10_000); # 10,000 particles in 3D

julia> sys = ParticleSystem(positions=x, cutoff=0.05, unitcell=[1,1,1], output=0.0)

julia> pairwise!((pair, u) -> u += 1/pair.d, sys)
792925.6234732079
```

Note that in the above example the `pairwise` method is actually performing a `pairwise` operation, where the output value `u` is summed up over all neighboring pairs of particles. 

## Installation

This is a [Julia](http://julialang.org) package. Install `Julia` first following the instructions in the [download page](https://julialang.org/downloads/).

Once Julia is installed, install the `CellListMap` package from the Julia REPL with:

```julia-repl
julia> import Pkg

julia> Pkg.add("CellListMap")
```

## Help!

Please ask for help if having any difficulty using the package. Reach us by:

- Asking a question on the [Julia Discourse forum](https://discourse.julialang.org/). Please
  mark `@lmiq` on your post, otherwise we may miss it! This may be very effective to get help from 
  many Julia users on questions that might not be directly related this package.
- [Opening an issue](https://github.com/m3g/CellListMap.jl/issues/new/choose) if you think you found a problem in the package.
  Even documentation problems can be reported.
- Joining us at Zulip-chat in the [m3g stream](https://julialang.zulipchat.com/#narrow/stream/435348-m3g) of the Julia Zulip forum.
- Sending an e-mail to: [lmartine@unicamp.br](mailto:lmartine@unicamp.br?subject="CellListMap.jl help").

