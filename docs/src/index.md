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

### High level interface for particle systems

Since version `0.8.30`, a simpler, higher level interface was introduced, that will facilitate the use of `CellListMap` without any loss in performance. The new interface is flexible enough for the majority of applications. See the [ParticleSystem interface](@ref) menu for details. 

### Cutoff-delimited neighbor lists

The user might be more comfortable in using the package to compute the list of neighboring particles. A custom interface for this application is provided though the [Neighbor lists](@ref) interface. 

Note that, in general, neighbor lists are used to compute other pairwise dependent properties, and these can be, in principle, computed directly with `CellListMap` without the need to explicitly compute or store the lists of neighbors. 

### Lower level interface

The [Low level interface](@ref) allows the customization and optimization of very demanding calculations (although the ParticleSystem interface does not have any performance limitation and is easier to use).

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

