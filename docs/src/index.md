# CellListMap.jl

[CellListMap.jl](https://github.com/m3g/CellListMap.jl) implements an efficient cell list scheme for the computation of interactions, neighbor lists, or any other property dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. 

It maps a generic function to be computed for each pair of particles, using periodic boundary conditions of any type. Parallel and serial implementations can be used. The user defined function will be evaluated only for the pairs closer to each other than the desired cutoff distance.

## Installation

This is a [Julia](http://julialang.org) package. Install `Julia` first following the instructions in the [download page](https://julialang.org/downloads/), or using [juliaup](https://github.com/JuliaLang/juliaup).

Once Julia is installed, install the `CellListMap` package from the Julia REPL with:

```julia-repl
julia> import Pkg

julia> Pkg.add("CellListMap")
```

## Help!

Please ask for help if having any difficulty using the package. Reach us by:

- Asking a question on the [Julia Discourse forum](https://discourse.julialang.org/). Please
  mark `@lmiq` on your post, otherwise we may miss it! This may be very effective to get help from 
  many Julia users on questions not directly related to physical-chemistry.
- [Opening an issue](https://github.com/m3g/ComplexMixtures.jl/issues/new/choose) if you think you found a problem in the package.
  Even documentation problems can be reported.
- Joining us at Zulip-chat in the [m3g stream](https://julialang.zulipchat.com/#narrow/stream/435348-m3g) of the Julia Zulip forum.
- Sending an e-mail to: [lmartine@unicamp.br](mailto:lmartine@unicamp.br?subject="ComplexMixtures.jl help").

## Overview

### High level interface for periodic system

Since version `0.7.22`, a new simpler, higher level interface was introduced, that will facilitate the use of `CellListMap` without any loss in performance. The new interface is flexible enough for the majority of applications. It may become the default interface in the future. See the [PeriodicSystems interface](@ref) menu for details. 

### Cutoff-delimited neighbor lists

The user might be more comfortable in using the package to compute the list of neighboring particles. A custom interface for this application is provided though the [Neighbor lists](@ref) interface. 

Note that, in general, neighbor lists are used to compute other pairwise dependent properties, and these can be, in principle, computed directly with `CellListMap` without the need to compute or store the lists of neighbors. 

### Lower level interface

The [Low level interface](@ref) allows the customization and optimization of very demanding calculations (although the PeriodicSystems interface does not have any performance limitation and is easier to use).

