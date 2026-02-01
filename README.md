
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/CellListMap.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/CellListMap.jl/dev)
[![Build Status](https://github.com/m3g/CellListMap.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/m3g/CellListMap.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/m3g/CellListMap.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/CellListMap.jl)
[![Aqua QA](https://JuliaTesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Package Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FCellListMap&query=total_requests&label=Downloads)](http://juliapkgstats.com/pkg/CellListMap)

<p align=center>
<img src="https://m3g.github.io/CellListMap.jl/stable/assets/logo.svg">
</p>

# CellListMap.jl

This package is for computing interactions or any other property that is dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. It maps a function to be computed pairwise using cell lists, using periodic boundary conditions of any type. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example pairwise potentials and forces, neighbor lists, minimum distances, an average distance or an histogram of distances, etc. This is done by passing the function to be evaluated as a parameter of the `foreachneighbor!` function, which receives a `NeighborPair` struct for each pair.

The user guide provides direct examples of each of these applications. 

<h3>
<br>
<p align=center>
USER GUIDE: <br> 
<a href=https://m3g.github.io/CellListMap.jl>https://m3g.github.io/CellListMap.jl</a>
</p>
<br>
</h3>

## Installation

Download and install Julia for your platform from [this http url](https://julialang.org/downloads/). 
Version 1.10 or greater is required.

Install it as usual for registered Julia packages:

```julia
julia> import Pkg

julia> Pkg.add("CellListMap")
```

## Brief overview

The main function is `foreachneighbor!`: 

If the analysis is performed on the pairs of a single vector `x` (`n*(n-1)/2` pairs), the function can be called with:
```julia
foreachneighbor!(f::Function,output,box::Box,cl::CellList)
```
while if two distinct sets of points are provided (`n*m` pairs), it is called with:
```julia
foreachneighbor!(f::Function,output,box::Box,cl::CellListPair)
```
where the `cl` variable of type `CellList` or `CellListPair` contains the cell lists built from the coordinates of the system, and `box` contains the system box properties.

These functions will run over every pair of particles which are closer than `box.cutoff` and compute the (squared) Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a user defined function `f` will be called with a `NeighborPair` struct and the output variable. The `NeighborPair` struct contains:

- `pair.i`, `pair.j`: indices of the particles
- `pair.x`, `pair.y`: coordinates (minimum-image adjusted)
- `pair.d2`: squared distance
- `pair.d`: distance (computed lazily)

The function `f` must have the form:
```julia
function f(pair, output)
    # update output using pair.i, pair.j, pair.x, pair.y, pair.d2, pair.d
    return output
end
```
and the user can define more or less parameters or additional data required to compute the function using closures, as shown in the examples.

Parallel calculations are the default if more than one thread is available. Use `parallel=false` as an optional argument to `foreachneighbor!` to run the serial version instead.

## Some benchmarks

The goal here is to provide a good implementation of cell lists. We compare it with the implementation of the nice cython/python [halotools](https://github.com/astropy/halotools) package, in the computation of an histogram of mean pairwise velocities. 

<img src="https://m3g.github.io/CellListMap.jl/stable/assets/b_cd.png">

<img src="https://m3g.github.io/CellListMap.jl/stable/assets/b_cv.png">

The full test is available [at this](https://github.com/lmiq/PairVelocities) repository. And we kindly thank [Carolina Cuesta](https://github.com/florpi) for providing the example. These benchmarks were run on an Intel(R) Core(TM) i7-12700, using 8 CPUs. 

## Citation

If you use this software and need to cite it, please use the following reference:

L. Martínez, *CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff.* Computer Physics Communications, 279, 108452, 2022. https://doi.org/10.1016/j.cpc.2022.108452



