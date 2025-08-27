
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/CellListMap.jl/stable)
[![Build Status](https://github.com/m3g/CellListMap.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/m3g/CellListMap.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/m3g/CellListMap.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/CellListMap.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Package Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FCellListMap&query=total_requests&label=Downloads)](http://juliapkgstats.com/pkg/CellListMap)

<p align=center>
<img src=https://raw.githubusercontent.com/m3g/CellListMap.jl/main/docs/src/assets/logo.svg>
</p>

# CellListMap.jl

This package is for computing interactions or any other property that is dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. It maps a function to be computed pairwise using cell lists, using periodic boundary conditions of any type. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example pairwise potentials and forces, neighbor lists, minimum distances, an average distance or an histogram of distances, etc. This is done by passing the function to be evaluated as a parameter of the `map_pairwise!` function. 

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

Download and install Julia for your platform from [this http url](https://julialang.org/downloads/). Version 1.9 or greater is required.

Install it as usual for registered Julia packages:

```julia
julia> import Pkg

julia> Pkg.add("CellListMap")
```

## Brief overview

The main function is `map_parwise!`: 

If the analysis is performed on the pairs of a single vector `x` (`n*(n-1)/2` pairs), the function can be called with:
```julia
map_pairwise!(f::Function,output,box::Box,cl::CellList)
```
while if two distinct sets of points are provided (`n*m` pairs), it is called with:
```julia
map_pairwise!(f::Function,output,box::Box,cl::CellListPair)
```
where the `cl` variable of type `CellList` or `CellListPair` contains the cell lists built from the coordinates of the system, and `box` contains the system box properties.

These functions will run over every pair of particles which are closer than `box.cutoff` and compute the (squared) Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a user defined function `f` of the coordinates of the two particles will be computed. 

The function `f` receives six arguments as input: 
```julia
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it.  The squared distance `d2` is computed   internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used for the desired computation. Thus, the function `f` that is passed to `map_pairwise!` must be always of the form:
```julia
function f(x,y,i,j,d2,output)
    # update output
    return output
end
```
and the user can define more or less parameters or additional data required to compute the function using closures, as shown in the examples.

Parallel calculations are the default if more than one thread is available. Use `parallel=false` as an optional argument to `map_pairwise!` to run the serial version instead.

## Some benchmarks

The goal here is to provide a good implementation of cell lists. We compare it with the implementation of the nice cython/python [halotools](https://github.com/astropy/halotools) package, in the computation of an histogram of mean pairwise velocities. 

<img src=https://github.com/lmiq/PairVelocities/blob/main/data/cd_v0.8.27-DEV.png>

<img src=https://github.com/lmiq/PairVelocities/blob/main/data/cv_v0.8.27-DEV.png>

The full test is available [at this](https://github.com/lmiq/PairVelocities) repository. And we kindly thank [Carolina Cuesta](https://github.com/florpi) for providing the example. These benchmarks were run on an Intel(R) Core(TM) i7-12700, using 8 CPUs. 

## Citation

If you use this software and need to cite it, please use the following reference:

L. Martínez, *CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff.* Computer Physics Communications, 279, 108452, 2022. https://doi.org/10.1016/j.cpc.2022.108452



