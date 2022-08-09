# CellListMap.jl

This package is for computing interactions or any other property dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. It maps a function to be computed pairwise using cell lists, using periodic boundary conditions of any type. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example an average distance or an histogram of distances, forces, potentials, minimum distances, etc., as the examples below illustrate. This is done by passing the function to be evaluated as a parameter of the `map_pairwise!` function. 

## Installation

```julia-repl
julia> import Pkg

julia> Pkg.add("CellListMap")
```

## Overview

### High level interface for periodic system

Since version `0.7.22`, a new simpler, higher level interface was introduced, that will facilitate the use of `CellListMap` without any loss in performance. The new interface is flexible enough for the majority of applications. It may become the default interface in the future. See the [PeriodicSystems interface] menu for details. 

### Standard lower level interface

The main function is `map_parwise!` (or `map_pairwise`): 

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

### Mutable and immutable outputs

`map_pairwise!` and `map_pairwise` (with the bang, or not) are aliases of the same function, which always returns the result value. It is a convention in Julia that functions ending with the `!` mutate the arguments, while those without do not. Here, this behavior is dependent on the type of input. If the output variable is immutable, its value won't be mutated, and the assignment of the result to the output value depends on explicit assignment. In these cases, it is customary to use the `map_pairwise` (without `!`) function name:
```julia
output = map_pairwise(function, output0, box, cl)
```
where `output0` represents the initial value of the immutable `output`. When, on the contrary, the output is a mutable variable (an array, for example), the `map_pairwise!` version is preferred for code clarity, and the reassignment is not needed (nor recommendable): 
```julia
output = zeros(10) # example of mutable output
map_pairwise!(function, output, box, cl)
```
