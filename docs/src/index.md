# CellListMap

This package is for computing interactions or any other property that is dependent on the distances between pairs of two- or three-dimensional particles, within a cutoff. It maps a function to be computed pairwise using cell lists, using periodic boundary conditions of any type. Parallel and serial implementations can be used. 

It allows the fast computation of any quantity from the pairs that are within the desired cutoff, for example an average distance or an histogram of distances, forces, potentials, minimum distances, etc., as the examples below illustrate. This is done by passing the function to be evaluated as a parameter of the `map_pairwise!` function. 

## Installation

```julia-repl
julia> import Pkg

julia> Pkg.add("CellListMap")

```

## Overview

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