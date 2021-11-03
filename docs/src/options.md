# Additional options

## Input coordinates as matrices

For compatibility with other software, the input coordinates can be provided as matrices. The matrices must have dimensions `(2,N)` or `(3,N)`, where `N` is the number of particles (because Julia is column-major, thus this has the same memory layout of an array of length `N` of static vectors). 

For example:
```julia-repl
julia> x = rand(3,100);

julia> box = Box([1,1,1],0.1);

julia> cl = CellList(x,box)
CellList{3, Float64}
  100 real particles.
  99 cells with real particles.
  162 particles in computing box, including images.

julia> map_pairwise!((x,y,i,j,d2,n) -> n += 1, 0, box, cl) # count neighbours
23
```

## Non-allocating type conversion 

Internally, `CellListMap`  will perform operations on real numbers. However, sometimes the coordinates are provided using some custom type for which conversion to standard floats is not defined. For example, it is common to use coordinates with units:

```julia-repl
julia> using Unitful

julia> x = [ rand(3)u"nm" for i in 1:100 ];

julia> x[1]
3-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 0.5555681432511039 nm
 0.3112134334494392 nm
 0.6849761663523335 nm
```

 In order to use the type of coordinates without allocations and other complications in `CellListMap`, just overload the `CellListMap.strip_value` function such that it converts a value of the given type to a float. For example, the `Unitful` package provides the `ustrip` function to convert values of the `Quantity` type to floats, by removing the units. We define, then:

```julia-repl
julia> CellListMap.strip_value(x::Quantity) = Unitful.ustrip(x)
```

such that it converts a single value of type `Quantity` to a standard float:
```julia-repl
julia> CellListMap.strip_value(x[1][1])
0.5555681432511039
```

With that, the `Unitful` quantities can be passed to `CellList` without modification:

```julia-repl
julia> box = Box([1,1,1],0.1);

julia> CellList(x,box)
CellList{3, Float64}
  100 real particles.
  93 cells with real particles.
  167 particles in computing box, including images.
```








