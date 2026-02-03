#
# Functions to copy, reset and reduce output variables, that must be implemented
# by the user for custom output types.
#

"""
    copy_output(x)

Defines how the `output` variable is copied. Identical to `Base.copy(x)`
and implemented for the types in `$(SupportedTypes)`.

Other custom output types must have their `copy_output` method implemented.

# Example

```julia
using CellListMap
# Custom data type
struct A x::Int end
# Custom output type (array of A)
output = [ A(0) for _ in 1:100 ]
# How to copy an array of `A`
CellListMap.copy_output(v::Vector{A}) = [ x for x in v ]

# Alternatively, in this case, one could have defined:
Base.copy(a::A) = a
CellListMap.copy_output(v::Vector{A}) = copy(v)
```

The user must guarantee that the copy is independent of the original array.
For many custom types it is possible to define 
```
CellListMap.copy_output(v::Vector{T}) where {T<:CustomType} = deepcopy(v)
```

"""
function copy_output(x)
    throw(
        ArgumentError(
            """\n
                No method matching `copy_output($(typeof(x)))`

                Please implement a method 
               
                CellListMap.copy_output(x::$(typeof(x)))

                with an appropriate way to copy the required output variable. Many times just
                defining `CellListMap.copy_output(x::$(typeof(x))) = deepcopy(x)` is ok. 
            """
        )
    )
end
copy_output(x::T) where {T <: SupportedTypes} = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T} = T[copy_output(el) for el in x]

"""
    reset_output(x)
    reset_output!(x)

Function that defines how to reset (or zero) the `output` variable. For `$(SupportedTypes)` it is 
implemented as `zero(x)`.

Other custom output types must have their `reset_output!` method implemented. 

The function *must* return the variable itself. If it is immutable,
a new instante of the variable must be created, with the reset value. 

!!! note
    By default, if
    `reset_output!` is defined for one element type, `reset_output!` is defined for arrays of that type
    by calling `reset_output!` for each element of the array.  The user must overload the `reset_output!` 
    function for the custom type array if that is not the desired behavior.

`reset_output` and `reset_output!` are aliases, and by convention `reset_output!` is preferred for mutable types.

# Example

In this example, we define a `reset_output` function that will set to `+Inf` the
minimum distance between particles (not always resetting means zeroing).

```jldoctest
julia> using CellListMap

julia> struct MinimumDistance d::Float64 end

julia> CellListMap.reset_output(x::MinimumDistance) = MinimumDistance(+Inf)

julia> x = MinimumDistance(1.0)
MinimumDistance(1.0)

julia> CellListMap.reset_output(x)
MinimumDistance(Inf)
```

See the `reducer` help entry for a complete example of how to use `reset_output`.

"""
function reset_output!(x)
    throw(
        ArgumentError(
            """\n
                No method matching `reset_output!($(typeof(x)))`

                Please add a method 
                
                CellListMap.reset_output!(x::$(typeof(x)))
                
                with the appropriate way to reset (zero) the data of the output variables.

                The reset_output! methods **must** return the output variable to
                conform with the interface, even if the variable is mutable. 
            """
        )
    )
end
reset_output!(x::T) where {T <: SupportedTypes} = zero(x)
function reset_output!(x::AbstractVecOrMat{T}) where {T}
    for i in eachindex(x)
        x[i] = reset_output!(x[i])
    end
    return x
end
const reset_output = reset_output!

"""
    reducer(x,y)
    reducer!(x,y)

Defines how to reduce (combine, or merge) to variables computed in parallel
to obtain a single instance of the variable with the reduced result. 

`reducer` and `reducer!` are aliases, and `reducer!` is preferred, by convention
for mutating functions.

The most common `reducer` is the sum, and this is how it is implemented for
`$(SupportedTypes)`. For example, when computing energies, or forces,
the total energy is the sum of the energies. The force on one particle is the sum of the
forces between the particle and every other particle. Thus, the implemented reducer is
the sum: 

```
reducer(x,y) = +(x,y)
```

However, in  many cases, reduction must be done differently. For instance, if the minimum
distance between particles is to be computed, it is interesting to define a custom type
and associated reducer. For example:

```
struct MinimumDistance d::Float64 end
reducer(x::MinimumDistance, y::MinimumDistance) = MinimumDistance(min(x.d, y.d))
```

The overloading of `reducer` allows the use of parallel computations for custom, 
complex data types, containing different types of variables, fields, or sizes.

The appropriate behavior of the reducer should be carefully inspected by the user
to avoid spurious results. 

# Example

In this example we show how to obtain the minimum distance among argon atoms
in a simulation box.

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
julia> using CellListMap, PDBTools

julia> positions = coor(read_pdb(CellListMap.argon_pdb_file));

julia> struct MinimumDistance d::Float64 end # Custom output type

julia> CellListMap.copy_output(d::MinimumDistance) = MinimumDistance(d.d) # Custom copy function for `Out`

julia> CellListMap.reset_output(d::MinimumDistance) = MinimumDistance(+Inf) # How to reset an array with elements of type `MinimumDistance`

julia> CellListMap.reducer(md1::MinimumDistance, md2::MinimumDistance) = MinimumDistance(min(md1.d, md2.d)) # Custom reduction function

julia> # Construct the system
       sys = ParticleSystem(;
           positions = positions,
           unitcell = [21,21,21],
           cutoff = 8.0,
           output = MinimumDistance(+Inf),
       );

julia> # Obtain the minimum distance between atoms:
       pairwise!((pair,output) -> pair.d < output.d ? MinimumDistance(pair.d) : output, sys)
MinimumDistance(2.1991993997816563)
```

"""
function reducer!(x, y)
    throw(
        ArgumentError(
            """\n
                No method matching `reducer!($(typeof(x)),$(typeof(y)))`

                Please implement a method 
                
                CellListMap.reducer(x::$(typeof(x)),y::$(typeof(y)))
                
                with the appropriate way to combine two instances of the type (summing, keeping
                the minimum, etc), such that threaded computations can be reduced.

            """
        )
    )
end
reducer!(x::T, y::T) where {T <: SupportedTypes} = +(x, y)
const reducer = reducer!