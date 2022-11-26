module PeriodicSystems

using TestItems
using DocStringExtensions
using StaticArrays

import ..CellListMap
import ..CellListMap: INTERNAL
import ..CellListMap: Box, update_box, unitcelltype
import ..CellListMap: CellListPair, Swapped, NotSwapped

export PeriodicSystem
export map_pairwise!, map_pairwise
export update_cutoff!
export update_unitcell!
export unitcelltype

export copy_output
export resize_output!
export reset_output!, reset_output
export reducer!, reducer

# Types of variables that have support for multi-threading without having 
# to explicit add methods to copy_output, reset_output!, and reducer functions.
const SupportedTypes = Union{Number,SVector,FieldVector}

"""

```
PeriodicSystem( 
    xpositions::Vector{<:AbstractVector},
    #or
    xpositions::Vector{<:AbstractVector},
    ypositions::Vector{<:AbstractVector},
    # and
    unitcell::AbstractVecOrMat,
    cutoff::Number,
    output::Any;
    output_name::Symbol,
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
```

Function that sets up the `PeriodicSystem` type given the positions of
the particles.

- Positions can be provided as vectors of 2D or 3D vectors 
  (preferentially static vectors from `StaticArrays`).

- If only the `xpositions` array is provided, a single set of coordinates 
  is considered, and the computation will be mapped for the `N(N-1)` 
  pairs of this set. 

- If the `xpositions` and `ypositions` arrays of coordinates are provided, 
  the computation will be mapped to the `N×M` pairs of particles, 
  being `N` and `M` the number of particles of each set of coordinates.

The unit cell (either a vector for `Orthorhombic` cells or a 
full unit cell matrix for `Triclinic` cells), the cutoff used for the
construction of the cell lists and the output variable of the calculations.

`output_name` can be set to a symbol that best identifies the output variable.
For instance, if `output_name=:forces`, the forces can be retrieved from the
structure using the `system.forces` notation.

The `parallel` and `nbatches` flags control the parallelization scheme of
computations (see https://m3g.github.io/CellListMap.jl/stable/parallelization/#Number-of-batches)).
By default the parallelization is turned on and `nbatches` is set with heuristics
that may provide good efficiency in most cases.

# Example

In these examples, we compute the sum of the squared distances between
the particles that are within the cutoff:

## Single set of particles

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> positions = rand(SVector{3,Float64}, 100)
       unitcell = [1,1,1]
       cutoff = 0.1
       output = 0.0;

julia> sys = PeriodicSystem(xpositions = positions, unitcell= [1.0,1.0,1.0], cutoff = 0.1, output = 0.0)
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [12, 12, 12]
      computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
      Total number of cells = 1728
    CellListMap.CellList{3, Float64}
      100 real particles.
      96 cells with real particles.
      164 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 8
    Type of output variable: Float64

julia> map_pairwise!((x,y,i,j,d2,output) -> output += d2, sys)
0.14556244865996287
```
## Two sets of particles

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> xpositions = rand(SVector{3,Float64}, 100)
       ypositions = rand(SVector{3,Float64}, 1000)
       unitcell = [1,1,1]
       cutoff = 0.1
       output = 0.0;

julia> sys = PeriodicSystem(
           xpositions = xpositions, 
           ypositions = ypositions,
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           )
PeriodicSystem2 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [12, 12, 12]
      computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
      Total number of cells = 1728
    CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64, CellListMap.Swapped}
       1000 particles in the reference vector.
       97 cells with real particles of target vector.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 8
    Type of output variable: Float64

julia> map_pairwise!((x,y,i,j,d2,output) -> output += d2, sys)
2.3568143238242314

```

"""
function PeriodicSystem(;
    positions::Union{Nothing,Vector{<:AbstractVector}}=nothing,
    xpositions::Union{Nothing,Vector{<:AbstractVector}}=nothing,
    ypositions::Union{Nothing,Vector{<:AbstractVector}}=nothing,
    unitcell::AbstractVecOrMat,
    cutoff::Number,
    output::Any,
    output_name::Symbol=:output,
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0),
    lcell=1,
    autoswap::Bool=true,
)
    if !isnothing(positions) && isnothing(xpositions)
        xpositions = positions
    elseif !isnothing(positions) && !isnothing(xpositions)
        throw(ArgumentError(
            """Either define `positions` OR `xpositions`, they are aliases one to the other."""
        ))
    end
    if !isnothing(xpositions) && isnothing(ypositions)
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(xpositions, _box; parallel=parallel, nbatches=nbatches)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem1{output_name}(xpositions, output, _box, _cell_list, _output_threaded, _aux, parallel)
    elseif !isnothing(xpositions) && !isnothing(ypositions)
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(xpositions, ypositions, _box; parallel=parallel, nbatches=nbatches, autoswap=autoswap)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem2{output_name}(xpositions, ypositions, output, _box, _cell_list, _output_threaded, _aux, parallel)
    else
        throw(ArgumentError(
            """Either define `xpositions` OR (`xpositions` AND `ypositions`), to build systems for self- or cross-pair computations, respectively."""
        ))
    end
    return sys
end


# Abstract type only for cleaner dispatch
abstract type AbstractPeriodicSystem{OutputName} end

import Base: getproperty, propertynames
getproperty(sys::AbstractPeriodicSystem, s::Symbol) = getproperty(sys, Val(s))
getproperty(sys::AbstractPeriodicSystem, s::Val{S}) where {S} = getfield(sys, S)
# publi properties
getproperty(sys::AbstractPeriodicSystem, ::Val{:unitcell}) = getfield(getfield(getfield(sys, :_box), :unit_cell), :matrix)
getproperty(sys::AbstractPeriodicSystem, ::Val{:cutoff}) = getfield(getfield(sys, :_box), :cutoff)
getproperty(sys::AbstractPeriodicSystem{OutputName}, ::Val{OutputName}) where {OutputName} = getfield(sys, :output)
propertynames(sys::AbstractPeriodicSystem{OutputName}) where {OutputName} =
    (:xpositions, :ypositions, :unitcell, :cutoff, :positions, :output, :parallel, OutputName)

import Base: setproperty!
# public properties
setproperty!(sys::AbstractPeriodicSystem, s::Symbol, x) = setproperty!(sys, Val(s), x)
setproperty!(sys::AbstractPeriodicSystem, ::Val{:unitcell}, x) = update_unitcell!(sys, x)
setproperty!(sys::AbstractPeriodicSystem, ::Val{:cutoff}, x) = update_cutoff!(sys, x)
setproperty!(sys::AbstractPeriodicSystem, ::Val{:parallel}, x) = setfield!(sys, :parallel, x)
# private properties
setproperty!(sys::AbstractPeriodicSystem, ::Val{:_box}, x) = setfield!(sys, :_box, x)
setproperty!(sys::AbstractPeriodicSystem, ::Val{:_cell_list}, x) = setfield!(sys, :_cell_list, x)
setproperty!(sys::AbstractPeriodicSystem, ::Val{:output}, x) = setfield!(sys, :output, x)

"""
    unitcelltype(sys::AbstractPeriodicSystem)

Returns the type of a unitcell from the `PeriodicSystem` structure.

"""
unitcelltype(sys::AbstractPeriodicSystem) = unitcelltype(sys._box)

@testitem "PeriodicSystems properties" begin

    using CellListMap.PeriodicSystems
    using StaticArrays
    sys = PeriodicSystem(
        positions=rand(SVector{3,Float64}, 1000),
        cutoff=0.1,
        unitcell=[1, 1, 1],
        output=0.0,
        output_name=:test
    )
    @test length(sys.positions) == 1000
    @test sys.cutoff == 0.1
    @test sys.unitcell == @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test sys.output == 0
    @test sys.test == 0
    @test sys.parallel == true

    sys.parallel = false
    @test sys.parallel == false
    sys.cutoff = 0.2
    @test sys.cutoff == 0.2
    sys.positions[1] = SVector(0.0, 0.0, 0.0)
    @test sys.positions[1] == SVector(0.0, 0.0, 0.0)
    sys.unitcell = [1.2, 1.2, 1.2]
    @test sys.unitcell == @SMatrix [1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2]

    # test the construction with pathologically few particles
    for x in [ SVector{3,Float64}[], Vector{Float64}[], [rand(SVector{3,Float64})], [rand(3)] ]
        sys = PeriodicSystem(
            positions=x,
            cutoff=0.1,
            unitcell=[1, 1, 1],
            output=0.0,
            output_name=:test
        )
        @test PeriodicSystems.map_pairwise((x,y,i,j,d2,out) -> out += d2, sys) == 0.0
    end

    # unitcell type
    x = rand(SVector{3,Float64},100)
    @test unitcelltype(PeriodicSystem(positions=x, cutoff=0.1, unitcell=[1,1,1], output=0.0)) == OrthorhombicCell
    @test unitcelltype(PeriodicSystem(positions=x, cutoff=0.1, unitcell=[1 0 0; 0 1 0; 0 0 1], output=0.0)) == TriclinicCell

end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `map_pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). 

The `xpositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrive or mutate `xpositions`, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fileds of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `PeriodicSystem1` structure
is done through the `PeriodicSystem(;xpositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct PeriodicSystem1{OutputName,V,O,B,C,A} <: AbstractPeriodicSystem{OutputName}
    xpositions::Vector{V}
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
PeriodicSystem1{OutputName}(v::Vector{V}, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool) where {OutputName,V,O,B,C,A} =
    PeriodicSystem1{OutputName,V,O,B,C,A}(v, o, b, c, vo, a, p)
getproperty(sys::PeriodicSystem1, ::Val{:positions}) = getfield(sys, :xpositions)

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `map_pairwise!` computations,
for systems with two set of positions (thus, replacing the loops over `N×M` 
pairs of particles, being `N` and `M` the number of particles of each set).

The `xpositions`, `ypositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrive or mutate positions, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fileds of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `PeriodicSystem1` structure
is done through the `PeriodicSystem(;xpositions, ypositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct PeriodicSystem2{OutputName,V,O,B,C,A} <: AbstractPeriodicSystem{OutputName}
    xpositions::Vector{V}
    ypositions::Vector{V}
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
PeriodicSystem2{OutputName}(vx::Vector{V}, vy::Vector{V}, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool) where {OutputName,V,O,B,C,A} =
    PeriodicSystem2{OutputName,V,O,B,C,A}(vx, vy, o, b, c, vo, a, p)

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "PeriodicSystem1{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println()
    show(io_sub, mime, sys._cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.nbatches)
    print("\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "PeriodicSystem2{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println()
    show(io_sub, mime, sys._cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.target.nbatches)
    print("\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

#
# Functions to copy, reset and reduce output variables, that must be implemented
# by the user for custom output types.
#
"""

```
copy_output(x)
```

Function that defines how the `output` variable is copied. Identical to `Base.copy(x)`
and implemented for the types in `$(SupportedTypes)`.

Other custom output types must have their `copy_output` method implemented.

# Example

```julia
using CellListMap.PeriodicSystems
# Custom data type
struct A x::Int end
# Custom output type (array of A)
output = [ A(0) for _ in 1:100 ]
# How to copy an array of `A`
PeriodicSystems.copy_output(v::Vector{A}) = [ x for x in v ]

# Alternativelly, in this case, one could have defined:
Base.copy(a::A) = a
PeriodicSystems.copy_output(v::Vector{A}) = copy(v)
```

The user must guarantee that the copy is independent of the original array.
For many custom types it is possible to define 
`PeriodicSystems.copy_output(v::Vector{T}) where {T<:CustomType} = deepcopy(v)`.

"""
function copy_output(x)
    error("""
        MethodError: no method matching `copy_output($(typeof(x)))`

        Please implement a method 
       
        PeriodicSystems.copy_output(x::$(typeof(x)))

        with an appropriate way to copy the required output variable. Many times just
        defining `output_copy(x::$(typeof(x))) = deepcopy(x)` is ok. 
    """)
end
copy_output(x::T) where {T<:SupportedTypes} = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T} = T[copy_output(el) for el in x]

"""

```
reset_output!(x)
```

Function that defines how to reset (or zero) the `output` variable. For `$(SupportedTypes)` it is 
implemented as `zero(x)`, and for `AbstractVecOrMat` containers of `Number`s or `SVector`s
it is implemented as `fill!(x, zero(eltype(x))`.

Other custom output types must have their `reset_output!` method implemented.

The `reset_output!` function *must* return the `output` variable, being it mutable
or immutable. `reset_output` is an alias for `reset_output!` that can be used for consistency if the
`output` variable is immutable.

# Example

In this example, we define a `reset_output` function that will set to `+Inf` the
minimum distance between particles (not always resetting means zeroing).

```julia
# Custom data type
struct MinimumDistance d::Float64 end
# How to reset the minimum distance
PeriodicSystems.reset_output!(x::MinimumDistance) = MinimumDistance(+Inf)
```

"""
function reset_output!(x)
    error("""
        MethodError: no method matching `reset_output!($(typeof(x)))`

        Please add a method 
        
        PeriodicSystems.reset_output!(x::$(typeof(x)))
        
        with the appropriate way to reset (zero) the data of the output variables.

        The `reset_output!` methods **must** return the output variable to
        conform with the interface, even if the variable is mutable. 
    """)
end
reset_output!(x::T) where {T<:SupportedTypes} = zero(x)
reset_output!(x::AbstractVecOrMat{T}) where {T} = fill!(x, reset_output!(x[begin]))
const reset_output = reset_output!

"""

```
_reset_all_output!(output, output_threaded)
```

$(INTERNAL)

Function that resets the output variable and the threaded copies of it.

"""
function _reset_all_output!(output, output_threaded)
    output = reset_output!(output)
    for i in eachindex(output_threaded)
        output_threaded[i] = reset_output!(output_threaded[i])
    end
    return output
end

"""

```
reducer!(x,y)
```

Function that defines how to reduce (combine, or merge) to variables computed in parallel
to obtain a single instance of the variable with the reduced result.

The most commont `reducer` is the sum, and this is how it is implemented for
`$(SupportedTypes)`. For example, when computin energies, or forces,
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

The appropriate behavior of the reducer should be carefuly inspected by the user
to avoid spurious results. 

"""
function reducer!(x, y)
    error("""
        MethodError: no method matching `reducer!($(typeof(x)),$(typeof(x)))`

        Please implement a method 
        
        PeriodicSystems.reducer(x::$(typeof(x)),y::$(typeof(x)))
        
        with the appropriate way to combine two instances of the type (summing, keeping
        the minimum, etc), such that threaded computations can be reduced.
    """)
end
reducer!(x::T, y::T) where {T<:SupportedTypes} = +(x, y)
const reducer = reducer!

"""

```
reduce_output!(reducer::Function, output, output_threaded)
```

$(INTERNAL)

# Extended help

Function that defines how to reduce the vector of `output` variables, after a threaded
computation. This function is implemented for `output` variables that are numbers, 
and vectors or arrays of number of static arrays, as the sum of the values of the 
threaded computations, which is the most common application, found in computing
forces, energies, etc. 

It may be interesting to implement custom `PeriodicSystems.reduce_output!` function for other types 
of output variables, considering:

- The arguments of the function must be the return `output` value and a vector 
  `output_threaded` of `output` variables, which is created (automatically) by copying
  the output the number of times necessary for the multi-threaded computation. 

- The function *must* return the `output` variable, independently of it being mutable
  or immutable.

`reduce_output` is an alias to `reduce_output!` that can be used for consistency if the `output`
variable is immutable.

# Example

In this example we show how to obtain the minimum distance between two sets
of particles. This requires a custom reduction function.

```julia
using CellListMap.PeriodicSystems, StaticArrays
# Custom output type
struct MinimumDistance
    d::Float64
end
# Custom copy function for `Out`
PeriodicSystems.copy_output(d::MinimumDistance) = MinimumDistance(d.d)
# How to reset an array with elements of type `MinimumDistance`
PeriodicSystems.reset_output!(d::MinimumDistance) = MinimumDistance(+Inf)
# Custom reduction function (keep the minimum distance)
function PeriodicSystems.reduce_output!(
    output::MinimumDistance, 
    output_threaded::Vector{MinimumDistance}
)
    output = reset_output!(output)
    for i in eachindex(output_threaded)
        if output_threaded[i].d < output.d
            output = output_threaded[i]
        end
    end
    return output
end
# Construct the system
sys = PeriodicSystem(;
    xpositions = rand(SVector{3,Float64}, 1000),
    ypositions = rand(SVector{3,Float64}, 1000),
    unitcell = [1,1,1],
    cutoff = 0.1,
    output = MinimumDistance(+Inf),
)

# Obtain the minimum distance between the sets
map_pairwise!((x,y,i,j,d2,output) -> sqrt(d2) < output.d ? MinimumDistance(sqrt(d2)) : output, sys)
# will output something like: MinimumDistance(0.00956913034767034)
```

"""
function reduce_output!(reducer::Function, output::T, output_threaded::Vector{T}) where {T}
    output = reset_output!(output)
    for ibatch in eachindex(output_threaded)
        output = reducer(output, output_threaded[ibatch])
    end
    return output
end
function reduce_output!(
    reducer::Function,
    output::AbstractVecOrMat{T},
    output_threaded::Vector{<:AbstractVecOrMat{T}}
) where {T}
    output = reset_output!(output)
    for ibatch in eachindex(output_threaded)
        for i in eachindex(output, output_threaded[ibatch])
            output[i] = reducer(output[i], output_threaded[ibatch][i])
        end
    end
    return output
end

"""

```
resize_output!(sys::AbstractPeriodicSystem, n::Int)
```

Resizes the output array and the auxiliary output arrays used
for multithreading, if needed because of the system change.

"""
function resize_output!(sys::AbstractPeriodicSystem, n::Int)
    resize!(sys.output, n)
    for i in eachindex(sys._output_threaded)
        resize!(sys._output_threaded[i], n)
    end
    return sys
end

#
# Function used to update the properties of the systems
#
"""
```
update_unitcell!(system, unitcell)
```

Function to update the unit cell of the system. The `unicell` must be of the 
same type (`OrthorhombicCell` or `TriclinicCell`) of the original `system` 
(changing the type of unit cell requires reconstructing the system).

The `unitcell` can be a `N×N` matrix or a vector of dimension `N`, where
`N` is the dimension of the sytem (2D or 3D).

This function can be used to update the system geometry in iterative schemes,
where the size of the simulation box changes during the simulation.

# Example

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> sys = PeriodicSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           );

julia> update_unitcell!(sys, [1.2, 1.1, 1.0])
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.2, 0.0, 0.0; 0.0, 1.1, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [13, 13, 12]
      computing cell sizes = [0.11, 0.1, 0.1] (lcell: 1)
      Total number of cells = 2028
    CellListMap.CellList{3, Float64}
      1000 real particles.
      633 cells with real particles.
      1703 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 12
    Type of output variable: Float64
```

"""
function update_unitcell!(sys, unitcell)
    sys._box = update_box(sys._box; unitcell=unitcell)
    return sys
end

@testitem "update_unitcell!" begin
    using BenchmarkTools
    using LinearAlgebra: diag
    using StaticArrays
    using CellListMap.PeriodicSystems
    x = rand(SVector{3,Float64}, 1000)
    sys1 = PeriodicSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_unitcell!(sys1, SVector(2, 2, 2))
    @test diag(sys1.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys1, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys2 = PeriodicSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_unitcell!(sys2, SVector(2, 2, 2))
    @test diag(sys2.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys2, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == 0
end

"""
```
update_cutoff!(system, cutoff)
```

Function to update the `cutoff`` of the system. 

This function can be used to update the system geometry in iterative schemes.

# Example

```julia-repl
julia> using CellListMap.PeriodicSystems, StaticArrays

julia> sys = PeriodicSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           );

julia> update_cutoff!(sys, 0.2)
PeriodicSystem1 of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.2
      number of computing cells on each dimension = [7, 7, 7]
      computing cell sizes = [0.2, 0.2, 0.2] (lcell: 1)
      Total number of cells = 343
    CellListMap.CellList{3, Float64}
      1000 real particles.
      620 cells with real particles.
      1746 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 12
    Type of output variable: Float64
```
"""
function update_cutoff!(sys, cutoff)
    sys._box = update_box(sys._box; cutoff=cutoff)
    return sys
end

@testitem "update_cutoff!" begin
    using BenchmarkTools
    using StaticArrays
    using CellListMap.PeriodicSystems
    x = rand(SVector{3,Float64}, 1000)
    sys1 = PeriodicSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_cutoff!(sys1, 0.2)
    @test sys1.cutoff == 0.2
    a = @ballocated update_cutoff!($sys1, 0.1) evals = 1 samples = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys2 = PeriodicSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_cutoff!(sys2, 0.2)
    @test sys2.cutoff == 0.2
    a = @ballocated update_cutoff!($sys2, 0.1) evals = 1 samples = 1
    @test a == 0
end

"""

```
UpdatePeriodicSystem!
```

$(INTERNAL)

Updates the cell lists for periodic systems.

"""
function UpdatePeriodicSystem!(sys::PeriodicSystem1, update_lists::Bool=true)
    if update_lists
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel=sys.parallel
        )
    end
    return sys
end

function _update_ref_positions!(cl::CellListPair{V,N,T,Swap}, sys) where {V,N,T,Swap <: NotSwapped} 
    resize!(cl.ref, length(sys.xpositions))
    cl.ref .= sys.xpositions
end
function _update_ref_positions!(cl::CellListPair{V,N,T,Swap}, sys) where {V,N,T,Swap <: Swapped} 
    error("update_lists === false requires autoswap == false for 2-set systems.")
end
function UpdatePeriodicSystem!(sys::PeriodicSystem2, update_lists::Bool=true)
    if update_lists
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys.ypositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel=sys.parallel
        )
    else
        # Always update the reference set positions (the cell lists of the target set are not updated)
        _update_ref_positions!(sys._cell_list, sys)
    end
    return sys
end

# this updates must be non-allocating in the serial case
@testitem "UpdatePeriodicSystem!" begin
    using BenchmarkTools
    using StaticArrays
    using CellListMap.PeriodicSystems
    x = rand(SVector{3,Float64}, 1000)
    sys = PeriodicSystem(xpositions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated PeriodicSystems.UpdatePeriodicSystem!($sys) samples = 1 evals = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys = PeriodicSystem(xpositions=x, ypositions=y, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated PeriodicSystems.UpdatePeriodicSystem!($sys) samples = 1 evals = 1
    @test a == 0
end

"""

```
map_pairwise!(f::Function, system; show_progress = true, update_lists = true)
```

Function that maps the `f` function into all pairs of particles of
`system` that are found to be within the `cutoff`. 

The function `f` must be of the general form:
```
function f(x,y,i,j,d2,output)
    # operate on particle coordinates, distance and indexes
    # update output
    return output
end
```
where `x` and `y` are the coordinates (adjusted for the minimum
image) of the two particles involved, `i` and `j` their indices in the
original arrays of positions, `d2` their squared Euclidean distance,
and `output` the current value of the `output` variable. The `output`
variable must be updated within this function with the contribution
of the two particles involved. 

Thread-safety is taken care automatically in parallel executions.

`map_pairwise` is an alias to `map_pairwise!` for syntax consistency
when the `output` variable is immutable.

If `update_lists` is `false`, the cell lists will not be recomputed,
this may be useful for computing a different function from the same
coordinates.

# Example

In this example we compute the sum of `1/(1+d)` where `d` is the
distance between particles of a set, for `d < cutoff`. 

```julia-repl
julia> sys = PeriodicSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           );

julia> map_pairwise((x,y,i,j,d2,output) -> output += 1 / (1 + sqrt(d2)), sys)
1870.0274887950268
```

"""
function map_pairwise!(
    f::F,
    sys::AbstractPeriodicSystem;
    update_lists::Bool=true,
    show_progress::Bool=false
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded)
    UpdatePeriodicSystem!(sys, update_lists)
    sys.output = CellListMap.map_pairwise!(
        f, sys.output, sys._box, sys._cell_list;
        output_threaded=sys._output_threaded,
        parallel=sys.parallel,
        reduce=(output, output_threaded) -> reduce_output!(reducer, output, output_threaded),
        show_progress=show_progress
    )
    return sys.output
end
const map_pairwise = map_pairwise!

end # module PeriodicSystems
