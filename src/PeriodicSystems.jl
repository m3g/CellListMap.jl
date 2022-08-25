module PeriodicSystems

using DocStringExtensions
using StaticArrays

import ..CellListMap
import ..CellListMap: INTERNAL

export PeriodicSystem
export map_pairwise!, map_pairwise
export update_cutoff!
export update_unitcell!

export copy_output
export resize_output!
export reset_output!, reset_output
export reducer!, reducer

"""

```
PeriodicSystem( 
    positions::Vector{<:AbstractVector},
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

- If the `positions` array is provided, a single set of coordinates 
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

julia> sys = PeriodicSystem(positions = positions, unitcell= [1.0,1.0,1.0], cutoff = 0.1, output = 0.0)
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
)
    if !isnothing(positions) && (isnothing(xpositions) && isnothing(ypositions))
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(positions, _box; parallel=parallel, nbatches=nbatches)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem1{output_name}(positions, output, _box, _cell_list, _output_threaded, _aux, parallel)
    elseif isnothing(positions) && (!isnothing(xpositions) && !isnothing(ypositions))
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(xpositions, ypositions, _box; parallel=parallel, nbatches=nbatches)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem2{output_name}(xpositions, ypositions, output, _box, _cell_list, _output_threaded, _aux, parallel)
    else
        throw(ArgumentError(
            """Either define `positions` OR (`xpositions` AND `ypositions`), to build
               systems for self- or cross-pair computations, respectively.
            """))
    end
    return sys
end

function copy_to_vector(positions)
    if positions isa AbstractVector{<:AbstractVector}
        posvec = [ SVector(ntuple(i -> v[i], length(v))) for v in positions ]
    elseif positions isa AbstractMatrix
        posvec = [ SVector(ntuple(i -> v[i], length(v))) for v in eachcol(positions) ]
    end
    return posvec
end

# Abstrct type only for cleaner dispatch
abstract type AbstractPeriodicSystem{OutputName} end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `map_pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). 

The `positions`, `output`, and `parallel` fields are considered part of the API,
and you can retrive or mutate `positions`, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fileds of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `PeriodicSystem1` structure
is done through the `PeriodicSystem(;positions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct PeriodicSystem1{OutputName,V,O,B,C,A} <: AbstractPeriodicSystem{OutputName}
    positions::Vector{V}
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
PeriodicSystem1{OutputName}(v::Vector{V},o::O,b::B,c::C,vo::Vector{O},a::A,p::Bool) where {OutputName,V,O,B,C,A} = 
    PeriodicSystem1{OutputName,V,O,B,C,A}(v,o,b,c,vo,a,p)

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
mutable struct PeriodicSystem2{OutputName,V,B,C,O,A} <: AbstractPeriodicSystem{OutputName}
    xpositions::Vector{V}
    ypositions::Vector{V}
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
PeriodicSystem2{OutputName}(vx::Vector{V},vy::Vector{V},o::O,b::B,c::C,vo::Vector{O},a::A,p::Bool) where {OutputName,V,B,C,O,A} = 
    PeriodicSystem2{OutputName,V,B,C,O,A}(vx,vy,o,b,c,vo,a,p)

#
# This method of getproperty allows the user to access the output
# variable by a custom name given in `output_name`.
#
# This overloading causes some allocations, which do not impair
# perfomance in the typical use case, but maybe we will get resized
# of this in the future, changing the interface.
#
import Base: getproperty, propertynames
function getproperty(sys::AbstractPeriodicSystem{OutputName}, s::Symbol) where {OutputName}
    if s == OutputName
        return getfield(sys, :output)
    end
    return getfield(sys, s)
end
function propertynames(sys::AbstractPeriodicSystem{OutputName}, private::Bool=true) where {OutputName}
    names = fieldnames(typeof(sys))
    ntuple(i -> i <= length(names) ? names[i] : OutputName, length(names) + 1)
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys._box.unit_cell.matrix, 1)
    println(io, "PeriodicSystem1{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println()
    show(io_sub, mime, sys._cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.nbatches)
    print("\n    Type of output variable ($OutputName)): $(typeof(sys.output))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys._box.unit_cell.matrix, 1)
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
and implemented for `Number` types, and for `AbstractVecOrMat` containers of `Number`s or `SVector`s.

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

        Please implement a method `PeriodicSystems.copy_output(x::$(typeof(x))` defining
        an appropriate way to copy the required output variable. Many times just
        defining `output_copy(x::$(typeof(x))) = deepcopy(x)` is ok. 
    """
    )
end
copy_output(x::Number) = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T} = copy(x)

"""

```
reset_output!(x)
```

Function that defines how to reset (or zero) the `output` variable. For `Number`s it is 
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

The `reset_output!` function **must** return the output variable, being
it mutable or immutable. The user must guarantee that the operation takes place in-place,
for mutable output variables.  

"""
function reset_output!(x)
    error("""
        MethodError: no method matching `reset_output!($(typeof(x)))`

        Please add a method `PeriodicSystems.reset_output!(x::$(typeof(x))`, defining
        the appropriate way to reset (zero) the data of the output variables.

        The `reset_output!` methods **must** return the output variable to
        conform with the interface, even if the variable is mutable. For example:

        ```
        struct A x::Float64 end
        PeriodicSystems.reset_output!(v::Vector{A}) = fill!(v, A(0.0))
        ```
    """
    )
end
reset_output!(x::Number) = zero(x)
reset_output!(x::SVector) = zero(x)
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

The most commont `reducer` is the sum. For example, when computin energies, or forces,
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
reducer!(x,y) = +(x,y)
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
           positions = rand(SVector{3,Float64},1000), 
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
    sys._box = CellListMap.Box(unitcell, sys._box.cutoff)
    return sys
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
           positions = rand(SVector{3,Float64},1000), 
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
function update_cutoff!(sys::PeriodicSystem1{V,<:CellListMap.Box{UnitCellType}}, cutoff) where {V,UnitCellType}
    sys._box = CellListMap.Box(sys._box.unit_cell.matrix, cutoff; UnitCellType=UnitCellType)
    return sys
end
function update_cutoff!(sys::PeriodicSystem2{V,<:CellListMap.Box{UnitCellType}}, cutoff) where {V,UnitCellType}
    sys._box = CellListMap.Box(sys._box.unit_cell.matrix, cutoff; UnitCellType=UnitCellType)
    return sys
end

"""

```
UpdatePeriodicSystem!
```

$(INTERNAL)

Updates the cell lists for periodic systems.

"""
function UpdatePeriodicSystem!(sys::PeriodicSystem1)
    sys._cell_list = CellListMap.UpdateCellList!(
        sys.positions, 
        sys._box, 
        sys._cell_list, 
        sys._aux; 
        parallel=sys.parallel
    )
    return sys
end

function UpdatePeriodicSystem!(sys::PeriodicSystem2)
    sys._cell_list = CellListMap.UpdateCellList!(
        sys.xpositions, 
        sys.ypositions, 
        sys._box, 
        sys._cell_list, 
        sys._aux; 
        parallel=sys.parallel
    )
    return sys
end

"""

```
map_pairwise!(f::Function, system; show_progress = true)
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

# Example

In this example we compute the sum of `1/(1+d)` where `d` is the
distance between particles of a set, for `d < cutoff`. 

```julia-repl
julia> sys = PeriodicSystem(
           positions = rand(SVector{3,Float64},1000), 
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
    show_progress::Bool=false
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded)
    UpdatePeriodicSystem!(sys)
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
