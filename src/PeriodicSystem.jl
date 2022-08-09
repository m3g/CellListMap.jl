module PeriodicSystems

import ..CellListMap
import ..CellListMap: INTERNAL
using DocStringExtensions
using StaticArrays

export PeriodicSystem
export map_parwise!, map_pairwise
export update_cutoff!
export update_unitcell!

"""

```
PeriodicSystem( 
    positions::AbstractVecOrMat,
    #or
    xpositions::AbstractVecOrMat,
    ypositions::AbstractVecOrMat,
    # and
    unitcell::AbstractVecOrMat,
    cutoff::Number,
    output::Any;
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
```

Function that sets up the `PeriodicSystem` type given a set of positions.

Particle positions can be provided as vectors of 2D or 3D vectors 
(preferentially static vectors from `StaticArrays`).

If the `positions` array is provided, a single set of coordinates is considered,
and the computation will be mapped for the `N(N-1)` pairs of this set. 

If the `xpositions` and `ypositions` arrays of coordinates are provided, the computation
will be mapped to the `N×M` pairs of particles, being `N` and `M` the number
of particles of each set of coordinates.

The unit cell (either a vector for `Orthorhombic` cells or a 
full unit cell matrix for `Triclinic` cells), the cutoff used for the
construction of the cell lists and the output variable of the calculations.

The `parallel` and `nbatches` flags control the parallelization scheme of
computations (see the user manual).

## Example

```julia-repl



```

"""
function PeriodicSystem(;
    positions::Union{Nothing,AbstractVecOrMat}=nothing,
    xpositions::Union{Nothing,AbstractVecOrMat}=nothing,
    ypositions::Union{Nothing,AbstractVecOrMat}=nothing,
    unitcell::AbstractVecOrMat,
    cutoff::Number,
    output::Any,
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
    
)
    if !isnothing(positions) && (isnothing(xpositions) && isnothing(ypositions))
        _box = Box(unitcell, cutoff)
        _cell_list = CellList(positions, _box; parallel=parallel, nbatches=nbatches)
        _aux = AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem1(positions, output, parallel, _box, _cell_list, _output_threaded, _aux)
    elseif isnothing(positions) && (!isnothing(xpositions) && !isnothing(ypositions))
        _box = Box(unitcell, cutoff)
        _cell_list = CellList(xpositions, ypositions, _box; parallel=parallel, nbatches=nbatches)
        _aux = AuxThreaded(cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = PeriodicSystem2(xpositions, ypositions, output, parallel, _box, _cell_list, _output_threaded, _aux)
    else
       throw(ArgumentError("""
           Either define `positions` **or** `xpositions` AND `ypositions`, to build
           systems for self- or cross-pair computations, respectively.  
       """)) 
    end
    return sys
end

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
auxliary function.

# Example

In this example we build the system and compute the sum of the squared distances
between particles within the cutoff. 

```julia-repl
julia> x = rand(SVector{3,Float64}, 100)
       unitcell = [1,1,1]
       cutoff = 0.1
       output = 0.0;

julia> sys = PeriodicSystem(positions = x, unitcell=[1,1,1], cutoff = 0.1, output = 0.0)
PeriodicSystem1 of dimension 3, composed of:
    Box{OrthorhombicCell, 3}
      unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
      cutoff = 0.1
      number of computing cells on each dimension = [12, 12, 12]
      computing cell sizes = [0.1, 0.1, 0.1] (lcell: 1)
      Total number of cells = 1728
    CellList{3, Float64}
      100 real particles.
      95 cells with real particles.
      161 particles in computing box, including images.
    Parallelization auxiliary data set for: 
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 8
    Type of output variable: Float64

julia> map_pairwise!((x,y,i,j,d2,output) -> output += d2, sys)
0.16331317695853984

```

"""
mutable struct PeriodicSystem1{V,B,C,O,A}
    positions::Vector{V}
    output::O
    parallel::Bool
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that stores the data required for cell list computations, for systems
with two sets of coordinates (cross-set computations).

"""
mutable struct PeriodicSystem2{V,B,C,O,A}
    xpositions::Vector{V}
    ypositions::Vector{V}
    output::O
    _parallel::Bool
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem1)
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys._box.unit_cell.matrix, 1)
    println(io, "PeriodicSystem1 of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println()
    show(io_sub, mime, sys._cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.nbatches)
    print("\n    Type of output variable: $(eltype(sys._output_threaded))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem2)
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys._box.unit_cell.matrix, 1)
    println(io, "PeriodicSystem2 of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println()
    show(io_sub, mime, sys._cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.target.nbatches)
    print("\n    Type of output variable: $(eltype(sys._output_threaded))")
end

"""

```
copy_output(x)
```

Function that defines how the `output` variable is copied. Identical to `Base.copy(x)`
and implemented 
for `Number` types, and for `AbstractVecOrMat` containers of `Number`s or `SVector`s.

Other custom output types must have their `copy_output` method implemented.

# Example

```julia
# Custom data type
struct A x::Int end
# Custom output type (array of A)
output = [ A(0) for _ in 1:100 ]
# How to copy an array of `A`
CellListMap.copy_output(v::Vector{A}) = [ x for x in v ]

# Alternativelly, in this case, one could have defined:
Base.copy(a::A) = a
CellListMap.copy_output(v::Vector{A}) = copy(v)
```

The user must guarantee that the copy is independent of the original array.
For many custom types it is possible to define 
`CellListMap.copy_output(v::Vector{T}) where {T<:CustomType} = deepcopy(v)`.

"""
function copy_output(x)
    error("""
        MethodError: no method matching `copy_output($(typeof(x)))`

        Please implement a method `CellListMap.copy_output(x::$(typeof(x))` defining
        an appropriate way to copy the required output variable. Many times just
        defining `output_copy(x::$(typeof(x))) = deepcopy(x)` is ok. 
    """
    )
end
copy_output(x::Number) = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T <:Union{Number,SVector}} = copy(x)

"""

```
reset_output!(x)
```

Function that how to reset (or zero) the `output` variable. For `Number`s it is 
implemented as `zero(x)`, and for `AbstractVecOrMat` containers of `Number`s or `SVector`s
it is implemented as `fill!(x, zero(eltype(x))`.

Other custom output types must have their `reset_output!` method implemented.

# Example

```julia
# Custom data type
struct A x::Int end
# Custom output type (array of A)
output = [ A(0) for _ in 1:100 ]
# How to reset an array with elements of type `A`
CellListMap.reset_output!(v::Vector{A}) = fill!(v, A(0))
```

The `reset_output!` function **must** return the output variable, being
it mutable or immutable. The user must guarantee that the operation takes place in-place,
for mutable output variables.  

"""
function reset_output!(x)
    error("""
        MethodError: no method matching `reset_output!($(typeof(x)))`

        Please add a method `CellListMap.reset_output!(x::$(typeof(x))`, defining
        the appropriate way to reset (zero) the data of the output variables.

        The `reset_output!` methods **must** return the output variable to
        conform with the interface, even if the variable is mutable. For example:

        ```
        struct A x::Float64 end
        CellListMap.reset_output!(v::Vector{A}) = fill!(v, A(0.0))
        ```
    """
    )
end
reset_output!(x::Number) = zero(x)
reset_output!(x::AbstractVecOrMat{T}) where {T <:Union{Number,SVector}} = fill!(x, zero(T))

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


function reduce_output!(output, output_threaded)
    return CellListMap.reduce(output, output_threaded)
end


#
# Function used to update the properties of the systems
#
"""

"""
update_unitcell!(sys, unitcell) = sys._box = Box(unitcell, sys._box.cutoff)

function update_cutoff!(sys::PeriodicSystem1{V,<:Box{UnitCellType}}, cutoff) where {V,UnitCellType} 
    sys._box = Box(sys._box.unit_cell.matrix, cutoff; UnitCellType = UnitCellType)
end
function update_cutoff!(sys::PeriodicSystem2{V,<:Box{UnitCellType}}, cutoff) where {V,UnitCellType} 
    sys._box = Box(sys._box.unit_cell.matrix, cutoff; UnitCellType = UnitCellType)
end

"""

```
UpdatePeriodicSystem!
```

$(INTERNAL)

Updates the cell lists for periodic systems.

"""
function UpdatePeriodicSystem!(sys::PeriodicSystem1)
    sys._cell_list = UpdateCellList!(sys.positions, sys._box, sys._cell_list, sys._aux)
    return sys
end

function UpdatePeriodicSystem!(sys::PeriodicSystem2)
    sys._cell_list = UpdateCellList!(sys.xpositions, sys.ypositions, sys._box, sys._cell_list, sys._aux)
    return sys
end

"""

"""
function map_pairwise!(
    f::F,
    sys::Union{PeriodicSystem1,PeriodicSystem2};
    show_progress::Bool=false
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded)
    UpdatePeriodicSystem!(sys)
    sys.output = map_pairwise!(
        f, sys.output, sys._box, sys._cell_list;
        output_threaded=sys._output_threaded,
        parallel=sys.parallel,
        show_progress=show_progress
    )
    return sys.output
end

end # module PeriodicSystems
