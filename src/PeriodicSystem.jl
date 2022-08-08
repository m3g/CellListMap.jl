"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that stores the data required for cell list computations, for systems
with only one set of coordinates (self-pairwise computations).

"""
mutable struct PeriodicSystem{V,B,C,O,A}
    positions::Vector{V}
    box::B
    cell_list::C
    output::O
    output_threaded::Vector{O}
    aux::A
    parallel::Bool
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that stores the data required for cell list computations, for systems
with two sets of coordinates (cross-set computations).

"""
mutable struct PeriodicSystemPair{V,B,C,O,A}
    positions1::Vector{V}
    positions2::Vector{V}
    box::B
    cell_list::C
    output::O
    output_threaded::Vector{O}
    aux::A
    parallel::Bool
end

"""

```
CellListMap.copy_output(x)
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
CellListMap.reset_output(x)
```

Function that how to reset (or zero) the `output` variable. Identical to `Base.copy(x)`
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
function reset_output!(x)
    error("""
        MethodError: no method matching `reset_output!($(typeof(x)))`

        Please add a method `CellListMap.reset_output!(x::$(typeof(x))`, defining
        the appropriate way to reset (zero) the data of the output variables.

        The `reset_output!` methods **must** return the output variable to
        conform with the interface, even if the variable is mutable. For example:

        ```
        struct A x::Float64 end
        function reset_output!(v::Vector{A})
            for i in eachindex(v)
                v[i] = A(0.0)
            end
            return v
        end
        ```
    """
    )
end
reset_output!(x::Number) = zero(x)
reset_output!(x::AbstractVector{<:Number}) = fill!(x, zero(eltype(x)))
reset_output!(x::AbstractVector{<:SVector}) = fill!(x, zeros(eltype(x)))
reset_output!(x::AbstractMatrix{<:Number}) = fill!(x, zero(eltype(x)))

function _reset_all_output!(output, output_threaded)
    output = reset_output!(output)
    for i in eachindex(output_threaded)
        output_threaded[i] = reset_output!(output_threaded[i])
    end
    return output
end

function PeriodicSystem( 
    positions::AbstractVecOrMat,
    sides::AbstractVecOrMat,
    cutoff::Number,
    output::Any;
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
    box = Box(sides, cutoff)
    cell_list = CellList(positions, box; parallel=parallel, nbatches=nbatches)
    aux = AuxThreaded(cell_list)
    output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(cell_list)]
    output = _reset_all_output!(output, output_threaded)
    sys = PeriodicSystem(positions, box, cell_list, output, output_threaded, aux, parallel)
    return sys
end

function PeriodicSystem( 
    positions1::AbstractVecOrMat,
    positions2::AbstractVecOrMat,
    sides::AbstractVecOrMat,
    cutoff::Number,
    output::Any;
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
    box = Box(sides, cutoff)
    cell_list = CellList(positions1, positions2, box; parallel=parallel, nbatches=nbatches)
    aux = AuxThreaded(cell_list)
    output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(cell_list)]
    output = _reset_all_output!(output, output_threaded)
    sys = PeriodicSystemPair(positions1, positions2, box, cell_list, output, output_threaded, aux, parallel)
    return sys
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystem)
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.box.unit_cell.matrix, 1)
    println(io, "PeriodicSystem of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys.box)
    println()
    show(io_sub, mime, sys.cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys.cell_list.nbatches)
    print("\n    Type of output variable: $(eltype(sys.output_threaded))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::PeriodicSystemPair)
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.box.unit_cell.matrix, 1)
    println(io, "PeriodicSystemPair of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys.box)
    println()
    show(io_sub, mime, sys.cell_list)
    println("\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys.cell_list.nbatches)
    print("\n    Type of output variable: $(eltype(sys.output_threaded))")
end

function set_positions(sys::PeriodicSystem{V}, coor::V, i::Int) where {V}
    sys.positions[i] = coor
    return sys
end

function set_positions(sys::PeriodicSystem{V}, positions::Vector{V}) where {V}
    set_positions.(sys.positions, positions)
    return sys
end

function set_positions1(sys::PeriodicSystemPair{V}, coor::V, i::Int) where {V}
    sys.positions1[i] = coor
    return sys
end

function set_positions1(sys::PeriodicSystemPair{V}, positions::Vector{V}) where {V}
    set_positions1.(sys.positions1, positions)
    return sys
end

function set_positions2(sys::PeriodicSystemPair{V}, coor::V, i::Int) where {V}
    sys.positions2[i] = coor
    return sys
end

function set_positions2(sys::PeriodicSystemPair{V}, positions::Vector{V}) where {V}
    set_positions2.(sys.positions2, positions)
    return sys
end


get_positions(sys::PeriodicSystem, i::Int) = sys.positions[i]
get_positions1(sys::PeriodicSystemPair, i::Int) = sys.positions1[i]
get_positions2(sys::PeriodicSystemPair, i::Int) = sys.positions2[i]

get_output(sys::PeriodicSystem, i::Int) = sys.output[i]
get_output(sys::PeriodicSystemPair, i::Int) = sys.output[i]

# If the positions, or box, where updated independently, just update the
# cell lists with the new data.
function UpdatePeriodicSystem!(sys::PeriodicSystem)
    sys.cell_list = UpdateCellList!(sys.positions, sys.box, sys.cell_list, sys.aux)
    return sys
end

function UpdatePeriodicSystem!(sys::PeriodicSystemPair)
    sys.cell_list = UpdateCellList!(sys.positions1, sys.positions2, sys.box, sys.cell_list, sys.aux)
    return sys
end

function map_pairwise!(
    f::F,
    sys::Union{PeriodicSystem,PeriodicSystemPair};
    show_progress::Bool=false
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys.output_threaded)
    UpdatePeriodicSystem!(sys)
    sys.output = map_pairwise!(
        f, sys.output, sys.box, sys.cell_list;
        output_threaded=sys.output_threaded,
        parallel=sys.parallel,
        show_progress=show_progress
    )
    return sys
end
