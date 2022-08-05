
mutable struct PeriodicSystem{V,B,C,O,A}
    coordinates::Vector{V}
    box::B
    cell_list::C
    output::O
    output_threaded::Vector{O}
    aux::A
    parallel::Bool
end

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
copy_output(x::Vector{<:Number}) = copy(x)
copy_output(x::Vector{<:SVector}) = copy(x)
copy_output(x::AbstractMatrix{<:Number}) = copy(x)

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

function PeriodicSystem(;
    coordinates=nothing,
    sides=nothing,
    cutoff=nothing,
    output=nothing,
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
    sys = PeriodicSystem(
        coordinates, sides, cutoff, output;
        parallel=parallel, nbatches=nbatches
    )
    return sys
end

function PeriodicSystem(
    coordinates::AbstractVecOrMat,
    sides::AbstractVector,
    cutoff::Real,
    output;
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
)
    box = Box(sides, cutoff)
    cell_list = CellList(coordinates, box; parallel=parallel, nbatches=nbatches)
    aux = AuxThreaded(cell_list)
    output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(cell_list)]
    output = _reset_all_output!(output, output_threaded)
    sys = PeriodicSystem(coordinates, box, cell_list, output, output_threaded, aux, parallel)
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

function set_coordinates!(sys::PeriodicSystem{V}, coor::V, i::Int) where {V}
    sys.coordinates[i] = coor
    return sys
end

function set_coordinates!(sys::PeriodicSystem{V}, coordinates::Vector{V}) where {V}
    for i in eachindex(sys.coordinates, coordinates)
        @inbounds sys.coordinates[i] = coor[i]
    end
    return sys
end

get_coordinates(sys::PeriodicSystem, i::Int) = sys.coordinates[i]
get_output(sys::PeriodicSystem, i::Int) = sys.output[i]

# If the coordinates, or box, where updated independently, just update the
# cell lists with the new data.
function UpdatePeriodicSystem!(sys::PeriodicSystem)
    sys.cell_list = UpdateCellList!(sys.coordinates, sys.box, sys.cell_list, sys.aux)
    return sys
end

function UpdatePeriodicSystem!(
    sys::PeriodicSystem, x;
    sides=nothing
) where {N}
    # update box if necessary
    if !isnothing(sides)
        sys.box = Box(sides, sys.box.cutoff)
    end
    sys.cell_list = UpdateCellList!(x, sys.box, sys.cell_list, sys.aux)
    return sys
end

function map_pairwise!(
    f::F,
    sys::PeriodicSystem;
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
