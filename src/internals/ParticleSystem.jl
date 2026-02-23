# Types of variables that have support for multi-threading without having
# to explicit add methods to copy_output, reset_output!, and reducer functions.
const SupportedTypes = Union{Number, SVector, FieldVector}

# Supported types for coordinates
const SupportedCoordinatesTypes = Union{Nothing, AbstractVector{<:AbstractVector}, AbstractMatrix}

#=
    unitcelltype(sys::AbstractParticleSystem)

Returns the type of a unitcell from the `ParticleSystem` structure.

=#
unitcelltype(sys::AbstractParticleSystem) = unitcelltype(sys._box)

ParticleSystem1{OutputName}(v::V, o::O, b::B, c::C, vo::AbstractVector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem1{OutputName, V, O, B, C, A, VC}(v, o, b, c, vo, a, p, vc)
getproperty(sys::ParticleSystem1, ::Val{:positions}) = getfield(sys, :xpositions)
getproperty(sys::ParticleSystem2, ::Val{:positions}) = getfield(sys, :xpositions)

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem1{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(sys._cell_list, :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem2{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(sys._cell_list, :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

#=
    _reset_all_output!(output, output_threaded)

Function that resets the output variable and the threaded copies of it.

=#
function _reset_all_output!(output, output_threaded; reset::Bool)
    if reset
        output = reset_output!(output)
    end
    for i in eachindex(output_threaded)
        output_threaded[i] = reset_output!(output_threaded[i])
    end
    return output
end

function reduce_output!(output::T, output_threaded::Vector{T}) where {T}
    for ibatch in eachindex(output_threaded)
        output = reducer(output, output_threaded[ibatch])
    end
    return output
end
function reduce_output!(
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

#=
    UpdateParticleSystem!

Updates the cell lists for periodic systems.

=#
function UpdateParticleSystem!(sys::ParticleSystem1)
    if sys.xpositions.updated[]
        if unitcelltype(sys) == NonPeriodicCell
            sys._box = Box(limits(sys.xpositions), sys.cutoff)
        end
        n_particles_changed = length(sys.xpositions) != sys._cell_list.n_real_particles
        sys._cell_list = UpdateCellList!(
            sys.xpositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel = sys.parallel,
            validate_coordinates = sys.validate_coordinates,
        )
        if n_particles_changed
            _old_nbatches = get_nbatches(sys)
            sys._cell_list = update_number_of_batches!(sys._cell_list; parallel = sys.parallel)
            _new_nbatches = get_nbatches(sys)
            if _old_nbatches != _new_nbatches
                sys._aux = AuxThreaded(sys._cell_list)
                sys._output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
            end
        end
        sys.xpositions.updated[] = false
    end
    return sys
end

# Check if the new limits are compatible with the current box: the new particle
# bounding box fits within the old one, so the box does not need to be updated.
function _limits_fit_in_box(new_limits::Limits{N}, box::Box) where {N}
    old_origin = box.origin
    # Recover the old limits extent by undoing the cutoff padding from the box sides.
    old_upper = old_origin .+ SVector(ntuple(i -> box.input_unit_cell.matrix[i,i], Val(N))) .- (210 * box.cutoff / 100)
    new_upper = new_limits.origin .+ new_limits.limits
    return all(new_limits.origin .>= old_origin) && all(new_upper .<= old_upper)
end

function UpdateParticleSystem!(sys::ParticleSystem2)
    if sys.xpositions.updated[] || sys.ypositions.updated[]
        if unitcelltype(sys) == NonPeriodicCell
            new_limits = limits(sys.xpositions, sys.ypositions)
            if !_limits_fit_in_box(new_limits, sys._box)
                sys._box = Box(new_limits, sys._box.cutoff)
                # Both cell lists must be rebuilt when the box changes,
                # because the origin and sides changed.
                sys.xpositions.updated[] = true
                sys.ypositions.updated[] = true
            end
        end
        n_particles_changed = (length(sys.xpositions) != sys._cell_list.ref_list.n_real_particles) ||
            (length(sys.ypositions) != sys._cell_list.target_list.n_real_particles)
        sys._cell_list = UpdateCellList!(
            sys.xpositions,
            sys.ypositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel = sys.parallel,
            validate_coordinates = sys.validate_coordinates,
        )
        if n_particles_changed
            _old_nbatches = get_nbatches(sys)
            sys._cell_list = update_number_of_batches!(sys._cell_list; parallel = sys.parallel)
            _new_nbatches = get_nbatches(sys)
            if _old_nbatches != _new_nbatches
                sys._aux = AuxThreaded(sys._cell_list)
                sys._output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
            end
        end
        sys.xpositions.updated[] = false
        sys.ypositions.updated[] = false
    end
    return sys
end

# Return the number of batches for ParticleSystems
"""
    get_nbatches(sys::AbstractParticleSystem)

Returns the number of batches for parallel computations of cell list construction
and function mapping.

"""
get_nbatches(sys::ParticleSystem1) = nbatches(sys._cell_list)
get_nbatches(sys::ParticleSystem2) = nbatches(sys._cell_list.ref_list)
