
# Types of variables that have support for multi-threading without having
# to explicit add methods to copy_output, reset_output!, and reducer functions.
const SupportedTypes = Union{Number, SVector, FieldVector}

# Supported types for coordinates
const SupportedCoordinatesTypes = Union{Nothing, AbstractVector{<:AbstractVector}, AbstractMatrix}

#
# Internal getters
#
box(sys::AbstractParticleSystem) = sys.private.box
cutoff(sys::AbstractParticleSystem) = sys.private.box.cutoff
unitcell(sys::AbstractParticleSystem) = sys.private.box.input_unit_cell.matrix
unitcelltype(sys::AbstractParticleSystem) = unitcelltype(sys.private.box)
dim(sys::AbstractParticleSystem) = size(unitcell(sys), 1)
output_type(sys) = typeof(sys.output)

aux_threaded(sys) = sys.private.aux
output_threaded(sys) = sys.private.output_threaded

celllist(sys::AbstractParticleSystem) = sys.private.cell_list
celllist(sys::ParticleSystem2, s::Symbol) = celllist(sys, Val(s))
celllist(sys, ::Val{:ref}) = sys.private.cell_list.ref_list
celllist(sys, ::Val{:target}) = sys.private.cell_list.target_list

n_real_particles(sys, s::Symbol) = n_real_particles(sys, Val(s)) 
n_real_particles(sys, ::Val{:ref}) = sys.private.cell_list.ref_list.n_real_particles
n_real_particles(sys, ::Val{:target}) = sys.private.cell_list.target_list.n_real_particles

updated(sys, s::Symbol) = updated(sys, Val(s))
updated(sys, ::Val{:x}) = sys.xpositions.updated[]
updated(sys, ::Val{:y}) = sys.ypositions.updated[]

# Constructor for ParticleSystem1
ParticleSystem1{OutputName}(vx::V, o::O, b::B, c::C, vo::AbstractVector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem1{OutputName, V, O, VC, PrivateParticleSystemData{B, C, O, A}}(vx, o, p, vc, PrivateParticleSystemData(b, c, vo, a))

# Constructor for ParticleSystem2
ParticleSystem2{OutputName}(vx::V, vy::V, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem2{OutputName, V, O, VC, PrivateParticleSystemData{B, C, O, A}}(vx, vy, o, p, vc, PrivateParticleSystemData(b, c, vo, a))

# Alias "positions" to "xpositions"
getproperty(sys::ParticleSystem1, ::Val{:positions}) = getfield(sys, :xpositions)
getproperty(sys::ParticleSystem2, ::Val{:positions}) = getfield(sys, :xpositions)

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    println(io, "ParticleSystem1{$OutputName} of dimension $(dim(sys)), composed of:")
    show(IOContext(io, :indent => indent + 4), mime, box(sys))
    println(io)
    show(io_sub, mime, celllist(sys))
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(celllist(sys), :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(output_type(sys))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    println(io, "ParticleSystem2{$OutputName} of dimension $(dim(sys)), composed of:")
    show(IOContext(io, :indent => indent + 4), mime, box(sys))
    println(io)
    show(io_sub, mime, celllist(sys))
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(celllist(sys), :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(output_type(sys))")
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

#
# Internal updaters
#
_update!(sys::AbstractParticleSystem, box::Box) = sys.private.box = box
function _update!(sys::ParticleSystem1, ::Type{CellList}) 
    n_particles_changed = length(sys.xpositions) != sys.private.cell_list.n_real_particles
    UpdateCellList!(
        sys.xpositions,
        box(sys),
        celllist(sys),
        aux_threaded(sys);
        parallel = sys.parallel,
        validate_coordinates = sys.validate_coordinates,
    )
    if n_particles_changed
        _old_nbatches = get_nbatches(sys)
        sys.private.cell_list = update_number_of_batches!(celllist(sys); parallel = sys.parallel)
        _new_nbatches = get_nbatches(sys)
        if _old_nbatches != _new_nbatches
            sys.private.aux = _create_aux(box(sys), celllist(sys))
            sys.private.output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
        end
    end
    sys.xpositions.updated[] = false
    return sys
end

#=
    UpdateParticleSystem!

Updates the cell lists for ParticleSystem1

=#
function UpdateParticleSystem!(sys::ParticleSystem1)
    if updated(sys, :x)
        unitcelltype(sys) == NonPeriodicCell && _update!(sys, Box(limits(sys.xpositions), cutoff(sys)))
        _update!(sys, CellList)
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

function _update!(sys::ParticleSystem2, ::Type{CellList})
    n_particles_changed = 
        (length(sys.xpositions) != n_real_particles(sys, :ref)) ||
        (length(sys.ypositions) != n_real_particles(sys, :target))
    sys.private.cell_list = UpdateCellList!(
        sys.xpositions,
        sys.ypositions,
        box(sys),
        celllist(sys),
        aux_threaded(sys);
        parallel = sys.parallel,
        validate_coordinates = sys.validate_coordinates,
    )
    if n_particles_changed
        _old_nbatches = get_nbatches(sys)
        sys.private.cell_list = update_number_of_batches!(celllist(sys); parallel = sys.parallel)
        _new_nbatches = get_nbatches(sys)
        if _old_nbatches != _new_nbatches
            sys.private.aux = _create_aux(box(sys), celllist(sys))
            sys.private.output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
        end
    end
    sys.xpositions.updated[] = false
    sys.ypositions.updated[] = false
    return sys
end

#=
    UpdateParticleSystem!

Updates the cell lists for ParticleSystem2

=#
function UpdateParticleSystem!(sys::ParticleSystem2)
    if updated(sys, :x) || updated(sys, :y) 
        if unitcelltype(sys) == NonPeriodicCell
            new_limits = limits(sys.xpositions, sys.ypositions)
            if !_limits_fit_in_box(new_limits, box(sys))
                _update!(sys, Box(new_limits, cutoff(sys)))
                # Both cell lists must be rebuilt when the box changes,
                # because the origin and sides changed.
                sys.xpositions.updated[] = true
                sys.ypositions.updated[] = true
            end
        end
        _update!(sys, CellList)
    end
    return sys
end

# Return the number of batches for ParticleSystems
"""
    get_nbatches(sys::AbstractParticleSystem)

Returns the number of batches for parallel computations of cell list construction
and function mapping.

"""
get_nbatches(sys::ParticleSystem1) = nbatches(celllist(sys))
get_nbatches(sys::ParticleSystem2) = nbatches(celllist(sys, :ref))
