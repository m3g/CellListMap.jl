"""
    resize_output!(sys::AbstractParticleSystem, n::Int)

Resizes the output array and the auxiliary output arrays used
for multithreading, if the number of particles of the system changed.

The function will error if `Base.resize!` is not defined for the
type of `system.output`. In this case, a `Base.resize!` method
must be implemented by the user.

!!! warn
    This function *must* be used whenever the output is dependent on
    the number of particles, and that changes, because it adjust the
    size of the copies of the output variable used for multi-threading.

"""
function resize_output!(sys::AbstractParticleSystem, n::Int)
    resize!(sys.output, n)
    for i in eachindex(sys._output_threaded)
        resize!(sys._output_threaded[i], n)
    end
    return sys
end

#
# Internal functions used to update the properties of the systems
#
function _update_unitcell!(sys, unitcell)
    isnothing(unitcell) && return sys
    if unitcelltype(sys) == NonPeriodicCell
        throw(
            ArgumentError(
                """\n
                    Manual updating of the unit cell of non-periodic systems is not allowed.

                """
            )
        )
    end
    sys._box = update_box(sys._box; unitcell)
    _update_sys_box!(sys)
    return sys
end
_update_sys_box!(sys::ParticleSystem1) = sys.xpositions.updated[] = true
function _update_sys_box!(sys::ParticleSystem2)
    sys.xpositions.updated[] = true
    sys.ypositions.updated[] = true
end

function _update_cutoff!(sys::ParticleSystem1, cutoff)
    _cutoff = isnothing(cutoff) ? sys.cutoff : cutoff
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions), _cutoff)
    else
        sys._box = update_box(sys._box; cutoff=_cutoff)
    end
    sys.xpositions.updated[] = true
    return sys
end
function _update_cutoff!(sys::ParticleSystem2, cutoff)
    _cutoff = isnothing(cutoff) ? sys.cutoff : cutoff
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions, sys.ypositions), _cutoff)
    else
        sys._box = update_box(sys._box; cutoff=_cutoff)
    end
    sys.xpositions.updated[] = true
    sys.ypositions.updated[] = true
    return sys
end

# Internal helpers: update a ParticleSystemPositions from any supported coordinate type,
# resizing if the number of particles changed. Each method copies element-wise into the
# existing storage without allocating a temporary Vector{SVector}.

# Same element type: plain copy into psp.x
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractVector{SVector{N,T}}) where {N,T}
    n_new = length(new_x)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    copyto!(psp.x, new_x)
    return psp
end

# Vector of plain (non-SVector) vectors: construct SVectors element-wise
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractVector{<:AbstractVector}) where {N,T}
    n_new = length(new_x)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    for (j, x_i) in enumerate(new_x)
        psp.x[j] = SVector{N,T}(x_i)
    end
    return psp
end

# (D, N) matrix: columns are particles
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractMatrix) where {N,T}
    n_new = size(new_x, 2)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    for (j, i) in enumerate(axes(new_x, 2))
        psp.x[j] = SVector{N,T}(@view(new_x[:, i]))
    end
    return psp
end

"""
    update!(
        sys::AbstractParticleSystem;
        xpositions = nothing,
        ypositions = nothing,
        cutoff = nothing,
        unitcell = nothing,
        parallel = nothing,
    )

Update one or more properties of `sys` in a single call. Only the keyword
arguments that are provided (i.e. not `nothing`) are updated.

- `xpositions`: new coordinates for the first (or only) set of particles.
  Accepts the same types as the `ParticleSystem` constructor: a
  `Vector{SVector}`, a vector of vectors, or an `(D, N)` matrix.
  The internal storage is resized automatically if the number of particles changes.

- `ypositions`: new coordinates for the second set of particles (only valid for
  two-set systems). Same accepted types as `xpositions`.

- `cutoff`: new cutoff distance.

- `unitcell`: new unit cell. Must be of the same cell type (orthorhombic or
  triclinic) as the original system. Manual updating of non-periodic systems is
  not allowed.

- `parallel`: whether to use multi-threading (`true` or `false`).

# Example

```jldoctest ;filter = r"( +Parallelization.*)" => ""
julia> using CellListMap, StaticArrays

julia> sys = ParticleSystem(
           xpositions = rand(SVector{3,Float64}, 100),
           unitcell = [1,1,1],
           cutoff = 0.1,
           output = 0.0,
       );

julia> new_x = rand(SVector{3,Float64}, 100);

julia> update!(sys; xpositions=new_x, cutoff=0.2, unitcell=[2,2,2], parallel=false);

```
"""
function update!(
    sys::AbstractParticleSystem;
    xpositions = nothing,
    ypositions = nothing,
    cutoff = nothing,
    unitcell = nothing,
    parallel = nothing,
)
    if !isnothing(ypositions) && !(sys isa ParticleSystem2)
        throw(ArgumentError("ypositions can only be set for a two-set particle system"))
    end
    !isnothing(xpositions) && _update_positions!(sys.xpositions, xpositions)
    !isnothing(ypositions) && _update_positions!(sys.ypositions, ypositions)
    !isnothing(cutoff) && _update_cutoff!(sys, cutoff)
    !isnothing(unitcell) && _update_unitcell!(sys, unitcell)
    !isnothing(parallel) && setfield!(sys, :parallel, parallel)
    return sys
end

#
# Deprecated public wrappers
#
"""
    update_unitcell!(system, unitcell)

!!! warning "Deprecated"
    `update_unitcell!` is deprecated. Use `update!(sys; unitcell=unitcell)` instead.

"""
function update_unitcell!(sys, unitcell)
    @warn """\n
        `update_unitcell!` is deprecated. Use `update!(sys; unitcell=unitcell)` instead.
    """ maxlog=1
    _update_unitcell!(sys, unitcell)
end

"""
    update_cutoff!(system, cutoff)

!!! warning "Deprecated"
    `update_cutoff!` is deprecated. Use `update!(sys; cutoff=cutoff)` instead.

"""
function update_cutoff!(sys, cutoff)
    @warn """\n
        `update_cutoff!` is deprecated. Use `update!(sys; cutoff=cutoff)` instead.
    """ maxlog=1
    _update_cutoff!(sys, cutoff)
end
