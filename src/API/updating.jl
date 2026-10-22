"""
    resize_output!(sys::AbstractParticleSystem, n::Int)

Resizes the output array and the auxiliary output arrays used
for multithreading, if the number of particles of the system changed.

The function will error if `Base.resize!` is not defined for the
type of `system.output`. In this case, a `Base.resize!` method
must be implemented by the user.

!!! warning
    This function *must* be used whenever the output is an array and
    must be resiszed (for example for force arrays dependent on 
    variable number of particles) because it adjust the
    size of the copies of the output variable used for multi-threading.

"""
function resize_output!(sys::AbstractParticleSystem, n::Int)
    resize!(sys.output, n)
    _output_threaded = output_threaded(sys)
    for i in eachindex(_output_threaded)
        resize!(_output_threaded[i], n)
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
    _update!(sys, update_box(box(sys); unitcell))
    _updated_sys_coordinates!(sys)
    return sys
end
_updated_sys_coordinates!(sys::ParticleSystem1) = sys.xpositions.updated[] = true
function _updated_sys_coordinates!(sys::ParticleSystem2)
    sys.xpositions.updated[] = true
    sys.ypositions.updated[] = true
end

function _update_cutoff!(sys::ParticleSystem1, cutoff)
    _cutoff = isnothing(cutoff) ? sys.cutoff : cutoff
    if unitcelltype(sys) == NonPeriodicCell
        _update!(sys, Box(limits(sys.xpositions), _cutoff))
    else
        _update!(sys, update_box(box(sys); cutoff=_cutoff))
    end
    sys.xpositions.updated[] = true
    return sys
end
function _update_cutoff!(sys::ParticleSystem2, cutoff)
    _cutoff = isnothing(cutoff) ? sys.cutoff : cutoff
    if unitcelltype(sys) == NonPeriodicCell
        _update!(sys, Box(limits(sys.xpositions, sys.ypositions), _cutoff))
    else
        _update!(sys, update_box(box(sys); cutoff=_cutoff))
    end
    sys.xpositions.updated[] = true
    sys.ypositions.updated[] = true
    return sys
end

# Internal helpers: update a ParticleSystemPositions from any supported coordinate type,
# resizing if the number of particles changed. Each method copies element-wise into the
# existing storage without allocating a temporary Vector{SVector}. The copy is split
# in chunks run in parallel if `parallel == true` and the number of particles is large.

# Minimum number of particles for the copy to be run in parallel (below this, the
# cost of spawning the tasks is comparable to the cost of the copy)
const PARALLEL_COPY_MIN_PARTICLES = 2^15

# Calls `f(range)` for consecutive chunks of `1:n`, in parallel if `parallel == true`
function _foreach_chunk(f::F, n::Int, parallel::Bool) where {F<:Function}
    if !parallel || nthreads() == 1 || n < PARALLEL_COPY_MIN_PARTICLES
        f(1:n)
    else
        # (An explicit vector of tasks, rather than `@sync`, which waits on untyped tasks.)
        tasks = [@spawn f($range) for range in index_chunks(1:n; n = nthreads(), split = Consecutive())]
        foreach(wait, tasks)
    end
    return nothing
end

# Same element type: plain copy into psp.x
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractVector{SVector{N,T}}; parallel::Bool = false) where {N,T}
    n_new = length(new_x)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    i0 = firstindex(new_x) - 1
    _foreach_chunk(n_new, parallel) do range
        copyto!(psp.x, first(range), new_x, i0 + first(range), length(range))
    end
    return psp
end

# Vector of plain (non-SVector) vectors: construct SVectors element-wise
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractVector{<:AbstractVector}; parallel::Bool = false) where {N,T}
    n_new = length(new_x)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    i0 = firstindex(new_x) - 1
    _foreach_chunk(n_new, parallel) do range
        for j in range
            psp.x[j] = SVector{N,T}(new_x[i0 + j])
        end
    end
    return psp
end

# (D, N) matrix: columns are particles
function _update_positions!(psp::ParticleSystemPositions{N,T}, new_x::AbstractMatrix; parallel::Bool = false) where {N,T}
    n_new = size(new_x, 2)
    if n_new != length(psp)
        resize!(psp, n_new)
    end
    psp.updated[] = true
    i0 = firstindex(new_x, 2) - 1
    _foreach_chunk(n_new, parallel) do range
        for j in range
            psp.x[j] = SVector{N,T}(@view(new_x[:, i0 + j]))
        end
    end
    return psp
end

"""
    update!(
        sys::AbstractParticleSystem;
        xpositions = nothing,
        #or
        positions = nothing,
        ypositions = nothing,
        cutoff = nothing,
        unitcell = nothing,
        parallel = nothing,
        rebuild = false,
    )

Update one or more properties of `sys` in a single call. Only the keyword
arguments that are provided (i.e. not `nothing`) are updated.

- `xpositions` (or `positions` as a shortcut for single-set systems): new
  coordinates for the first (or only) set of particles. Accepts the same types
  as the `ParticleSystem` constructor: a `Vector{SVector}`, a vector of vectors,
  or an `(D, N)` matrix. The internal storage is resized automatically if the
  number of particles changes.

- `ypositions`: new coordinates for the second set of particles (only valid for
  two-set systems). Same accepted types as `xpositions`.

- `cutoff`: new cutoff distance.

- `unitcell`: new unit cell. Must be of the same cell type (orthorhombic or
  triclinic) as the original system. Manual updating of non-periodic systems is
  not allowed.

- `parallel`: whether to use multi-threading (`true` or `false`).

- `rebuild`: if `true`, the cell lists are rebuilt immediately, reflecting any
  pending changes (from this call or from prior direct mutations of
  `sys.xpositions`/`sys.ypositions`). By default (`false`), the rebuild is
  deferred to the next call of `pairwise!` or `neighborlist!`, as usual.

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

julia> update!(sys; rebuild=true); # force the cell lists to be rebuilt now

```
"""
function update!(
    sys::AbstractParticleSystem;
    positions = nothing,
    xpositions = nothing,
    ypositions = nothing,
    cutoff = nothing,
    unitcell = nothing,
    parallel = nothing,
    rebuild = false,
)
    if !isnothing(positions) && !isnothing(xpositions)
        throw(ArgumentError("Either `positions` OR `xpositions` must be provided, not both."))
    end
    xpositions = isnothing(positions) ? xpositions : positions
    if !isnothing(ypositions) && !(sys isa ParticleSystem2)
        throw(ArgumentError("ypositions can only be set for a two-set particle system"))
    end
    _parallel = isnothing(parallel) ? sys.parallel : parallel
    !isnothing(xpositions) && _update_positions!(sys.xpositions, xpositions; parallel = _parallel)
    !isnothing(ypositions) && _update_positions!(sys.ypositions, ypositions; parallel = _parallel)
    !isnothing(cutoff) && _update_cutoff!(sys, cutoff)
    !isnothing(unitcell) && _update_unitcell!(sys, unitcell)
    !isnothing(parallel) && setfield!(sys, :parallel, parallel)
    rebuild && UpdateParticleSystem!(sys)
    return sys
end

