#=
    get_dim(unitcell, x, [y])

Try to determine the dimension (2 or 3) of the constructed system from the data
provided, error if imposible. 

=#
function get_dim(unitcell, x, y=nothing)
    Ds = ( 
        !isnothing(unitcell) ?  size(unitcell, 1) : 0,
        _element_length(x),
        _element_length(y),
    )
    D = maximum(Ds)
    if D == 0
        throw(ArgumentError("""\n
            Could not infer dimension (2 or 3) from the unitcell or coordinate arrays. 
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    if !(D in (2,3))
        throw(ArgumentError("""\n
            Dimension must be 2 or 3.
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    if any(d -> d != 0 && d != D, Ds)
        throw(ArgumentError("""\n
            Incompatible dimensions between unitcell and/or coordinates. 
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    return D
end

# When the coordinates are vectors of static vectors the dimension is known
# from their type, so it is returned as a compile-time constant (the unitcell,
# if given as a dynamic array, is only checked for consistency at runtime).
function get_dim(unitcell, x::AbstractVector{<:SVector{N}}, y::Union{Nothing,AbstractVector{<:SVector{N}}}=nothing) where {N}
    if !(N in (2,3))
        throw(ArgumentError("Dimension must be 2 or 3. Got: x: $N"))
    end
    if !isnothing(unitcell) && size(unitcell, 1) != N
        throw(ArgumentError("""\n
            Incompatible dimensions between unitcell and/or coordinates. 
            Got: unitcell: $(size(unitcell, 1)); x: $N

        """))
    end
    return N
end
_element_length(::Nothing) = 0
_element_length(::AbstractVector{<:SVector{N}}) where {N} = N
_element_length(x::AbstractVector{<:AbstractVector}) = isempty(x) ? 0 : length(first(x))
_element_length(x::AbstractMatrix) = size(x,1)
_element_length(::AbstractVector) = throw(ArgumentError(" Coordinates must be vectors of vectors, or matrices."))

get_dim(sys::AbstractParticleSystem) = get_dim(sys.xpositions)
get_dim(::ParticleSystemPositions{N}) where {N} = N

#
# Auxiliary functions to control the exhibition of the progress meter
#
# `show_progress = nothing` (the default) is known at compile time, such that the
# progress meter code is not even compiled (required for `juliac --trim`).
_progress(::Nothing, n) = nothing
_progress(show_progress::Bool, n) = show_progress ? Progress(n, dt = 1) : nothing
_next!(::Nothing) = nothing
_next!(p) = next!(p)

#=
    _n_workqueue_chunks(_nbatches, n_work_items)

# Extended help

Number of small `Consecutive` chunks to split a range into for on-demand
work-queue scheduling — used both for the map computation (`batch`/
`_pairwise_parallel!` in `self.jl`/`cross.jl`, which splits the *cell* range)
and for cell list construction (`UpdateCellList!`'s Phase 1 in
`CellLists.jl`/`NonPeriodicCells.jl`, which splits the *particle* range):
enough that a fixed pool of `_nbatches` worker tasks can rebalance across
cores of unequal speed (e.g. performance/efficiency hybrid CPUs — see
Finding 2 in `PERFORMANCE_NOTES_pairwise_scaling.md`), without needlessly
increasing task spawn/scheduling overhead. Only affects scheduling
granularity, not the number of accumulator slots (still `_nbatches`), so
this does not change memory usage.

`n_work_items` is a proxy for the *actual amount of work*, used only to size
the cap below — it is deliberately **not** required to be the same count as
whatever range is being split into chunks. In particular, the map step
splits the cell range, but sizes its cap using particle count: cell count
alone would badly understate real work for a dense system (few, heavily
populated cells — e.g. a water-density system can have ~100 particles/cell),
capping oversubscription far below where it would actually help. Callers
that split what they size by (construction splits and sizes by particle
count) just pass the same count for both.

Capped so each chunk represents at least `WORKQUEUE_MIN_CHUNK_SIZE` work
items on average — for small problems, oversubscribing by a fixed
`_nbatches` multiple regardless of size produces chunks so small that
per-chunk overhead (the atomic claim, task loop) dominates the (trivial)
actual work in each chunk; this still gives every worker multiple chunks to
pull from once there's enough work to matter for load balancing.

The size-based cap is never allowed to drop the chunk count below
`_nbatches` itself — capping below that starves most workers of any chunk at
all (some workers get zero chunks and sit idle for the whole call), which is
worse than no oversubscription — see PERFORMANCE_NOTES_pairwise_scaling.md.

=#
const WORKQUEUE_OVERSUBSCRIPTION = 16
const WORKQUEUE_MIN_CHUNK_SIZE = 256
function _n_workqueue_chunks(_nbatches, n_work_items)
    max_chunks_by_size = max(_nbatches, n_work_items ÷ WORKQUEUE_MIN_CHUNK_SIZE)
    return clamp(WORKQUEUE_OVERSUBSCRIPTION * _nbatches, _nbatches, max_chunks_by_size)
end
