#
# Dedicated computation path for NonPeriodicCell systems.
#
# Key differences from the periodic path:
#   - No coordinate wrapping (positions stored as-is)
#   - No phantom/image particle creation (replicate_particle! not called)
#   - Cell list construction uses thread-local histograms instead of atomic CAS loops
#   - Inner loop and vicinal cell interactions carry no real/image particle checks
#

# Auxiliary structure for non-periodic cell list construction.
# thread_counts[ibatch] is a per-thread histogram over cell linear indices;
# it is repurposed as per-thread write offsets during the scatter phase.
mutable struct AuxNonPeriodic{N,T}
    idxs::Vector{UnitRange{Int}}
    thread_counts::Vector{Vector{Int}}
    total_np::Vector{Int}
end

mutable struct AuxNonPeriodicPair{N,T}
    ref_list::AuxNonPeriodic{N,T}
    target_list::AuxNonPeriodic{N,T}
end

function Base.show(io::IO, ::MIME"text/plain", aux::AuxNonPeriodic)
    _println(io, typeof(aux))
    return _print(io, " Auxiliary arrays for nbatches = ", length(aux.idxs))
end

function AuxNonPeriodic(cl::CellList{N,T}) where {N,T}
    _nbatches = nbatches(cl, :build)
    ncells = cl.number_of_cells
    idxs = Vector{UnitRange{Int}}(undef, _nbatches)
    thread_counts = [zeros(Int, ncells) for _ in 1:_nbatches]
    total_np = zeros(Int, ncells)
    set_idxs!(idxs, cl.n_real_particles, _nbatches)
    return AuxNonPeriodic{N,T}(idxs, thread_counts, total_np)
end

function AuxNonPeriodicPair(cl_pair::CellListPair{N,T}) where {N,T}
    return AuxNonPeriodicPair{N,T}(
        AuxNonPeriodic(cl_pair.ref_list),
        AuxNonPeriodic(cl_pair.target_list),
    )
end

#
# Dispatch helper: creates the appropriate auxiliary structure for the box type.
# Used in ParticleSystem constructors and update methods so that NonPeriodicCell
# systems get AuxNonPeriodic while periodic systems get the existing AuxThreaded.
#
_create_aux(::Box{NonPeriodicCell}, cl::CellList) = AuxNonPeriodic(cl)
_create_aux(::Box{NonPeriodicCell}, cl::CellListPair) = AuxNonPeriodicPair(cl)

#
# add_particles! — no wrapping, no image replication.
#
function add_particles!(x, box::Box{NonPeriodicCell}, ishift, cl::CellList{N,T}) where {N,T}
    for ip in eachindex(x)
        xp = x[ip]
        p = SVector{N,T}(ntuple(i -> xp[i], Val(N)))
        add_particle_to_celllist!(ishift + ip, p, box, cl)
    end
    return cl
end

# Resize projected_particles to hold the largest cell; called after every cell list build.
function _update_projected_particles!(cl::CellList)
    maxnp = 0
    for i in 1:cl.n_cells_with_particles
        maxnp = max(maxnp, cl.cells[i].n_particles)
    end
    for i in eachindex(cl.projected_particles)
        if maxnp > length(cl.projected_particles[i])
            resize!(cl.projected_particles[i], maxnp)
        end
    end
end

#
# UpdateCellList! for a single CellList.
#
# The Nothing overload resolves the ambiguity that would otherwise arise between
# the generic CellLists.jl method (aux::Union{Nothing,AuxThreaded}) and the
# AuxNonPeriodic method below: both unions contain Nothing, making them
# incomparable to Julia's dispatch even though NonPeriodicCell is more specific.
#
function UpdateCellList!(
    x::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl::CellList{N,T},
    ::Nothing;
    parallel::Bool=true,
    validate_coordinates::F=_validate_coordinates,
) where {N,T,F<:Function}
    isnothing(validate_coordinates) || validate_coordinates(x)
    reset!(cl, box, length(x))
    add_particles!(x, box, 0, cl)
    _update_projected_particles!(cl)
    x.updated[] = false
    return cl
end

#
# Algorithm (parallel path):
#   Phase 1 — each thread builds a local histogram over cell linear indices (no atomics).
#   Phase 2a — serial reduction of thread histograms into total_np.
#   Phase 2b — serial cell initialisation using total_np.
#   Phase 3 — serial computation of per-thread write offsets (thread_counts repurposed).
#   Phase 4 — parallel scatter: each thread writes its particles using pre-computed offsets.
#
function UpdateCellList!(
    x::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl::CellList{N,T},
    aux::AuxNonPeriodic{N,T};
    parallel::Bool=true,
    validate_coordinates::F=_validate_coordinates,
) where {N,T,F<:Function}

    isnothing(validate_coordinates) || validate_coordinates(x)

    _nbatches = nbatches(cl, :build)
    if !parallel || _nbatches == 1
        reset!(cl, box, length(x))
        add_particles!(x, box, 0, cl)
    else

        reset!(cl, box, 0)
        cl.n_real_particles = length(x)
        set_idxs!(aux.idxs, length(x), _nbatches)

        number_of_cells = cl.number_of_cells
        for ibatch in eachindex(aux.thread_counts)
            if length(aux.thread_counts[ibatch]) < number_of_cells
                resize!(aux.thread_counts[ibatch], number_of_cells)
            end
            fill!(view(aux.thread_counts[ibatch], 1:number_of_cells), 0)
        end
        if length(aux.total_np) < number_of_cells
            resize!(aux.total_np, number_of_cells)
        end
        fill!(view(aux.total_np, 1:number_of_cells), 0)

        # Phase 1: thread-local histograms — no atomic operations
        @sync for ibatch in eachindex(aux.idxs)
            @spawn begin
                prange = aux.idxs[ibatch]
                isempty(prange) && return
                counts = aux.thread_counts[ibatch]
                for ip in prange
                    xp = x[ip]
                    p = SVector{N,T}(ntuple(i -> xp[i], Val(N)))
                    ci = real_particle_border_case(particle_cell(p, box), box)
                    li = cell_linear_index(box.nc, ci)
                    counts[li] += 1
                end
            end
        end

        # Phase 2a: serial reduction of histograms
        for ibatch in eachindex(aux.thread_counts)
            counts = aux.thread_counts[ibatch]
            for li in 1:number_of_cells
                aux.total_np[li] += counts[li]
            end
        end

        # Phase 2b: serial cell initialisation
        for li in 1:number_of_cells
            np = aux.total_np[li]
            np == 0 && continue
            cl.n_cells_with_particles += 1
            cl.n_particles += np
            cl.cell_indices[li] = cl.n_cells_with_particles
            cell_index = cl.n_cells_with_particles
            cartesian_idx = cell_cartesian_indices(box.nc, li)
            center_c = cell_center(cartesian_idx, box)
            if cell_index > length(cl.cells)
                push!(
                    cl.cells, Cell{N,T}(
                        linear_index=li,
                        cartesian_index=cartesian_idx,
                        center=center_c,
                        contains_real=true,
                        n_particles=np,
                        particles=Vector{ParticleWithIndex{N,T}}(undef, np),
                    )
                )
            else
                cell = cl.cells[cell_index]
                @set! cell.linear_index = li
                @set! cell.cartesian_index = cartesian_idx
                @set! cell.center = center_c
                @set! cell.contains_real = true
                @set! cell.n_particles = np
                cl.cells[cell_index] = cell
                if np > length(cl.cells[cell_index].particles)
                    resize!(cl.cells[cell_index].particles, np)
                end
            end
            # All cells in a non-periodic system contain only real particles
            cl.n_cells_with_real_particles += 1
            if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                push!(cl.cell_indices_real, cell_index)
            else
                cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
            end
        end

        # Phase 3: compute per-thread write offsets (serial).
        # thread_counts[ibatch][li] is repurposed: original count → starting write offset
        # for thread ibatch at cell li.  total_np is reused as a running offset counter.
        fill!(view(aux.total_np, 1:number_of_cells), 0)
        for ibatch in eachindex(aux.thread_counts)
            counts = aux.thread_counts[ibatch]
            for li in 1:number_of_cells
                np = counts[li]
                np == 0 && continue
                counts[li] = aux.total_np[li]   # starting offset for this batch
                aux.total_np[li] += np
            end
        end

        # Phase 4: parallel scatter — write ranges are disjoint by construction
        @sync for ibatch in eachindex(aux.idxs)
            @spawn begin
                prange = aux.idxs[ibatch]
                isempty(prange) && return
                offsets = aux.thread_counts[ibatch]
                for ip in prange
                    xp = x[ip]
                    p = SVector{N,T}(ntuple(i -> xp[i], Val(N)))
                    ci = real_particle_border_case(particle_cell(p, box), box)
                    li = cell_linear_index(box.nc, ci)
                    target_cell_idx = cl.cell_indices[li]
                    offset = offsets[li]
                    offsets[li] += 1
                    cl.cells[target_cell_idx].particles[offset+1] = ParticleWithIndex(ip, true, p)
                end
            end
        end

    end # if/else parallel

    _update_projected_particles!(cl)
    x.updated[] = false
    return cl
end

# No-aux version: allocates AuxNonPeriodic and delegates.
function UpdateCellList!(
    x::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl::CellList;
    parallel::Bool=true,
    kargs...
)
    if parallel
        aux = AuxNonPeriodic(cl)
        cl = UpdateCellList!(x, box, cl, aux; parallel, kargs...)
    else
        cl = UpdateCellList!(x, box, cl, nothing; parallel, kargs...)
    end
    cl = update_number_of_batches!(cl; parallel)
    return cl
end

#
# UpdateCellList! for a CellListPair — delegates each list to the single-list method.
#
# Nothing overload resolves the same ambiguity as for the single-list version.
#
function UpdateCellList!(
    x::ParticleSystemPositions,
    y::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl_pair::CellListPair{N,T},
    ::Nothing;
    parallel::Bool=true,
    validate_coordinates::F=_validate_coordinates,
) where {N,T,F<:Function}
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    ref_list = cl_pair.ref_list
    target_list = cl_pair.target_list
    if x.updated[]
        UpdateCellList!(x, box, ref_list, nothing; validate_coordinates)
        x.updated[] = false
    end
    if y.updated[]
        UpdateCellList!(y, box, target_list, nothing; validate_coordinates)
        y.updated[] = false
    end
    return CellListPair{N,T}(ref_list, target_list)
end

function UpdateCellList!(
    x::ParticleSystemPositions,
    y::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl_pair::CellListPair{N,T},
    aux::AuxNonPeriodicPair{N,T};
    parallel::Bool=true,
    validate_coordinates::F=_validate_coordinates,
) where {N,T,F<:Function}
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    ref_list = cl_pair.ref_list
    target_list = cl_pair.target_list
    if x.updated[]
        UpdateCellList!(x, box, ref_list, aux.ref_list; parallel, validate_coordinates)
        x.updated[] = false
    end
    if y.updated[]
        UpdateCellList!(y, box, target_list, aux.target_list; parallel, validate_coordinates)
        y.updated[] = false
    end
    return CellListPair{N,T}(ref_list, target_list)
end

# No-aux CellListPair version: allocates AuxNonPeriodicPair and delegates.
function UpdateCellList!(
    x::ParticleSystemPositions,
    y::ParticleSystemPositions,
    box::Box{NonPeriodicCell},
    cl_pair::CellListPair;
    parallel::Bool=true,
    kargs...
)
    cl_pair = if parallel
        aux = AuxNonPeriodicPair(cl_pair)
        UpdateCellList!(x, y, box, cl_pair, aux; parallel, kargs...)
    else
        UpdateCellList!(x, y, box, cl_pair, nothing; parallel, kargs...)
    end
    cl_pair = update_number_of_batches!(cl_pair; parallel)
    return cl_pair
end

#
# Inner loop dispatch: forward neighbor iteration avoids double-counting without
# needing the real/image index checks used in the periodic path.
#
inner_loop!(f::F, box::Box{NonPeriodicCell}, cellᵢ, cl::CellList, output, ibatch) where {F<:Function} =
    inner_loop!(f, neighbor_cells_forward, box, cellᵢ, cl, output, ibatch)

#
# Self interactions within a single cell: upper-triangle loop, no real/image checks.
#
function _current_cell_interactions!(box::Box{NonPeriodicCell}, f::F, cell, output) where {F<:Function}
    (; cutoff_sqr, inv_rotation) = box
    for i in 1:(cell.n_particles-1)
        pᵢ = cell.particles[i]
        xpᵢ = pᵢ.coordinates
        xpᵢ_rot = inv_rotation * xpᵢ
        for j in (i+1):cell.n_particles
            pⱼ = cell.particles[j]
            d2 = sum(abs2, xpᵢ - pⱼ.coordinates)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, xpᵢ_rot, inv_rotation * pⱼ.coordinates, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end

#
# Cross interactions between two cells (CellListPair): all pairs, no real/image checks.
#
function _current_cell_interactions!(
    box::Box{NonPeriodicCell}, f::F, cellᵢ::Cell, cellⱼ::Cell, output
) where {F<:Function}
    (; cutoff_sqr, inv_rotation) = box
    for i in 1:cellᵢ.n_particles
        pᵢ = cellᵢ.particles[i]
        xpᵢ = pᵢ.coordinates
        xpᵢ_rot = inv_rotation * xpᵢ
        for j in 1:cellⱼ.n_particles
            pⱼ = cellⱼ.particles[j]
            d2 = sum(abs2, xpᵢ - pⱼ.coordinates)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, xpᵢ_rot, inv_rotation * pⱼ.coordinates, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end

#
# Vicinal cell interactions: no real/image checks needed; Skip is irrelevant because
# forward-only neighbor iteration already prevents double-counting for self computations.
#
function _vinicial_cells!(
    f::F, box::Box{NonPeriodicCell}, cellᵢ, pp, Δc, output, ::Val
) where {F<:Function}
    (; cutoff, cutoff_sqr, inv_rotation) = box
    for i in 1:cellᵢ.n_particles
        pᵢ = cellᵢ.particles[i]
        xpᵢ = pᵢ.coordinates
        xpᵢ_rot = inv_rotation * xpᵢ
        xproj = dot(xpᵢ - cellᵢ.center, Δc)
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)
        for j in 1:n
            pⱼ = pp[j]
            d2 = sum(abs2, xpᵢ - pⱼ.coordinates)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, xpᵢ_rot, inv_rotation * pⱼ.coordinates, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end
