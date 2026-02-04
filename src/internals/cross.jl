#=
    pairwise!(f::Function, output, box::Box, cl::CellListPair)

Evaluate function f for pairs in two independent sets of particles, for which the `CellListPair`
object was constructed.

=#
function _pairwise!(
        f::F1, output, box::Box, cl::CellListPair{N, T};
        # Parallelization options
        parallel::Bool = true,
        output_threaded = nothing,
        reduce::F2 = reduce,
        show_progress::Bool = false
    ) where {F1, F2, N, T} # F1, F2 Needed for specialization for these functions
    fswap(pair, output) = f(NeighborPair(pair.j, pair.i, pair.y, pair.x, pair.d2), output)
    if !cl.swap
        if parallel
            output = _pairwise_parallel!(
                f, output, box, cl;
                output_threaded = output_threaded,
                reduce = reduce,
                show_progress = show_progress
            )
        else
            output = _pairwise_serial!(f, output, box, cl, show_progress = show_progress)
        end
    else
        if parallel
            output = _pairwise_parallel!(
                fswap, output, box, cl;
                output_threaded = output_threaded,
                reduce = reduce,
                show_progress = show_progress
            )
        else
            output = _pairwise_serial!(fswap, output, box, cl, show_progress = show_progress)
        end
    end
    return output
end

#
# Serial version for cross-interaction computations
#
function _pairwise_serial!(
        f::F, output, box::Box, cl::CellListPair;
        show_progress::Bool = false
    ) where {F <: Function}
    (; n_cells_with_real_particles) = cl.small_set
    p = show_progress ? Progress(n_cells_with_real_particles, dt = 1) : nothing
    ibatch = 1
    for i in 1:n_cells_with_real_particles
        cellᵢ = cl.small_set.cells[cl.small_set.cell_indices_real[i]]
        output = inner_loop!(f, box, cellᵢ, cl, output, ibatch)
        _next!(p)
    end
    return output
end

#
# Parallel version for cross computations
#
function batch(f::F, ibatch, cell_indices, output_threaded, box, cl::CellListPair, p) where {F}
    for i in cell_indices
        cellᵢ = cl.small_set.cells[cl.small_set.cell_indices_real[i]]
        output_threaded[ibatch] = inner_loop!(f, box, cellᵢ, cl, output_threaded[ibatch], ibatch)
        _next!(p)
    end
    return
end

# Note: the reason why in the next loop @spawn if followed by interpolated variables
# is to avoid allocations caused by the capturing of variables by the closures created
# by the macro. This may not be needed in the future, if the corresponding issue is solved.
# See: https://discourse.julialang.org/t/type-instability-because-of-threads-boxing-variables/78395
function _pairwise_parallel!(
        f::F1, output, box::Box, cl::CellListPair{N, T};
        output_threaded = nothing,
        reduce::F2 = reduce,
        show_progress::Bool = false
    ) where {F1, F2, N, T}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for i in 1:_nbatches]
    end
    (; n_cells_with_real_particles) = cl.small_set
    p = show_progress ? Progress(n_cells_with_real_particles, dt = 1) : nothing
    @sync for (ibatch, cell_indices) in enumerate(index_chunks(1:n_cells_with_real_particles; n = _nbatches, split = RoundRobin()))
        @spawn batch($f, $ibatch, $cell_indices, $output_threaded, $box, $cl, $p)
    end
    return reduce(output, output_threaded)
end

#
# The interactions do not skip the i>=j in any type of cell.
#
function inner_loop!(
        f::F,
        box,
        cellᵢ,
        cl::CellListPair{N, T},
        output,
        ibatch
    ) where {F <: Function, N, T}
    (; cutoff_sqr, inv_rotation, nc) = box
    jc_linear = cell_linear_index(nc, cellᵢ.cartesian_index)
    if cl.large_set.cell_indices[jc_linear] != 0
        cellⱼ = cl.large_set.cells[cl.large_set.cell_indices[jc_linear]]
        output = _current_cell_interactions!(box, f, cellᵢ, cellⱼ, output)
    end
    for jcell in neighbor_cells(box)
        jc_linear = cell_linear_index(nc, cellᵢ.cartesian_index + jcell)
        if cl.large_set.cell_indices[jc_linear] != 0
            cellⱼ = cl.large_set.cells[cl.large_set.cell_indices[jc_linear]]
            output = _vicinal_cell_interactions!(f, box, cellᵢ, cellⱼ, cl.large_set, output, ibatch)
        end
    end
    return output
end

#
# Interactions in the current cell
#
# Providing two cells for this function indicates that this is a cross-interaction, thus we need
# to loop over all pairs of particles.
#
function _current_cell_interactions!(box::Box, f::F, cellᵢ::Cell, cellⱼ::Cell, output) where {F <: Function}
    (; cutoff_sqr, inv_rotation) = box
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        xpᵢ = pᵢ.coordinates
        pᵢ.real || continue
        for j in 1:cellⱼ.n_particles
            @inbounds pⱼ = cellⱼ.particles[j]
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, inv_rotation * xpᵢ, inv_rotation * xpⱼ, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end

#
# The splitting of serial and parallel versions was done to avoid allocations
# associated to the code of creating of `output_threaded` in the serial version, despite
# the fact that it is not used and initialized with `nothing`.
#
function _serial_pairwise_x_vs_sys!(
        f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N, T};
        show_progress::Bool = false,
    ) where {F1 <: Function, N, T}
    p = show_progress ? Progress(length(x), dt = 1) : nothing
    for i in eachindex(x)
        xs = SVector{N, T}(x[i])
        output = single_particle_vs_list!(f, output, box, i, xs, cl)
        _next!(p)
    end
    return output
end

function _batch_x_vs_sys!(f::F, x, x_atom_indices, ibatch, output_threaded, box, cl, p) where {F <: Function}
    for i in x_atom_indices
        output_threaded[ibatch] = single_particle_vs_list!(f, output_threaded[ibatch], box, i, x[i], cl)
        _next!(p)
    end
    return
end

function _parallel_pairwise_x_vs_sys!(
        f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N, T};
        show_progress::Bool = false, output_threaded = nothing, reduce::F2 = reduce,
    ) where {F1 <: Function, F2 <: Function, N, T}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for _ in 1:_nbatches]
    end
    p = show_progress ? Progress(length(x), dt = 1) : nothing
    @sync for (ibatch, x_atom_indices) in enumerate(index_chunks(1:length(x); n = _nbatches, split = Consecutive()))
        @spawn _batch_x_vs_sys!($f, $x, $x_atom_indices, $ibatch, $output_threaded, $box, $cl, $p)
    end
    return reduce(output, output_threaded)
end

function _pairwise!(
        f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N, T};
        parallel::Bool = true, show_progress::Bool = false, output_threaded = nothing, reduce::F2 = reduce,
    ) where {F1 <: Function, F2 <: Function, N, T}
    output = if parallel
        _parallel_pairwise_x_vs_sys!(f, output, box, x, cl; show_progress, output_threaded, reduce)
    else
        _serial_pairwise_x_vs_sys!(f, output, box, x, cl; show_progress)
    end
    return output
end

function _pairwise!(
        f::F1, output, box::Box, x::AbstractMatrix, cl::CellList{N};
        parallel::Bool = true, show_progress::Bool = false, output_threaded = nothing, reduce::F2 = reduce,
    ) where {N, F1 <: Function, F2 <: Function}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N, eltype(x)}, x)
    return _pairwise!(f, output, box, x_re, cl; parallel, show_progress, output_threaded, reduce)
end

function single_particle_vs_list!(
        f::F, output, box::Box{<:PeriodicCellType},
        i::Integer, x::SVector{N, T},
        cl::CellList{N, T}
    ) where {F, N, T}
    (; nc, cutoff_sqr, inv_rotation) = box
    xpᵢ = box.rotation * wrap_to_first(x, box.input_unit_cell.matrix)
    ic = particle_cell(xpᵢ, box)
    for neighbor_cell in current_and_neighbor_cells(box)
        jc_cartesian = neighbor_cell + ic
        jc_linear = cell_linear_index(nc, jc_cartesian)
        # If cellⱼ is empty, cycle
        if cl.cell_indices[jc_linear] == 0
            continue
        end
        cellⱼ = cl.cells[cl.cell_indices[jc_linear]]
        # loop over particles of cellⱼ
        for j in 1:cellⱼ.n_particles
            @inbounds pⱼ = cellⱼ.particles[j]
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                pair = NeighborPair(i, pⱼ.index, inv_rotation * xpᵢ, inv_rotation * xpⱼ, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end

function single_particle_vs_list!(
        f::F, output, box::Box{NonPeriodicCell},
        i::Integer, x::SVector{N, T},
        cl::CellList{N, T}
    ) where {F, N, T}
    (; nc, cutoff_sqr, inv_rotation) = box
    xpᵢ = x
    ic = particle_cell(xpᵢ, box)
    for neighbor_cell in current_and_neighbor_cells(box)
        jc_cartesian = neighbor_cell + ic
        # if cell is outside computing box, cycle
        if !all(ntuple(i -> 1 .<= jc_cartesian[i] .<= nc[i], N))
            continue
        end
        # If cellⱼ is empty, cycle
        jc_linear = cell_linear_index(nc, jc_cartesian)
        if cl.cell_indices[jc_linear] == 0
            continue
        end
        cellⱼ = cl.cells[cl.cell_indices[jc_linear]]
        # loop over particles of cellⱼ
        for j in 1:cellⱼ.n_particles
            @inbounds pⱼ = cellⱼ.particles[j]
            if pⱼ.real
                xpⱼ = pⱼ.coordinates
                d2 = norm_sqr(xpᵢ - xpⱼ)
                if d2 <= cutoff_sqr
                    pair = NeighborPair(i, pⱼ.index, xpᵢ, inv_rotation * xpⱼ, d2)
                    output = f(pair, output)
                end
            end
        end
    end
    return output
end