#=
    map_pairwise!(
        f::Function,
        output,
        box::Box,
        cl::CellList
        ;parallel::Bool=true,
        show_progress::Bool=false
    )

This function will run over every pair of particles which are closer than
`box.cutoff` and compute the Euclidean distance between the particles,
considering the periodic boundary conditions given in the `Box` structure.
If the distance is smaller than the cutoff, a function `f` of the
coordinates of the two particles will be computed.

The function `f` receives two arguments as input:
    f(pair, output)
Where `pair` is a `NeighborPair` struct with fields:
- `pair.i`, `pair.j`: indices of the particles
- `pair.x`, `pair.y`: coordinates (minimum-image adjusted)
- `pair.d`: distance between particles
- `pair.d2`: squared distance

The function must return the updated `output` variable.

=#
function map_pairwise!(f::F, output, box::Box, cl::CellList;
    # Parallelization options
    parallel::Bool=true,
    output_threaded=nothing,
    reduce::Function=reduce,
    show_progress::Bool=false,
) where {F} # Needed for specialization for this function (avoids some allocations)
    if parallel
        output = map_pairwise_parallel!(
            f, output, box, cl;
            output_threaded=output_threaded,
            reduce=reduce,
            show_progress=show_progress
        )
    else
        output = map_pairwise_serial!(f, output, box, cl, show_progress=show_progress)
    end
    return output
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
    f::F, output, box::Box, cl::CellList{N,T};
    show_progress::Bool=false
) where {F,N,T}
    @unpack n_cells_with_real_particles = cl
    p = show_progress ? Progress(n_cells_with_real_particles, dt=1) : nothing
    ibatch = 1
    for i in 1:n_cells_with_real_particles
        cellᵢ = cl.cells[cl.cell_indices_real[i]]
        output = inner_loop!(f, box, cellᵢ, cl, output, ibatch)
        _next!(p)
    end
    return output
end

#
# Parallel version for self-pairwise computations
#
function batch(f::F, ibatch, cell_indices, output_threaded, box, cl::CellList, p) where {F}
    for i in cell_indices
        cellᵢ = cl.cells[cl.cell_indices_real[i]]
        output_threaded[ibatch] = inner_loop!(f, box, cellᵢ, cl, output_threaded[ibatch], ibatch)
        _next!(p)
    end
end

# Note: the reason why in the next loop @spawn if followed by interpolated variables
# is to avoid allocations caused by the capturing of variables by the closures created
# by the macro. This may not be needed in the future, if the corresponding issue is solved.
# See: https://discourse.julialang.org/t/type-instability-because-of-threads-boxing-variables/78395
function map_pairwise_parallel!(
    f::F1, output, box::Box, cl::CellList{N,T};
    output_threaded=nothing,
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2,N,T}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for i in 1:_nbatches]
    end
    @unpack n_cells_with_real_particles = cl
    p = show_progress ? Progress(n_cells_with_real_particles, dt=1) : nothing
    @sync for (ibatch, cell_indices) in enumerate(index_chunks(1:n_cells_with_real_particles; n=_nbatches, split=RoundRobin()))
        @spawn batch($f, $ibatch, $cell_indices, $output_threaded, $box, $cl, $p)
    end
    return reduce(output, output_threaded)
end

#
# Inner loop for self-interaction computations (cl::CellList)
#

#
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
inner_loop!(f::F, box::Box{<:OrthorhombicCellType}, cellᵢ, cl::CellList, output, ibatch) where {F<:Function} =
    inner_loop!(f, neighbor_cells_forward, box, cellᵢ, cl, output, ibatch)
inner_loop!(f::F, box::Box{<:TriclinicCell}, cellᵢ, cl::CellList, output, ibatch) where {F<:Function} =
    inner_loop!(f, neighbor_cells, box, cellᵢ, cl, output, ibatch)

#
# The call to the current_cell function
# has a single cell as input, and vicinal cell interactions skip the i>=j for
# triclinic cells.
#
function inner_loop!(
    f::F,
    _neighbor_cells::NC, # depends on cell type
    box,
    cellᵢ,
    cl::CellList{N,T},
    output,
    ibatch
) where {F<:Function,NC<:Function,N,T}
    @unpack cutoff_sqr, inv_rotation, nc = box
    output = _current_cell_interactions!(box, f, cellᵢ, output)
    for jcell in _neighbor_cells(box)
        jc_linear = cell_linear_index(nc, cellᵢ.cartesian_index + jcell)
        if cl.cell_indices[jc_linear] != 0
            cellⱼ = cl.cells[cl.cell_indices[jc_linear]]
            output = _vicinal_cell_interactions!(f, box, cellᵢ, cellⱼ, cl, output, ibatch; skip=true)
        end
    end
    return output
end

#
# Interactions in the current cell
#

# Call with single cell: this implies that this is a self-computation, and thus we loop over the
# upper triangle only in the case of the Orthorhombic cell
function _current_cell_interactions!(box::Box{<:OrthorhombicCellType}, f::F, cell, output) where {F<:Function}
    @unpack cutoff_sqr, inv_rotation = box
    # loop over list of non-repeated particles of cell ic.
    for i in 1:cell.n_particles-1
        @inbounds pᵢ = cell.particles[i]
        xpᵢ = pᵢ.coordinates
        for j in i+1:cell.n_particles
            @inbounds pⱼ = cell.particles[j]
            (pᵢ.real | pⱼ.real) || continue # voltar
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

# And loop over all pairs but skipping when i >= j when the box is triclinic
function _current_cell_interactions!(box::Box{TriclinicCell}, f::F, cell, output) where {F<:Function}
    @unpack cutoff_sqr, inv_rotation = box
    # loop over all pairs, skip when i >= j, skip if neither particle is real
    for i in 1:cell.n_particles
        @inbounds pᵢ = cell.particles[i]
        xpᵢ = pᵢ.coordinates
        pᵢ.real || continue
        for j in 1:cell.n_particles
            @inbounds pⱼ = cell.particles[j]
            (pᵢ.index >= pⱼ.index) && continue
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
