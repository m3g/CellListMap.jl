"""
    map_pairwise!(f::Function,output,box::Box,cl::CellListPair)

The same but to evaluate some function between pairs of the particles of the vectors.

"""
function map_pairwise!(f::F1, output, box::Box, cl::CellListPair{N,T};
    # Parallelization options
    parallel::Bool=true,
    output_threaded=nothing,
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2,N,T} # F1, F2 Needed for specialization for these functions
    fswap(x,y,i,j,d2,output) = f(y,x,j,i,d2,output) 
    if !cl.swap
        if parallel
            output = map_pairwise_parallel!(
                f,output,box,cl;
                output_threaded=output_threaded,
                reduce=reduce,
                show_progress=show_progress
            )
        else
            output = map_pairwise_serial!(f,output,box,cl,show_progress=show_progress)
        end
    else
        if parallel
            output = map_pairwise_parallel!(
                fswap,output,box,cl;
                output_threaded=output_threaded,
                reduce=reduce,
                show_progress=show_progress
            )
        else
            output = map_pairwise_serial!(fswap,output,box,cl,show_progress=show_progress)
        end
    end
    return output
end

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(
    f::F, output, box::Box, cl::CellListPair;
    show_progress::Bool=false
) where {F<:Function}
    @unpack n_cells_with_real_particles = cl.small_set
    p = show_progress ? Progress(n_cells_with_real_particles, dt=1) : nothing
    ibatch = 1
    for i in 1:n_cells_with_real_particles
        cellᵢ = cl.small_set[cl.cell_indices_real[i]]
        output = inner_loop!(f, box, cellᵢ, cl.large_set, output, ibatch)
        _next!(p)
    end
    return output
end

#
# Parallel version for cross computations
#
function batch(f::F, ibatch, cell_indices, output_threaded, box, cl::CellListPair, p) where {F}
    for i in cell_indices
        cellᵢ = cl.small_set[cl.cell_indices_real[i]]
        output_threaded[ibatch] = inner_loop!(f, box, cellᵢ, cl, output_threaded[ibatch], ibatch)
        _next!(p)
    end
end

# Note: the reason why in the next loop @spawn if followed by interpolated variables
# is to avoid allocations caused by the capturing of variables by the closures created
# by the macro. This may not be needed in the future, if the corresponding issue is solved.
# See: https://discourse.julialang.org/t/type-instability-because-of-threads-boxing-variables/78395
function map_pairwise_parallel!(
    f::F1, output, box::Box, cl::CellListPair{N,T};
    output_threaded=nothing,
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2,N,T}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for i in 1:_nbatches]
    end
    @unpack n_cells_with_real_particles = cl.small_set
    p = show_progress ? Progress(n_cells_with_real_particles, dt=1) : nothing
    @sync for (ibatch, cell_indices) in enumerate(index_chunks(1:n_cells_with_real_particles; n=_nbatches, split=RoundRobin()))
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
    cl::CellListPair{N,T},
    output,
    ibatch
) where {F<:Function,N,T}
    @unpack cutoff_sqr, inv_rotation, nc = box
    output = _current_cell_interactions!(box, f, cellᵢ, cellⱼ, output)
    for jcell in _neighbor_cells(box)
        jc_linear = cell_linear_index(nc, cellᵢ.cartesian_index + jcell)
        if cl.large_set.cell_indices[jc_linear] != 0
            cellⱼ = cl.large_set.cells[cl.cell_indices[jc_linear]]
            output = _vicinal_cell_interactions!(f, box, cellᵢ, cellⱼ, cl, output, ibatch)
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
function _current_cell_interactions!(box::Box{TriclinicCell}, f::F, cellᵢ, cellⱼ, output) where {F<:Function}
    @unpack cutoff_sqr, inv_rotation = box
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cell.particles[i]
        xpᵢ = pᵢ.coordinates
        pᵢ.real || continue
        for j in 1:cellⱼ.n_particles
            @inbounds pⱼ = cell.particles[j]
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, pᵢ.index, pⱼ.index, d2, output)
            end
        end
    end
    return output
end