"""
    map_pairwise!(f::Function, output, box::Box, cl::CellListPair)

Evaluate function f for pairs in two independent sets of particles, for which the `CellListPair`
object was constructed. 

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
function _current_cell_interactions!(box::Box, f::F, cellᵢ::Cell, cellⱼ::Cell, output) where {F<:Function}
    @unpack cutoff_sqr, inv_rotation = box
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        xpᵢ = pᵢ.coordinates
        pᵢ.real || continue
        for j in 1:cellⱼ.n_particles
            @inbounds pⱼ = cellⱼ.particles[j]
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, pᵢ.index, pⱼ.index, d2, output)
            end
        end
    end
    return output
end

#
# Cross-computations when only one cell list was computed
#
"""
    map_pairwise!(f::Function, x::AbstractVector{<:AbstractVector}, sys::ParticleSystem1; kargs...)
    map_pairwise!(f::Function, x::AbstractMatrix, sys::ParticleSystem1; kargs...)

Evaluate function f for pairs in two independent sets of particles, where the `sys::ParicleSystem1` object
contains the previously computed cell lists of one set of particles, and the second set is given by the
array of positions `x`.

This function can be advantageous over computing the interactions with `CellListPair`, because here the 
cell lists are only computed for one set. This is advantageous in two situations:

    1. The second set of particles is not changing, and the first set is changing. Thus, the cell lists
       of the second set can be computed only once.
    2. One of the sets is much smaller than the other. In this case, computing the cell lists of the largest
       set might be too expensive. Construct the `ParticleSystem` object for the smallest set, and use this
       function to compute the interactions with the largest set.

## Keyword arguments:

- `show_progress::Bool=false`: Show progress bar.
- `update_lists::Bool=true`: Update the cell lists or not. If the positions of the `ParticleSystem1` object
   have not changed, it is not necessary to update the cell lists.

## Example

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 1000);

julia> sys = ParticleSystem(positions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0);

julia> y = rand(SVector{3,Float64}, 100);

julia> map_pairwise((x,y,i,j,d2,output) -> output + sqrt(d2), y, sys; update_lists=false) # Compute the sum of the distances of x and y
31.121496300032163

julia> z = rand(SVector{3,Float64}, 200);

julia> map_pairwise((x,y,i,j,d2,output) -> output + sqrt(d2), z, sys; update_lists=false) # Compute the sum of the distances x and z
63.57860511891242
```

### Note that, in this case, if the computation is run serially, it is completely non-allocating:

```jldoctest
julia> using CellListMap, StaticArrays, BenchmarkTools

julia> sys = ParticleSystem(positions=rand(SVector{3,Float64}, 1000), unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false);

julia> y = rand(SVector{3,Float64}, 100);

julia> f(x,y,i,j,d2,output) = output + sqrt(d2);

julia> @ballocated map_pairwise(\$f, \$y, \$sys; update_lists=false) samples=1 evals=1
0
```

"""
function map_pairwise!(
    f::F, x::AbstractVecOrMat, sys::ParticleSystem1; 
    show_progress::Bool=false, update_lists::Bool=true,
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded)
    UpdateParticleSystem!(sys, update_lists)
    sys.output = map_pairwise!(
        f, sys.output, sys._box, x, sys._cell_list; 
        output_threaded=sys._output_threaded,
        reduce=(output, output_threaded) -> reduce_output!(reducer, output, output_threaded),
        parallel=sys.parallel, show_progress,
    )
    return sys.output
end

#
# The splitting of serial and parallel versions was done to avoid allocations
# associated to the code of creating of `output_threaded` in the serial version, despite
# the fact that it is not used and initialized with `nothing`.
#
function _serial_map_pairwise_x_vs_sys!(
    f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N,T}; 
    show_progress::Bool=false, 
) where {F1<:Function, N, T}
    p = show_progress ? Progress(length(x), dt=1) : nothing
    for i in eachindex(x)
        xs = SVector{N,T}(x[i])
        output = single_particle_vs_list!(f, output, box, i, xs, cl)
        _next!(p)
    end
    return output
end

function _batch_x_vs_sys!(f::F, x, x_atom_indices, ibatch, output_threaded, box, cl, p) where {F<:Function} 
    for i in x_atom_indices
        output_threaded[ibatch] = single_particle_vs_list!(f, output_threaded[ibatch], box, i, x[i], cl)
        _next!(p)
    end
end

function _parallel_map_pairwise_x_vs_sys!(
    f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N,T}; 
    show_progress::Bool=false, output_threaded=nothing, reduce::F2=reduce,
) where {F1<:Function, F2<:Function, N, T}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for _ in 1:_nbatches]
    end
    p = show_progress ? Progress(length(x), dt=1) : nothing
    @sync for (ibatch, x_atom_indices) in enumerate(index_chunks(1:length(x); n=_nbatches, split=Consecutive()))
        @spawn _batch_x_vs_sys!($f, $x, $x_atom_indices, $ibatch, $output_threaded, $box, $cl, $p)
    end
    return reduce(output, output_threaded)
end

function map_pairwise!(
    f::F1, output, box::Box, x::AbstractVector{<:AbstractVector}, cl::CellList{N,T}; 
    parallel::Bool=true, show_progress::Bool=false, output_threaded=nothing, reduce::F2=reduce,
) where {F1<:Function, F2<:Function, N, T}
    output = if parallel
        _parallel_map_pairwise_x_vs_sys!(f, output, box, x, cl; show_progress, output_threaded, reduce)
    else
        _serial_map_pairwise_x_vs_sys!(f, output, box, x, cl; show_progress)
    end
    return output
end

function map_pairwise!(
    f::F1, output, box::Box, x::AbstractMatrix, cl::CellList{N}; 
    parallel::Bool=true, show_progress::Bool=false, output_threaded=nothing, reduce::F2=reduce,
) where {N, F1<:Function, F2<:Function}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return map_pairwise!(f, output, box, x_re, cl; parallel, show_progress, output_threaded, reduce)
end

function single_particle_vs_list!(
    f::F, output, box::Box{<:PeriodicCellType}, 
    i::Integer, x::SVector{N,T},
    cl::CellList{N,T};
) where {F,N,T}
    @unpack nc, cutoff_sqr, inv_rotation = box
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
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, i, pⱼ.index, d2, output)
            end
        end
    end
    return output
end

function single_particle_vs_list!(
    f::F, output, box::Box{NonPeriodicCell}, 
    i::Integer, x::SVector{N,T},
    cl::CellList{N,T};
) where {F,N,T}
    @unpack nc, cutoff_sqr, inv_rotation = box
    xpᵢ = x
    ic = particle_cell(xpᵢ, box)
    for neighbor_cell in current_and_neighbor_cells(box)
        jc_cartesian = neighbor_cell + ic
        # if cell is outside computing box, cycle
        if !all(ntuple(i-> 1 .<= jc_cartesian[i] .<= nc[i], N))
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
                    output = f(xpᵢ, inv_rotation * xpⱼ, i, pⱼ.index, d2, output)
                end
            end
        end
    end
    return output
end