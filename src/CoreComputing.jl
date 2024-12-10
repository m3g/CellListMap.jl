#=
    reduce(output, output_threaded)

Most common reduction function, which sums the elements of the output. 
Here, `output_threaded` is a vector containing `nbatches(cl)` copies of
the `output` variable (a scalar or an array). Custom reduction functions 
must replace this one if the reduction operation is not a simple sum. 
The `output_threaded` array is, by default, created automatically by copying
the given `output` variable `nbatches(cl)` times. 

## Examples

Scalar reduction: 

```julia-repl
julia> output = 0.; output_threaded = [ 1, 2 ];

julia> CellListMap.reduce(output,output_threaded)
3
```
 
Array reduction:

```julia-repl
julia> output = [0,0]; output_threaded = [ [1,1], [2,2] ];

julia> CellListMap.reduce(output,output_threaded)
2-element Vector{Int64}:
 3
 3

julia> output
2-element Vector{Int64}:
 3
 3
```

=#
reduce(output::T, output_threaded::Vector{T}) where {T} = sum(output_threaded)
function reduce(output::T, output_threaded::Vector{T}) where {T<:AbstractArray}
    for ibatch in eachindex(output_threaded)
        @. output += output_threaded[ibatch]
    end
    return output
end
function reduce(output, output_threaded)
    T = typeof(output)
    throw(ArgumentError("""\n
    No method matching reduce(::$(typeof(output)),::$(typeof(output_threaded)))

    Please provide a method that appropriately reduces a `Vector{$T}`, with
    the signature:

    ```
    custom_reduce(output::$T, output_threaded::Vector{$T})
    ```

    The reduction function **must** return the `output` variable, even 
    if it is mutable.  

    See: https://m3g.github.io/CellListMap.jl/stable/parallelization/#Custom-reduction-functions

    """))
end

#=
    partition!(by, x::AbstractVector)

# Extended help

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

=#
function partition!(by, x::AbstractVector)
    iswap = 1
    @inbounds for i in eachindex(x)
        if by(x[i])
            if iswap != i
                x[iswap], x[i] = x[i], x[iswap]
            end
            iswap += 1
        end
    end
    return iswap - 1
end

#
# Auxiliary functions to control the exibition of the progress meter
#
_next!(::Nothing) = nothing
_next!(p) = ProgressMeter.next!(p)

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
function batch(f::F, ibatch, cell_indices, output_threaded, box, cl, p) where {F}
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
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(
    f::F, output, box::Box, cl::CellListPair{N,T};
    show_progress::Bool=false
) where {F,N,T}
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

# Parallel version for self-pairwise computations
#
function batch(f::F, ibatch, cell_indices, output_threaded, box, cl, p) where {F}
    for i in cell_indices
        cellᵢ = cl.small_set[cl.cell_indices_real[i]]
        output_threaded[ibatch] = inner_loop!(f, box, cellᵢ, cl.large_set, output_threaded[ibatch], ibatch)
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
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
inner_loop!(f::F, box::Box{<:OrthorhombicCellType}, cellᵢ, cl::CellList, output, ibatch) where {F<:Function} =
    inner_loop!(f, neighbor_cells_forward, box, cellᵢ, cl, output, ibatch)
inner_loop!(f::F, box::Box{<:TriclinicCell}, cellᵢ, cl::CellList, output, ibatch) where {F<:Function} =
    inner_loop!(f, neighbor_cells, box, cellᵢ, cl, output, ibatch)

#
# Inner loop for self-interaction computations: the call to the current_cell function
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
# Inner loop for cross-interaction computations: the call to the current_cell function
# has two cells as input, and vicinal cell interactions do not skip the i>=j in any type
# of cell. 
#
function inner_loop!(
    f::F,
    _neighbor_cells::NC, # depends on cell type
    box,
    cellᵢ,
    cl::CellListPair{N,T},
    output,
    ibatch
) where {F<:Function,NC<:Function,N,T}
    @unpack cutoff_sqr, inv_rotation, nc = box
    output = _current_cell_interactions!(box, f, cellᵢ, cellⱼ, output)
    for jcell in _neighbor_cells(box)
        jc_linear = cell_linear_index(nc, cellᵢ.cartesian_index + jcell)
        if cl.cell_indices[jc_linear] != 0
            cellⱼ = cl.cells[cl.cell_indices[jc_linear]]
            output = _vicinal_cell_interactions!(CrossVicinal, f, box, cellᵢ, cellⱼ, cl, output, ibatch)
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
    # loop over list of non-repeated particles of cell ic. All particles are real
    for i in 1:cell.n_particles-1
        @inbounds pᵢ = cell.particles[i]
        xpᵢ = pᵢ.coordinates
        for j in i+1:cell.n_particles
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
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, pᵢ.index, pⱼ.index, d2, output)
            end
        end
    end
    return output
end

# Providing two cells for this function indicates that this is a cross-interaction, thus we need
# to loop over all pairs of particles.
function _current_cell_interactions!(box::Box{TriclinicCell}, f::F, cellᵢ, cellⱼ, output) where {F<:Function}
    @unpack cutoff_sqr, inv_rotation = box
    # loop over all pairs, skip when i >= j, skip if neither particle is real
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

#
# Compute interactions between vicinal cells
#
function _vicinal_cell_interactions!(f::F, box::Box, cellᵢ, cellⱼ, cl::CellList{N,T}, output, ibatch; skip=nothing) where {F<:Function,N,T}
    # project particles in vector connecting cell centers
    Δc = cellⱼ.center - cellᵢ.center
    Δc_norm = norm(Δc)
    Δc = Δc / Δc_norm
    pp = project_particles!(cl.projected_particles[ibatch], cellⱼ, cellᵢ, Δc, Δc_norm, box)
    if length(pp) > 0
        output = _vinicial_cells!(f, box, cellᵢ, pp, Δc, output; skip)
    end
    return output
end

#
# The criteria form skipping computations is different then in Orthorhombic or Triclinic boxes. The 
# first parameter (the type Self, or Cross, computation, is not needed here, because the symmetry
# allows to never compute repeated interactions anyway. 
#
function _vinicial_cells!(f::F, box::Box{<:OrthorhombicCellType}, cellᵢ, pp, Δc, output; skip=nothing) where {F<:Function}
    @unpack cutoff, cutoff_sqr, inv_rotation = box
    # Loop over particles of cell icell
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        # project particle in vector connecting cell centers
        xpᵢ = pᵢ.coordinates
        xproj = dot(xpᵢ - cellᵢ.center, Δc)
        # Partition pp array according to the current projections
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)
        # Compute the interactions 
        for j in 1:n
            @inbounds pⱼ = pp[j]
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, pᵢ.index, pⱼ.index, d2, output)
            end
        end
    end
    return output
end

# Here skip determines if the interactions are self or cross, in such a way
# that, for self-computations, we need to skip the interactions when i >= j.
function _vinicial_cells!(f::F, box::Box{<:TriclinicCell}, cellᵢ, pp, Δc, output; skip=nothing) where {F<:Function}
    @unpack cutoff, cutoff_sqr, inv_rotation = box
    # Loop over particles of cell icell
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        # project particle in vector connecting cell centers
        xpᵢ = pᵢ.coordinates
        xproj = dot(xpᵢ - cellᵢ.center, Δc)
        # Partition pp array according to the current projections
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)
        # Compute the interactions 
        pᵢ.real || continue
        for j in 1:n
            @inbounds pⱼ = pp[j]
            if skip === true
                pᵢ.index >= pⱼ.index && continue
            end
            xpⱼ = pⱼ.coordinates
            d2 = norm_sqr(xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                output = f(inv_rotation * xpᵢ, inv_rotation * xpⱼ, pᵢ.index, pⱼ.index, d2, output)
            end
        end
    end
    return output
end

#=
    project_particles!(projected_particles,cellⱼ,cellᵢ,Δc,Δc_norm,box)

# Extended help

Projects all particles of the cell `cellⱼ` into unnitary vector `Δc` with direction 
connecting the centers of `cellⱼ` and `cellᵢ`. Modifies `projected_particles`, and 
returns a view of `projected particles, where only the particles for which
the projection on the direction of the cell centers still allows the particle
to be within the cutoff distance of any point of the other cell.

=#
function project_particles!(
    projected_particles, cellⱼ, cellᵢ,
    Δc, Δc_norm, box::Box{UnitCellType,N}
) where {UnitCellType,N}
    if box.lcell == 1
        margin = box.cutoff + Δc_norm / 2 # half of the distance between centers
    else
        margin = box.cutoff * (1 + sqrt(N) / 2) # half of the diagonal of the cutoff-box
    end
    iproj = 0
    @inbounds for j in 1:cellⱼ.n_particles
        pⱼ = cellⱼ.particles[j]
        xproj = dot(pⱼ.coordinates - cellᵢ.center, Δc)
        if abs(xproj) <= margin
            iproj += 1
            projected_particles[iproj] = ProjectedParticle(pⱼ.index, xproj, pⱼ.coordinates)
        end
    end
    pp = @view(projected_particles[1:iproj])
    return pp
end

