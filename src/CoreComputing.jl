#
# Parallel thread spliiter
#
splitter(first,nbatches,n) = first:nbatches:n

"""

```
reduce(output, output_threaded)
```

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

"""
reduce(output::Number, output_threaded::Vector{<:Number}) = sum(output_threaded)
function reduce(output::AbstractArray, output_threaded::AbstractVector{<:AbstractArray}) 
    @. output = output_threaded[1]
    for ibatch in 2:length(output_threaded)
        @. output += output_threaded[ibatch]
    end
    return output
end

"""

```
partition!(x::AbstractVector,by)
```

$(INTERNAL)

# Extended help

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

"""
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
        cell??? = cl.cells[cl.cell_indices_real[i]]
        output = inner_loop!(f, box, cell???, cl, output, ibatch) 
        _next!(p)
    end
    return output
end

#
# Parallel version for self-pairwise computations
#
function batch(f::F, ibatch, nbatches, n_cells_with_real_particles, output_threaded, box, cl, p) where F
    for i in splitter(ibatch, nbatches, n_cells_with_real_particles)
        cell??? = cl.cells[cl.cell_indices_real[i]]
        output_threaded[ibatch] = inner_loop!(f, box, cell???, cl, output_threaded[ibatch], ibatch) 
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
    nbatches = cl.nbatches.map_computation
    if isnothing(output_threaded) 
        output_threaded = [ deepcopy(output) for i in 1:nbatches ] 
    end
    @unpack n_cells_with_real_particles = cl
    nbatches = cl.nbatches.map_computation
    p = show_progress ? Progress(n_cells_with_real_particles, dt=1) : nothing
    @sync for ibatch in 1:nbatches
        Threads.@spawn batch($f, $ibatch, $nbatches, $n_cells_with_real_particles, $output_threaded, $box, $cl, $p)
    end 
    return reduce(output, output_threaded)
end

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(
    f::F, output, box::Box, 
    cl::CellListPair{N,T}; 
    show_progress::Bool=false
) where {F,N,T}
    p = show_progress ? Progress(length(cl.ref), dt=1) : nothing
    for i in eachindex(cl.ref)
        output = inner_loop!(f, output, i, box, cl)
        _next!(p)
    end
    return output
end

#
# Parallel version for cross-interaction computations
#
function batch(f::F, ibatch, nbatches, output_threaded, box, cl, p) where F
    for i in splitter(ibatch, nbatches, length(cl.ref))
        output_threaded[ibatch] = inner_loop!(f, output_threaded[ibatch], i, box, cl) 
        _next!(p)
    end
end

function map_pairwise_parallel!(
    f::F1, output, box::Box, 
    cl::CellListPair{N,T};
    output_threaded=nothing,
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2,N,T}
    nbatches = cl.target.nbatches.map_computation
    if isnothing(output_threaded) 
        output_threaded = [ deepcopy(output) for i in 1:nbatches ]
    end
    p = show_progress ? Progress(length(cl.ref), dt=1) : nothing
    @sync for ibatch in 1:nbatches
        Threads.@spawn batch($f, $ibatch, $nbatches, $output_threaded, $box, $cl, $p)
    end 
    return reduce(output, output_threaded)
end

#
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
function inner_loop!(
    f,box::Box{TriclinicCell},cell???,
    cl::CellList{N,T},
    output,
    ibatch
) where {N,T}
    @unpack cutoff, cutoff_sq, nc = box

    for neighbor_cell in neighbor_cells(box)
        jc_cartesian = cell???.cartesian_index + neighbor_cell
        jc_linear = cell_linear_index(nc,jc_cartesian)
        # if cell??? is empty, cycle
        if cl.cell_indices[jc_linear] == 0
            continue
        end
        cell??? = cl.cells[cl.cell_indices[jc_linear]]

        # same cell
        if cell???.linear_index == cell???.linear_index
            for i in 1:cell???.n_particles
                @inbounds p??? = cell???.particles[i]
                (!p???.real) && continue
                xp??? = p???.coordinates
                for j in 1:cell???.n_particles
                    @inbounds p??? = cell???.particles[j]
                    if p???.index < p???.index
                        xp??? = p???.coordinates
                        d2 = norm_sqr(xp??? - xp???)
                        if d2 <= cutoff_sq
                            output = f(xp???, xp???, p???.index, p???.index, d2, output)
                        end
                    end
                end
            end
        # neighbor cells
        else
            # Vector connecting cell centers
            ??c = cell???.center - cell???.center 
            ??c_norm = norm(??c)
            ??c = ??c / ??c_norm
            pp = project_particles!(cl.projected_particles[ibatch],cell???,cell???,??c,??c_norm,box)
            for i in 1:cell???.n_particles
                @inbounds p??? = cell???.particles[i]
                (!p???.real) && continue
                xp??? = p???.coordinates

                # Partition pp array according to the current projections. This
                # is faster than it looks, because after the first partitioning, the
                # array will be almost correct, and essentially the algorithm
                # will only run over most elements. Avoiding this partitioning or
                # trying to sort the array before always resulted to be much more
                # expensive. 
                xproj = dot(xp??? - cell???.center, ??c)
                n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)

                for j in 1:n
                    @inbounds p??? = pp[j]
                    if p???.index < p???.index
                        xp??? = p???.coordinates
                        d2 = norm_sqr(xp??? - xp???)
                        if d2 <= cutoff_sq
                            output = f(xp???, xp???, p???.index, p???.index, d2, output)
                        end
                    end
                end
            end
        end
    end

    return output
end

#
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
function inner_loop!(
    f,box::Box{OrthorhombicCell},cell???,
    cl::CellList{N,T},
    output,
    ibatch 
) where {N,T}
    @unpack cutoff_sq = box

    # loop over list of non-repeated particles of cell ic
    for i in 1:cell???.n_particles - 1
        @inbounds p??? = cell???.particles[i]
        xp??? = p???.coordinates
        for j in i + 1:cell???.n_particles
            @inbounds p??? = cell???.particles[j]
            xp??? = p???.coordinates
            d2 = norm_sqr(xp??? - xp???)
            if d2 <= cutoff_sq
                output = f(xp???, xp???, p???.index, p???.index, d2, output)
            end
        end
    end

    for jcell in neighbor_cells_forward(box)
        jc_linear = cell_linear_index(box.nc,cell???.cartesian_index + jcell)
        if cl.cell_indices[jc_linear] != 0
            cell??? = cl.cells[cl.cell_indices[jc_linear]]
            output = cell_output!(f, box, cell???, cell???, cl, output, ibatch)
        end
    end

    return output
end

function cell_output!(
    f,
    box::Box{OrthorhombicCell},
    cell???,
    cell???,
    cl::CellList{N,T},
    output,
    ibatch
) where {N,T}
    @unpack cutoff, cutoff_sq, nc = box

    # Vector connecting cell centers
    ??c = cell???.center - cell???.center 
    ??c_norm = norm(??c)
    ??c = ??c / ??c_norm
    pp = project_particles!(cl.projected_particles[ibatch],cell???,cell???,??c,??c_norm,box)
    if length(pp) == 0
        return output
    end

    # Loop over particles of cell icell
    for i in 1:cell???.n_particles
        @inbounds p??? = cell???.particles[i]
        xp??? = p???.coordinates
        xproj = dot(xp??? - cell???.center, ??c)
    
        # Partition pp array according to the current projections
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)

        # Compute the interactions 
        for j in 1:n 
            @inbounds p??? = pp[j]
            xp??? = p???.coordinates
            d2 = norm_sqr(xp??? - xp???)
            if d2 <= cutoff_sq
                output = f(xp???, xp???, p???.index, p???.index, d2, output)
            end
        end
    end

    return output
end

"""

```
project_particles!(projected_particles,cell???,cell???,??c,??c_norm,box)
```

$(INTERNAL)

# Extended help

Projects all particles of the cell `cell???` into unnitary vector `??c` with direction 
connecting the centers of `cell???` and `cell???`. Modifies `projected_particles`, and 
returns a view of `projected particles, where only the particles for which
the projection on the direction of the cell centers still allows the particle
to be within the cutoff distance of any point of the other cell.

"""
function project_particles!(
    projected_particles,cell???,cell???,
    ??c,??c_norm,box::Box{UnitCellType,N}
) where {UnitCellType,N}
    if box.lcell == 1
        margin = box.cutoff + ??c_norm/2 # half of the distance between centers
    else
        margin = box.cutoff*(1 + sqrt(N)/2) # half of the diagonal of the cutoff-box
    end
    iproj = 0
    @inbounds for j in 1:cell???.n_particles
        p??? = cell???.particles[j]
        xproj = dot(p???.coordinates - cell???.center, ??c)
        if abs(xproj) <= margin
            iproj += 1
            projected_particles[iproj] = ProjectedParticle(p???.index, xproj, p???.coordinates) 
        end
    end
    pp = @view(projected_particles[1:iproj])
    return pp
end

#
# Inner loop of cross-interaction computations. If the two sets
# are large, this might(?) be improved by computing two separate
# cell lists and using projection and partitioning.
#
function inner_loop!(
    f,output,i,box,
    cl::CellListPair{N,T}
) where {N,T}
    @unpack nc, cutoff_sq = box
    xp??? = wrap_to_first(cl.ref[i], box)
    ic = particle_cell(xp???, box)
    for neighbor_cell in neighbor_cells(box)
        jc_cartesian = neighbor_cell + ic
        jc_linear = cell_linear_index(nc,jc_cartesian) 
        # If cell??? is empty, cycle
        if cl.target.cell_indices[jc_linear] == 0
            continue
        end
        cell??? = cl.target.cells[cl.target.cell_indices[jc_linear]]
        # loop over particles of cell???
        for j in 1:cell???.n_particles
            @inbounds p??? = cell???.particles[j]
            xp??? = p???.coordinates
            d2 = norm_sqr(xp??? - xp???)
            if d2 <= cutoff_sq
                output = f(xp???, xp???, i, p???.index, d2, output)
            end
        end                                   
    end
    return output
end
