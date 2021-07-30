
"""

```
UpdateCellList!(
  x::AbstractVector{SVector{N,T}},y::AbstractVector{SVector{N,T}},
  box::Box,cl:CellListPair,parallel=true
) where {N,T}
```

Function that will update a previously allocated `CellListPair` structure, given new updated particle positions, for example.

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> y = [ 250*rand(SVector{3,Float64}) for i in 1:10000 ];

julia> cl = CellList(x,y,box);

julia> cl = UpdateCellList!(x,y,box,cl); # update lists

```

"""
function UpdateCellList!(
  x::AbstractVector{SVector{N,T}},
  y::AbstractVector{SVector{N,T}},
  box::Box,cl_pair::CellListPair;
  parallel::Bool=true
) where {N,T}

  if length(x) <= length(y)
    UpdateCellList!(y,box,cl_pair.large,parallel=parallel)
  else
    UpdateCellList!(x,box,cl_pair.large,parallel=parallel)
  end

  return cl_pair
end

"""

```
partition!(x::AbstractVector,by)
```

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

"""
function partition!(x::AbstractVector,by)
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

# Chunck splitter, we jump over cells such that possible 
# heterogeneities have greater changes of getting split into different chuncks
@inline splitter(it,n) = it:nthreads():n

#
# Functions to reduce the output of common options (vectors of numbers 
# and vectors of vectors)
#
reduce(output::Number, output_threaded::Vector{<:Number}) = sum(output_threaded)
function reduce(output::AbstractVector, output_threaded::AbstractVector{<:AbstractVector}) 
  for i in 1:nthreads()
    @. output += output_threaded[i]
  end
  return output
end

"""

```
cell_center(c::CartesianIndex{N},box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}
```

Computes the geometric center of a computin cell, to be used in the projection
of points. Returns a `SVector{N,T}`

"""
@inline cell_center(c::CartesianIndex{N},box::Box{UnitCellType,N,T}) where {UnitCellType,N,T} =
  SVector{N,T}(ntuple(i -> box.cell_size[i]*(c[i] - 0.5 - box.lcell), N))

#
# Currently the methods for cross-interaction computations are the
# same for all types of systems
#

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, 
  cl::CellListPair{SystemType,N,T}; 
  show_progress=show_progress
) where {F,SystemType,N,T}
  show_progress && (p = Progress(length(cl.small),dt=0))
  for i in eachindex(cl.small)
    output = inner_loop!(f,output,i,box,cl)
    show_progress && next!(p)
  end
  return output
end

#
# Parallel version for cross-interaction computations
#
function map_pairwise_parallel!(
  f::F1, output, box::Box, 
  cl::CellListPair{SystemType,N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress=show_progress
) where {F1,F2,SystemType,N,T}
  show_progress && (p = Progress(length(cl.small),dt=1))
  @threads for it in 1:nthreads()
    for i in splitter(it,length(cl.small))
      output_threaded[it] = inner_loop!(f,output_threaded[it],i,box,cl) 
      show_progress && next!(p)
    end
  end 
  output = reduce(output,output_threaded)
  return output
end

#
# Inner loop of cross-interaction computations
#
function inner_loop!(
  f,output,i,box,
  cl::CellListPair{SystemType,N,T}
) where {SystemType,N,T}
  @unpack nc, cutoff_sq = box
  xpᵢ = wrap_to_first(cl.small[i],box)
  ic = particle_cell(xpᵢ,box)
  for neighbour_cell in neighbour_cells_all(box)
    jc = cell_linear_index(nc,neighbour_cell+ic)
    pⱼ = cl.large.fp[jc]
    j = pⱼ.index
    # loop over particles of cell jc
    while j > 0
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        j_orig = pⱼ.index_original 
        output = f(xpᵢ,xpⱼ,i,j_orig,d2,output)
      end
      pⱼ = cl.large.np[j]
      j = pⱼ.index
    end                                   
  end
  return output
end



