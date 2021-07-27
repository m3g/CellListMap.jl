
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
function partition!(x::AbstractVector,by) where T
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
splitter(it,n) = it:nthreads():n

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
cell_center(c::CartesianIndex{N},box::Box{N,T,M}) where {N,T,M}
```

Computes the geometric center of a computin cell, to be used in the projection
of points. Returns a `SVector{N,T}`

"""
@inline cell_center(c::CartesianIndex{N},box::Box{N,T,M}) where {N,T,M} =
  SVector{N,T}(ntuple(i -> box.cell_size*(c[i] - 0.5 - box.lcell), N))


