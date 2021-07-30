#
# TinySystems: even constructing the cell list may be avoided
# also, these functions are very simple and serve for testing
# purposes
#
"""

```
UpdateCellList!(
    x::AbstractVector{SVector{N,T}},
    box::Box,cl:CellList{TinySystem,N,T},
    parallel=true
) where {N,T}
```

Function that will update a previously allocated `CellList` structure, 
given new updated particle positions very small systems.

## Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> cl = CellList(x,box);

julia> box = Box([260,260,260],10);

julia> x = [ 260*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> cl = UpdateCellList!(x,box,cl); # update lists

```

"""
function UpdateCellList!(
  x::AbstractVector{SVector{N,T}},
  box::Box,
  cl::CellList{TinySystem,N,T};
  parallel::Bool=true
) where {N,T}
  @unpack np, ncp, projected_particles = cl

  # This copy is of course redundant, but for the moment
  # it is what it takes to use the same CellList structure 
  ncp[1] = length(x)
  for (ip,coordinates) in pairs(x)
    np[ip] = AtomWithIndex(ip,ip,coordinates)
  end
  resize!(projected_particles[1],ncp[1])

  return cl
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellList{TinySystem,N,T}; 
  show_progress::Bool=false
) where {F,N,T}
  show_progress && println("TinySystems do not show progress.")
  @unpack ncp, np = cl
  @unpack unit_cell, cutoff_sq = box
  for i in 1:ncp[1]-1
    xᵢ = np[i].coordinates
    for j in i+1:ncp[1]
      xⱼ = wrap_relative_to(np[j].coordinates,xᵢ,unit_cell)
      d2 = norm_sqr(xᵢ - xⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,xⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

"""

```
iuppert(k::Integer,n::Integer)
```

Returns the cartesian indexes of of the elements of the upper
triangular part of a square matrix of dimension `n×n`, given
the linear index `k`.

### Examples

```julia-repl
julia> iupper(1,3)
(1,2)

julia> iupper(3,3)
(2,2)

julia> iuppert(3,4)
(1,4)

```

"""
@inline function iuppert(k::Integer,n::Integer)
  i = n - 1 - floor(Int,sqrt(-8*k + 4*n*(n-1) + 1)/2 - 0.5)
  j = k + i + ( (n-i+1)*(n-i) - n*(n-1) )÷2
  return i, j
end

function map_pairwise_parallel!(
  f::F1, output, box::Box, cl::CellList{TinySystem,N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2,N,T}
  @unpack ncp, np = cl
  @unpack unit_cell, cutoff_sq = box
  show_progress && println("TinySystems do not show progress.")
  npairs = ncp[1]*(ncp[1]-1)÷2
  @threads for it in 1:nthreads() 
    for ipair in it:nthreads():npairs
      i, j = iuppert(ipair,ncp[1])  
      xᵢ = np[i].coordinates
      xⱼ = wrap_relative_to(np[j].coordinates,xᵢ,unit_cell)
      d2 = norm_sqr(xᵢ - xⱼ)
      if d2 <= cutoff_sq
        output_threaded[it] = f(xᵢ,xⱼ,i,j,d2,output_threaded[it])
      end
    end
  end

  output = reduce(output,output_threaded)
  return output
end
