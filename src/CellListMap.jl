module CellListMap

using Base.Threads
using Parameters
using StaticArrays
using DocStringExtensions
using ProgressMeter
using Setfield
using LoopVectorization

export Box
export CellList, UpdateCellList!
export map_pairwise!

include("./StructTypes.jl")
include("./LargeDenseSystems.jl")

"""

```
cell_center(c::CartesianIndex{N},cutoff::T) where {N,T}
```

Computes the geometric center of a cell, to be used in the projection
of points. Returns a `SVector{N,T}`

"""
cell_center(c::CartesianIndex{N},cutoff::T) where {N,T} =
  SVector{N,T}(ntuple(i-> cutoff*c[i]/2, N))

#
# Function that checks if the particule is outside the computation bounding box
#
function out_of_bounding_box(x::SVector{N,T},box::Box{N,T,M}) where {N,T,M}
  for i in 1:N
    (x[i] >= box.unit_cell_max[i] + box.cutoff) && return true
    (x[i] < -box.cutoff) && return true
  end
  return false
end
out_of_bounding_box(p::AtomWithIndex,box::Box{N,T,M}) where {N,T,M} =
  out_of_bounding_box(p.coordinates,box)

"""

```
replicate_particle!(ip,p::T,box,cl) where {T <: SVector{2,S} where S}
```

Replicates the particle as many times as necessary to fill the computing box.

"""
function replicate_particle!(ip,p::T,box,cl) where {T <: SVector{2,S} where S}
  @unpack ranges = box
  for i in ranges[1]
    for j in ranges[2]
      i == 0 && j == 0 && continue
      x = translation_image(p,box.unit_cell,(i,j))
      if ! out_of_bounding_box(x,box)
        cl = add_particle_to_celllist!(ip,x,box,cl;real_particle=false) 
      end
    end
  end 
  return cl
end

"""

```
replicate_particle!(ip,p::T,box,cl) where {T <: SVector{3,S} where S}
```

Replicates the particle as many times as necessary to fill the computing box.

"""
function replicate_particle!(ip,p::T,box,cl) where {T <: SVector{3,S} where S}
  @unpack ranges = box
  for i in ranges[1]
    for j in ranges[2]
      for k in ranges[3]
        i == 0 && j == 0 && k == 0 && continue
        x = translation_image(p,box.unit_cell,(i,j,k))
        if ! out_of_bounding_box(x,box)
          cl = add_particle_to_celllist!(ip,x,box,cl;real_particle=false) 
        end
      end
    end
  end 
  return cl
end

function ranges_of_replicas(cutoff,nc,unit_cell,unit_cell_max::SVector{3,T}) where T
  V = SVector{3,T}
  c = cutoff
  um = unit_cell_max
  cell_vertices = SVector{8,V}( 
    V(        -c,        -c,        -c ), 
    V( um[1] + c,        -c,        -c ), 
    V( um[1] + c, um[2] + c,        -c ), 
    V( um[1] + c,        -c, um[3] + c ), 
    V( um[1] + c, um[2] + c, um[3] + c ), 
    V(        -c, um[2] + c,        -c ),
    V(        -c, um[2] + c, um[3] + c ),
    V(        -c,        -c, um[3] + c )
  )
  r_min, r_max = _ranges_of_replicas(
    nc,
    SVector{3,Int}(-1,-1,-1),
    unit_cell,
    cell_vertices
  )
  ranges = SVector{3,UnitRange{Int}}(
    r_min[1]:r_max[1],
    r_min[2]:r_max[2],
    r_min[3]:r_max[3]
  )
  return ranges
end

function ranges_of_replicas(cutoff,nc,unit_cell,unit_cell_max::SVector{2,T}) where T
  V = SVector{2,T}
  c = cutoff
  um = unit_cell_max
  cell_vertices = SVector{4,V}( 
    V(        -c,        -c ), 
    V( um[1] + c,        -c ), 
    V(        -c, um[2] + c ),
    V( um[1] + c, um[2] + c ) 
  )
  r_min, r_max = _ranges_of_replicas(
    nc,
    SVector{2,Int}(-1,-1),
    unit_cell,
    cell_vertices
  )
  ranges = SVector{2,UnitRange{Int}}(
    r_min[1]:r_max[1],
    r_min[2]:r_max[2]
  )
  return ranges
end

function _ranges_of_replicas(r_min,r_max,unit_cell,cell_vertices)
  for vert in cell_vertices
    r = unit_cell \ vert 
    ri = @. ceil(Int,abs(r))
    for (i,el) in pairs(r)
      if el < 0
        @set! ri[i] = -ri[i]
      end
    end
    r_min = min.(ri,r_min)
    r_max = max.(ri,r_max)
  end
  return r_min, r_max
end

"""

```
wrap_cell_fraction(x,cell)
```

`x` is a vector of dimension `N` and `cell` a matrix of dimension `NxN`

"""
@inline function wrap_cell_fraction(x,cell)
  p = rem.(cell\x,1)
  for i in eachindex(p)
    if p[i] < 0
      @set! p[i] += 1
    end
  end
  return p
end

"""

```
wrap_to_first(x::SVector{N,T},cell) where {N,T}
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates. 

"""
@inline function wrap_to_first(x,cell)
  p = wrap_cell_fraction(x,cell)
  return cell*p
end
"""

```
wrap_to_first(x::SVector{N,T},box::Box)
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates, given the `Box` structure.

"""
@inline wrap_to_first(x,box::Box) = wrap_to_first(x,box.unit_cell)

"""

```
wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

"""
function wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
  p = SVector{length(x),eltype(x)}(rem.(cell\(x-xref),1))
  for i in eachindex(p)
    if p[i] > 0.5 
      @set! p[i] -= 1
    elseif p[i] < -0.5
      @set! p[i] += 1
    end
  end
  return cell*p + xref
end

"""

```
translation_image(x::SVector{N,T},unit_cell,indices) where {N,T}
```

Translate vector `x` according to the `unit_cell` lattice vectors and the `indices`
provided.

"""
translation_image(x::SVector{N,T},unit_cell,indices) where {T,N} =
  x + unit_cell*SVector{N,T}(ntuple(i -> indices[i],N))

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
neighbour_cells(box::Box{N,T,M}) where {N,M}
```

Function that returns the iterator of the cartesian indices of the cells that must be 
evaluated (forward, i. e. to avoid repeated interactions) 
if the cells have sides of length `box.cutoff/box.lcell`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells(box::Box{3,T,9}) where T
  @unpack lcell = box
  nb = Iterators.flatten((
    CartesianIndices((1:lcell,-lcell:lcell,-lcell:lcell)),
    CartesianIndices((0:0,1:lcell,-lcell:lcell)),
    CartesianIndices((0:0,0:0,1:lcell))
  ))
  return nb
end
function neighbour_cells(box::Box{2,T,4}) where T
  @unpack lcell = box
  nb = Iterators.flatten((
    CartesianIndices((1:lcell,-lcell:lcell)),
    CartesianIndices((0:0,1:lcell))
  ))
  return nb
end

"""

```
neighbour_cells_all(box::Box{N,T,M}) where {N,M}
```

Function that returns the iterator of the cartesian indices of all neighbouring
cells of a cell if the cells have sides of `box.cutoff/box.lcell`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells_all(box::Box{3,T,9}) where T  
  @unpack lcell = box
  return CartesianIndices((-lcell:lcell,-lcell:lcell,-lcell:lcell))
end
function neighbour_cells_all(box::Box{2,T,4}) where T  
  @unpack lcell = box
  return CartesianIndices((-lcell:lcell,-lcell:lcell))
end

"""

```
norm(v::SVector{N,T}) where {N,T}
```

Computes the norm of a static vector.

"""
@inline function norm(v::AbstractVector{T}) where T
  norm = zero(T)
  for x in v
    norm += x^2
  end
  return sqrt(norm)
end

"""

```
dot(x::SVector{N,T},y::SVector{N,T}) where {N,T}
```

Computes the norm of a static vector.

"""
@inline function dot(x::AbstractVector{T},y::AbstractVector{T}) where T
  @assert length(x) == length(y) 
  dot = zero(T)
  for i in eachindex(x)
    @inbounds dot += x[i]*y[i] 
  end
  return dot
end

"""

```
distance_sq(x,y)
```

Function to compute squared Euclidean distances between two n-dimensional vectors.

"""
@inline function distance_sq(x::AbstractVector{T}, y::AbstractVector{T}) where T
  @assert length(x) == length(y)
  d = zero(T)
  @inbounds for i in eachindex(x)
    d += (x[i]-y[i])^2
  end
  return d
end

"""

```
distance(x,y)
```

Function to compute Euclidean distances between two n-dimensional vectors.

"""
@inline distance(x::AbstractVector{T}, y::AbstractVector{T}) where T =
  sqrt(distance_sq(x,y))

"""

```
particle_cell(x::SVector{N,T}, box::Box) where {N,T}
```

Returns the coordinates of the computing cell to which a particle belongs, given its coordinates
and the cutoff/cell. 

"""
@inline particle_cell(x::SVector{N,T}, box::Box) where {N,T} =
  CartesianIndex(ntuple(i -> floor(Int,(x[i] .+ box.cutoff)/(box.cutoff/box.lcell) + 1), N))

"""

```
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}
```

Given the linear index of the cell in the cell list, returns the cartesian indices 
of the cell (for arbitrary dimension N).

"""
@inline cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} = 
  CartesianIndices(ntuple(i -> nc[i],N))[i1D]

"""

```
cell_linear_index(nc::SVector{N,Int}, indices) where N
```

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

"""
@inline cell_linear_index(nc::SVector{N,Int}, indices) where N =
  LinearIndices(ntuple(i -> nc[i],N))[ntuple(i->indices[i],N)...]

"""

```
map_pairwise!(f::Function,output,box::Box,cl::CellList;parallel::Bool=true,show_progress::Bool=false)
```

This function will run over every pair of particles which are closer than `box.cutoff` and compute
the Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a function `f` of the coordinates
of the two particles will be computed. 

The function `f` receives six arguments as input: 
```
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it. With that, it is possible to compute an average property of the distance of the particles or, for example, build a histogram. The squared distance `d2` is computed internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used for the desired computation. 

## Example

Computing the mean absolute difference in `x` position between random particles, remembering the number of pairs of `n` particles is `n(n-1)/2`. The function does not use the indices or the distance, such that we remove them from the parameters by using a closure.

```julia-repl
julia> n = 100_000;

julia> box = Box([250,250,250],10);

julia> x = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:n ];

julia> cl = CellList(x,box);

julia> f(x,y,sum_dx) = sum_dx + abs(x[1] - y[1])

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> avg_dx = normalization * map_parwise!((x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, box, cl)

```

"""
function map_pairwise!(f::F, output, box::Box, cl::CellList; 
  # Parallelization options
  parallel::Bool=true,
  output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
  reduce::Function=reduce,
  show_progress::Bool=false,
) where {F} # Needed for specialization for this function (avoids some allocations)
  if parallel && nthreads() > 1
    output = map_pairwise_parallel!(f,output,box,cl;
      output_threaded=output_threaded,
      reduce=reduce,
      show_progress=show_progress
    )
  else
    output = map_pairwise_serial!(f,output,box,cl,show_progress=show_progress)
  end
  return output
end

"""

```
map_pairwise!(f::Function,output,box::Box,cl::CellListPair)
```

The same but to evaluate some function between pairs of the particles of the vectors.

"""
function map_pairwise!(f::F1, output, box::Box, cl::CellListPair;
  # Parallelization options
  parallel::Bool=true,
  output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2} # Needed for specialization for this function (avoids some allocations) 
  if cl.swap 
    fswap(x,y,i,j,d2,output) = f(y,x,j,i,d2,output) 
  else
    fswap = f
  end
  if parallel && nthreads() > 1
    output = map_pairwise_parallel!(fswap,output,box,cl;
      output_threaded=output_threaded,
      reduce=reduce,
      show_progress=show_progress
    )
  else
    output = map_pairwise_serial!(fswap,output,box,cl,show_progress=show_progress)
  end
  return output
end

# Chunck splitter, we jump over cells such that possible 
# heterogeneities have greater changes of getting split into different chuncks
splitter(it,n) = it:nthreads():n

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

#
# Test and example functions
#
include("./testing.jl")
include("./examples.jl")
include("./halotools.jl")

end # module



