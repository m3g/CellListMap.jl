module CellListMap

using Base.Threads
using Parameters
using StaticArrays
using DocStringExtensions

export Box
export CellList, UpdateCellList!
export map_pairwise!

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains some data required to compute the linked cells. To
be initialized with the box size and cutoff. 

## Example

```julia-repl
julia> sides = [250,250,250];

julia> cutoff = 10;

julia> box = Box(sides,cutoff)
Box{3, Float64}
  sides: [250.0, 250.0, 250.0]
  cutoff: 10.0
  number of cells on each dimension: [25, 25, 25] (lcell: 1)
  Total number of cells: 15625

```

"""
Base.@kwdef struct Box{N,T}
  sides::SVector{N,T}
  nc::SVector{N,Int}
  cell_side::SVector{N,T}
  cutoff::T
  cutoff_sq::T
  lcell::Int
end
function Box(sides::AbstractVector, cutoff, T::DataType, lcell::Int=1)
  N = length(sides)
  nc = SVector{N,Int}(max.(1,floor.(Int,sides/(cutoff/lcell))))
  l = SVector{N,T}(sides ./ nc)
  box = Box{N,T}(SVector{N,T}(sides),nc,l,cutoff,cutoff^2,lcell)
  return box
end
Box(sides::AbstractVector,cutoff;T::DataType=Float64,lcell::Int=1) =
  Box(sides,cutoff,T,lcell)

function Base.show(io::IO,::MIME"text/plain",box::Box)
  println(typeof(box))
  println("  sides: ", box.sides) 
  println("  cutoff: ", box.cutoff)
  println("  number of cells on each dimension: ",box.nc, " (lcell: ",box.lcell,")")
  print("  Total number of cells: ", prod(box.nc))
end

"""

$(TYPEDEF)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort.

"""
struct AtomWithIndex{N,T}
  index::Int
  coordinates::SVector{N,T}
end

"""

$(TYPEDEF)

This structure contains the cell linear index and the information about if this cell
is in the border of the box (such that its neighbouring cells need to be wrapped) 

"""
struct CellWithBorderInfo
  icell::Int
  inborder::Bool
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef struct CellList{N,T}
  ncwp::Vector{Int} # One-element vector to contain the mutable number of cells with particles
  cwp::Vector{CellWithBorderInfo} # Indexes of the unique cells With Particles
  fp::Vector{AtomWithIndex{N,T}} # First particle of cell 
  np::Vector{AtomWithIndex{N,T}} # Next particle of cell
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
  println(typeof(cl))
  print("  with $(cl.ncwp[1]) cells with particles.")
end

# Structure that will cointain the cell lists of two independent sets of
# particles for cross-computation of interactions
struct CellListPair{V,N,T}
  small::V
  large::CellList{N,T}
end      
function Base.show(io::IO,::MIME"text/plain",cl::CellListPair)
  print(typeof(cl),"\n")
  print("   $(length(cl.small)) particles in the smallest vector.\n")
  print("   $(cl.large.ncwp[1]) cells with particles.")
end
  
"""

```
CellList(x::AbstractVector{SVector{N,T}},box::Box,parallel=true) where {N,T}
```

Function that will initialize a `CellList` structure from scracth, given a vector
or particle coordinates (as `SVector`s) and a `Box`, which contain the size ofthe
system, cutoff, etc.  

### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:100000 ];

julia> cl = CellListMap.CellList(x,box)
CellList{3, Float64}
  with 15597 cells with particles.

```

"""
function CellList(x::AbstractVector{SVector{N,T}},box::Box;parallel::Bool=true) where {N,T} 
  number_of_particles = length(x)
  number_of_cells = ceil(Int,1.1*prod(box.nc)) # some margin in case of box size variations
  ncwp = zeros(Int,1)
  cwp = Vector{CellWithBorderInfo}(undef,number_of_cells)
  fp = Vector{AtomWithIndex{N,T}}(undef,number_of_cells)
  np = Vector{AtomWithIndex{N,T}}(undef,number_of_particles)

  cl = CellList{N,T}(ncwp,cwp,fp,np)
  return UpdateCellList!(x,box,cl,parallel=parallel)
end

"""

```
CellList(x::AbstractVector{SVector{N,T}},y::AbstractVector{SVector{N,T}},box::Box;parallel=true) where {N,T}
```

Function that will initialize a `CellListPair` structure from scracth, given two vectors
of particle coordinates and a `Box`, which contain the size ofthe system, cutoff, etc. The cell lists
will be constructed for the largest vector, and a reference to the smallest vector is annotated.


### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> y = [ 250*rand(SVector{3,Float64}) for i in 1:10000 ];

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   1000 particles in the smallest vector.
   7452 cells with particles.

```

"""
function CellList(
  x::AbstractVector{SVector{N,T}},
  y::AbstractVector{SVector{N,T}},
  box::Box;
  parallel::Bool=true
) where {N,T} 

  if length(x) <= length(y)
    y_cl = CellList(y,box,parallel=parallel)
    cl_pair = CellListPair(x,y_cl)
  else
    x_cl = CellList(x,box,parallel=parallel)
    cl_pair = CellListPair(y,x_cl)
  end

  return cl_pair
end

"""

```
UpdateCellList!(x::AbstractVector{SVector{N,T}},box::Box,cl:CellList{N,T},parallel=true) where {N,T}
```

Function that will update a previously allocated `CellList` structure, given new updated particle positions, for example.

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
  cl::CellList{N,T};
  parallel::Bool=true
) where {N,T}
  @unpack ncwp, cwp, fp, np = cl

  number_of_cells = prod(box.nc)
  if number_of_cells > length(cwp) 
    number_of_cells = ceil(Int,1.1*number_of_cells) # some margin in case of box size variations
    resize!(cwp,number_of_cells)
    resize!(fp,number_of_cells)
  end

  ncwp[1] = 0
  if parallel
    @threads for i in eachindex(cwp)
      cwp[i] = CellWithBorderInfo(0,false)
      fp[i] = AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N)))
    end
    @threads for i in eachindex(np)
      np[i] = AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N)))
    end
  else
    fill!(cwp,CellWithBorderInfo(0,false))
    fill!(fp,AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N))))
    fill!(np,AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N)))) 
  end
  # Not worth paralellizing probably (would need to take care of concurrency)
  for (ip,xip) in pairs(x)
    set_celllist_index!(ip,xip,box,cl)
  end

  return cl
end

#
# Set one index of a cell list
#
function set_celllist_index!(ip,xip,box,cl)
  @unpack ncwp, cwp, fp, np = cl
  p = AtomWithIndex(ip,wrapone(xip,box.sides))
  icell_cartesian = particle_cell(p.coordinates,box)
  icell = cell_linear_index(box.nc,icell_cartesian)
  if fp[icell].index == 0
    ncwp[1] += 1
    cwp[ncwp[1]] = CellWithBorderInfo(icell,cell_in_border(icell_cartesian,box))
  end
  np[ip] = fp[icell]
  fp[icell] = p
end

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
cell_in_border(icell_cartesian,box)
```

Function that checks if a cell is in the border of the periodic cell box

"""
function cell_in_border(icell_cartesian,box)
  if icell_cartesian[1] <= box.lcell ||
     icell_cartesian[2] <= box.lcell ||
     icell_cartesian[3] <= box.lcell ||
     icell_cartesian[1] > box.nc[1]-box.lcell ||
     icell_cartesian[2] > box.nc[2]-box.lcell ||
     icell_cartesian[3] > box.nc[3]-box.lcell
    return true
  else
    return false
  end
end

"""

```
neighbour_cells(lcell::Int)
```

Function that returns the iterator of the cartesian indices of the cell that must be 
evaluated if the cells have sides of length `cutoff/lcell`. 

"""
function neighbour_cells(lcell)
  nb = Iterators.flatten((
    CartesianIndices((1:lcell,-lcell:lcell,-lcell:lcell)),
    CartesianIndices((0:0,1:lcell,-lcell:lcell)),
    CartesianIndices((0:0,0:0,1:lcell))
  ))
  return nb
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

Returns the coordinates of the cell to which a particle belongs, given its coordinates
and the sides of the periodic box (for arbitrary dimension N).

"""
function particle_cell(x::SVector{N,T}, box::Box) where {N,T}
  # Wrap to origin
  xwrapped = wrapone(x,box.sides)
  cell = CartesianIndex(
    ntuple(i -> floor(Int,(xwrapped[i]+box.sides[i]/2)/box.cell_side[i])+1, N)
  )
  return cell
end

"""

```
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}
```

Given the linear index of the cell in the cell list, returns the cartesian indexes
of the cell (for arbitrary dimension N).

"""
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} = 
  CartesianIndices(ntuple(i -> nc[i],N))[i1D]

"""

```
icell1D(nc::SVector{N,Int}, indexes) where N
```

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

"""
cell_linear_index(nc::SVector{N,Int}, indexes) where N =
  LinearIndices(ntuple(i -> nc[i],N))[ntuple(i->indexes[i],N)...]

"""

```
function wrap!(x::AbstractVector, sides::AbstractVector, center::AbstractVector)
```

Functions that wrap the coordinates They modify the coordinates of the input vector.  
Wrap to a given center of coordinates

"""
function wrap!(x::AbstractVector, sides::AbstractVector, center::AbstractVector)
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides,center)
  end
  return nothing
end

@inline function wrapone(x::AbstractVector, sides::AbstractVector, center::AbstractVector)
  s = @. (x-center)%sides
  s = @. wrapx(s,sides) + center
  return s
end

@inline function wrapx(x,s)
  if x > s/2
    x = x - s
  elseif x < -s/2
    x = x + s
  end
  return x
end

"""

```
wrap!(x::AbstractVector, sides::AbstractVector)
```

Wrap to origin (slightly cheaper).

"""
function wrap!(x::AbstractVector, sides::AbstractVector)
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides)
  end
  return nothing
end

@inline function wrapone(x::AbstractVector, sides::AbstractVector)
  s = @. x%sides
  s = @. wrapx(s,sides)
  return s
end

"""

```
wrap_cell(nc::SVector{N,Int}, indexes) where N
```

Given the dimension `N` of the system, return the periodic cell which correspondst to
it, if the cell is outside the main box.

"""
@inline function wrap_cell(nc::SVector{N,Int}, indexes) where N
  cell_indexes = ntuple(N) do i
    ind = indexes[i]
    if ind < 1
      ind = nc[i] + ind
    elseif ind > nc[i]
      ind = ind - nc[i]
    end
    return ind
  end
  return cell_indexes
end

"""

```
map_pairwise!(f::Function,output,box::Box,cl::CellList;parallel::Bool=true)
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

julia> x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:n ];

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
  reduce::Function=reduce
) where {F} # Needed for specialization for this function (avoids some allocations)
  if parallel && nthreads() > 1
    output = map_pairwise_parallel!(f,output,box,cl;
      output_threaded=output_threaded,
      reduce=reduce
    )
  else
    output = map_pairwise_serial!(f,output,box,cl)
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
  reduce::F2=reduce
) where {F1,F2} # Needed for specialization for this function (avoids some allocations) 
  if parallel && nthreads() > 1
    output = map_pairwise_parallel!(f,output,box,cl;
      output_threaded=output_threaded,
      reduce=reduce
    )
  else
    output = map_pairwise_serial!(f,output,box,cl)
  end
  return output
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(f::F, output, box::Box, cl::CellList) where {F}
  for icell in 1:cl.ncwp[1]
    output = inner_loop!(f,box,icell,cl,output) 
  end 
  return output
end

#
# Parallel version for self-pairwise computations
#
function map_pairwise_parallel!(f::F1, output, box::Box, cl::CellList;
  output_threaded=output_threaded,
  reduce::F2=reduce
) where {F1,F2}
  @threads for icell in 1:cl.ncwp[1]
    it = threadid()
    output_threaded[it] = inner_loop!(f,box,icell,cl,output_threaded[it]) 
  end 
  output = reduce(output,output_threaded)
  return output
end

function inner_loop!(f,box,icell,cl::CellList,output)
  @unpack sides, nc, cutoff_sq = box
  cell = cl.cwp[icell]
  ic = cell.icell
  ic_cartesian = cell_cartesian_indices(nc,ic)

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[ic]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.np[pᵢ.index] 
    j = pⱼ.index
    while j > 0
      if cell.inborder
        xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
      else
        xpⱼ = pⱼ.coordinates
      end
      d2 = distance_sq(xpᵢ,xpⱼ)
      if d2 <= cutoff_sq
        output = f(xpᵢ,xpⱼ,i,j,d2,output)
      end
      pⱼ = cl.np[pⱼ.index]
      j = pⱼ.index
    end
    pᵢ = cl.np[pᵢ.index]
    i = pᵢ.index
  end
   
  for jcell in neighbour_cells(box.lcell)
    output = cell_output!(f,box,cell,cl,output,ic_cartesian+jcell)
  end

  return output
end

#
# loops over the particles of a neighbour cell
#
function cell_output!(f,box,cell,cl,output,jc_cartesian)
  @unpack sides, nc, cutoff_sq = box
  ic = cell.icell
  if cell.inborder
    jc_cartesian_wrapped = wrap_cell(nc,jc_cartesian)
  else
    jc_cartesian_wrapped = jc_cartesian
  end
  jc = cell_linear_index(nc,jc_cartesian_wrapped)

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[ic]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.fp[jc]
    j = pⱼ.index
    while j > 0
      if cell.inborder
        xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
      else
        xpⱼ = pⱼ.coordinates
      end
      d2 = distance_sq(xpᵢ,xpⱼ)
      if d2 <= cutoff_sq
        output = f(xpᵢ,xpⱼ,i,j,d2,output)
      end
      pⱼ = cl.np[pⱼ.index]
      j = pⱼ.index
    end
    pᵢ = cl.np[pᵢ.index]
    i = pᵢ.index
  end

  return output
end

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(f::F, output, box::Box, cl::CellListPair) where {F}
  for i in eachindex(cl.small)
    output = inner_loop!(f,output,i,box,cl)
  end
  return output
end

#
# Parallel version for cross-interaction computations
#
function map_pairwise_parallel!(f::F1, output, box::Box, cl::CellListPair;
  output_threaded=output_threaded,
  reduce::F2=reduce
) where {F1,F2}
  @threads for i in eachindex(cl.small)
    it = threadid()
    output_threaded[it] = inner_loop!(f,output_threaded[it],i,box,cl) 
  end 
  output = reduce(output,output_threaded)
  return output
end

#
# Inner loop of cross-interaction computations
#
function inner_loop!(f,output,i,box,cl::CellListPair)
  @unpack sides, nc, cutoff_sq = box
  xpᵢ = wrapone(cl.small[i],box.sides)
  ic = particle_cell(xpᵢ,box)
  inborder = cell_in_border(ic,box)
  for neighbour_cell in CartesianIndices((-1:1, -1:1, -1:1))   
    if inborder
      jc_cartesian_wrapped = wrap_cell(nc,neighbour_cell+ic)
    else
      jc_cartesian_wrapped = neighbour_cell+ic
    end
    jc = cell_linear_index(nc,jc_cartesian_wrapped)
    pⱼ = cl.large.fp[jc]
    j = pⱼ.index
    # loop over particles of cell jc
    while j > 0
      if inborder
        xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
      else
        xpⱼ = pⱼ.coordinates
      end
      d2 = distance_sq(xpᵢ,xpⱼ)
      if d2 <= cutoff_sq
        output = f(xpᵢ,xpⱼ,i,j,d2,output)
      end
      pⱼ = cl.large.np[j]
      j = pⱼ.index
    end                                   
  end
  return output
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
include("./naive.jl")
include("./examples.jl")
include("./halotools.jl")

end # module



