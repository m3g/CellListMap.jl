module CellListMap

using Base.Threads
using Parameters
using StaticArrays
using DocStringExtensions

export LinkedLists, Box, initlists!
export map_pairwise!

"""

```
LinkedLists(N)
```

Structure that contains the vectors storing the first particle and next particle of the linked
cells. To be initialized with the number of particles.

## Example

```julia-repl
julia> lc = LinkedLists(100_000)

julia> lc = LinkedLists(100_000)
LinkedLists{100000}
  firstatom: Array{Int64}((100000,)) [0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  nextatom: Array{Int64}((100000,)) [0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  icell: Array{Int64}((100000,)) [0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

```

"""  
@with_kw struct LinkedLists{N}
  firstatom::Vector{Int} = zeros(Int,N)
  nextatom::Vector{Int} = zeros(Int,N)
  icell::Vector{Int} = zeros(Int,N)
end
LinkedLists(N::Int) = LinkedLists{N}()


"""

$(TYPEDEF)

Structure that contains some data required to compute the linked cells. To
be initialized with the box size and cutoff. An optional parameter `lcell` 
can be provided as the last argument to define the cell size relative to the
cutoff (default 2, meaning half of the cutoff).

## Example

```julia-repl
julia> sides = [250,250,250];

julia> cutoff = 10.;

julia> box = Box(sides,cutoff)

julia> box = Box(sides,cutoff)
Box{3}
  sides: SVector{3, Int64}
  nc: SVector{3, Int64}
  l: SVector{3, Float64}
  lcell: Int64 2
  cutoff: Float64 10.0
  cutoff_sq: Float64 100.0

```


"""
Base.@kwdef struct Box{N,T}
  sides::SVector{N,T}
  nc::SVector{N,Int}
  cell_side::SVector{N,T}
  lcell::Int
  cutoff::T
  cutoff_sq::T
end
function Box(sides::AbstractVector, cutoff, lcell::Int, T::DataType)
  N = length(sides)
  nc = SVector{N,Int}(max.(1,trunc.(Int,sides/(cutoff/lcell))))
  l = SVector{N,T}(sides ./ nc)
  return Box{N,T}(SVector{N,T}(sides),nc,l,lcell,cutoff,cutoff^2)
end
Box(sides::AbstractVector,cutoff;lcell=Int=2,T::DataType=Float64) =
  Box(sides,cutoff,lcell,T)

"""

$(TYPEDEF)

Copies particle coordinates and associated index, to build contiguous atom lists
in memory when building the cell lists. This strategy duplicates the atomic coordinates
data, but is probably worth the effort.

"""
struct AtomWithIndex{N,T}
  index::Int
  coordinates::SVector{N,T}
end

"""

```
CellLists(x::AbstractVector{SVector{N,T}},box::Box{N}) where {N,T}
```

Function that will initialize a `CellLists` structure from scracth, given a vector
or particle coordinates (as `SVector`s) and a `Box`, which contain the size ofthe
system, cutoff, etc.  

"""
Base.@kwdef struct CellLists{N,T}
  ncwp::Int # Number of cells with particles
  cwp::Vector{Int} # Indexes of the unique cells With Particles
  ncp::Vector{Int} # Number of cell particles
  fp::Vector{AtomWithIndex{N,T}} # First particle of cell 
  np::Vector{AtomWithIndex{N,T}} # Next particle of cell
end

function CellLists(x::AbstractVector{SVector{N,T}},box::Box{N}) where {N,T}
  #
  # Initialize the list with the first particle (just to avoid redundant testing of empty arrays)
  #
  number_of_particles = length(x)
  number_of_cells = prod(box.nc)
  ncwp = 0
  cwp = zeros(Int,number_of_cells)
  ncp = zeros(Int,number_of_cells)
  fp = fill(AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N))),number_of_cells)
  np = fill(AtomWithIndex{N,T}(0,SVector{N,T}(ntuple(i->zero(T),N))),number_of_particles)

  # Initialize cell, firstatom and nexatom
  for (ip,xip) in pairs(x)
    p = AtomWithIndex(ip,xip)
    icell_cartesian = particle_cell(xip,box)
    icell = cell_linear_index(box.nc,icell_cartesian)
    if fp[icell].index == 0
      ncwp += 1
      cwp[ncwp] = icell
      ncp[icell] = 1
    else
      ncp[icell] += 1
    end
    np[ip] = fp[icell]
    fp[icell] = p
  end

  return CellLists{N,T}(ncwp,cwp,ncp,fp,np)
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
particle_cell(x::AbstractVector{T}, box::Box{N}) where N
```

Returns the coordinates of the cell to which a particle belongs, given its coordinates
and the sides of the periodic box (for arbitrary dimension N).

"""
function particle_cell(x::AbstractVector, box::Box{N}) where N
  # Wrap to origin
  xwrapped = wrapone(x,box.sides)
  cell = SVector{N,Int}(
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
initlists!(x::AbstractVector{<:AbstractVector}, box::Box, lc::LinkedLists)
```

Function that initializes the linked cells by computing to each cell each atom
belongs and filling up the firstatom and nexatom arrays.
Modifies the data in the lc structure

## Example

```julia-repl
julia> lc = LinkedLists(100_000);

julia> box = Box([250,250,250],10);

julia> x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:100_000 ];

julia> initlists!(x,box,lc)

julia> lc.firstatom
125000-element Vector{Int64}:
 51137
 62998
     ⋮
 41753
 53909

```

"""
function initlists!(
  x::AbstractVector{<:AbstractVector}, 
  box::Box, lc::LinkedLists;
  parallel=true
) 

  # Count the number of boxes and checks if there is a problem with dimensions
  nboxes = prod(box.nc)
  if length(lc.firstatom) < nboxes
    resize!(lc.firstatom,nboxes)
  end

  if parallel
    # Reset arrays
    @threads for i in 1:nboxes
      lc.firstatom[i] = 0
    end
    @threads for i in eachindex(lc.nextatom)
      lc.nextatom[i] = 0
    end
    # Initialize cell, firstatom and nexatom
    @threads for iat in eachindex(x)
      icell_cartesian = particle_cell(x[iat],box)
      lc.icell[iat] = cell_linear_index(box.nc,icell_cartesian)
    end
    # This cannot be parallelizad, at least not easily
    for iat in eachindex(x)
      lc.nextatom[iat] = lc.firstatom[lc.icell[iat]]
      lc.firstatom[lc.icell[iat]] = iat
    end
  else
    # Reset arrays
    for i in 1:nboxes
      lc.firstatom[i] = 0
    end
    for i in eachindex(lc.nextatom)
      lc.nextatom[i] = 0
    end
    # Initialize cell, firstatom and nexatom
    for iat in eachindex(x)
      icell_cartesian = particle_cell(x[iat],box)
      icell = cell_linear_index(box.nc,icell_cartesian)
      lc.nextatom[iat] = lc.firstatom[icell]
      lc.firstatom[icell] = iat
      lc.icell[iat] = icell
    end
  end

  return nothing
end

"""

```
map_pairwise!(f::Function,output,x::AbstractVector,box::Box,lc::LinkedLists)
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

Computing the mean difference in `x` position between random particles, remembering the number of pairs of `n` particles is `n(n-1)/2`. The function does not use the indices or the distance, such that we remove them from the parameters by using a closure.

```julia-repl
julia> n = 100_000;

julia> lc = LinkedLists(n);

julia> box = Box([250,250,250],10);

julia> x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:n ];

julia> initlists!(x,box,lc)

julia> f(x,y,sum_dx) = sum_dx + x[1] - y[1] 

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> avg_dx = normalization * map_parwise!((x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, x, box, lc)

```

Computing the histogram of the distances between particles (considering the same particles as in the above example). Again,
the function does not use the indices, but uses the distance, which are removed from the function call using a closure:

```
julia> function build_histogram!(x,y,d2,hist)
         d = sqrt(d2)
         ibin = floor(Int,d) + 1
         hist[ibin] += 1
         return hist
       end;

julia> hist = zeros(Int,10);

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> hist = normalization * map_pairwise!((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),hist,x,box,lc)

```

In this test we compute the "gravitational potential", pretending that each particle
has a different mass. In this case, the closure is used to pass the masses to the
function that computes the potential.

```julia
# masses
mass = rand(N)

# Function to be evalulated for each pair: build distance histogram
function potential(x,y,i,j,d2,u,mass)
  d = sqrt(d2)
  u = u - 9.8*mass[i]*mass[j]/d
  return u
end

# Run pairwise computation
u = map_pairwise!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc)
```

The example above can be run with `CellLists.test3()`.


"""
map_pairwise!(
  f::Function, output, 
  x::AbstractVector,
  box::Box, lc::LinkedLists; 
  # Parallellization options
  parallel::Bool=true,
  output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
  reduce::Function=reduce
) =
  map_pairwise!(f,output,x,x,box,lc; self=true, 
    parallel=parallel,
    output_threaded=output_threaded,
    reduce=reduce
  )

"""

``
map_pairwise!(f::Function, output, 
  x::AbstractVector, y::AbstractVector, 
  box::Box, lc::LinkedLists
)
```

The same as the function `map_pairwise!`, but to compute interactions between two disjoint sets of particles `x` and `y`. 

"""
function map_pairwise!(
  f::Function, output, 
  x::AbstractVector, y::AbstractVector,
  box::Box, lc::LinkedLists; self::Bool=false, 
  # Parallelization options
  parallel::Bool=true,
  output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
  reduce::Function=reduce
)
  if parallel && nthreads() > 1
    output = map_pairwise_parallel!(
      (x,y,i,j,d2,output)->f(x,y,i,j,d2,output),
      output,x,y,box,lc;self=self,
      output_threaded=output_threaded,
      reduce=reduce
    )
  else
    output = map_pairwise_serial!((x,y,i,j,d2,output)->f(x,y,i,j,d2,output),output,x,y,box,lc;self=self)
  end
  return output
end

#
# Parallel version
#
function map_pairwise_parallel!(
  f::Function, output, 
  x::AbstractVector, y::AbstractVector, 
  box::Box, lc::LinkedLists; 
  self=false, 
  output_threaded=output_threaded,
  reduce::Function=reduce
)

  @threads for i in eachindex(x)
    it = threadid()
    output_threaded[it] = inner_map_loop!(
      (x,y,i,j,d2,output) -> f(x,y,i,j,d2,output),
      output_threaded[it],i,x[i],y,box,lc,self
    )
  end

  output = reduce(output,output_threaded)
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
# Serial version
#
function map_pairwise_serial!(
  f::Function, output, 
  x::AbstractVector, y::AbstractVector, 
  box::Box, lc::LinkedLists; 
  self=false
)

  for i in eachindex(x)
    output = inner_map_loop!(
      (x,y,i,j,d2,output) -> f(x,y,i,j,d2,output),
      output,i,x[i],y,box,lc,self
    )
  end
  return output
end

#
# Inner loop that runs over cells for each index of the `x` vector. 
#
function inner_map_loop!(f::Function,output,i,xᵢ,y,box::Box,lc,self)
  @unpack sides, nc, lcell, cutoff_sq = box
  # Check the cell of this atom
  ipc, jpc, kpc = particle_cell(xᵢ,box)
  # Loop over vicinal cells to compute distances to solvent atoms, and
  # add data to dc structure (includes current cell)
  for ic in ipc-lcell:ipc+lcell
    for jc in jpc-lcell:jpc+lcell
      for kc in kpc-lcell:kpc+lcell
        # Wrap cell if needed
        icell_cartesian = wrap_cell(nc,(ic,jc,kc))
        # get linear index of this cell
        icell = cell_linear_index(nc,icell_cartesian) 
        # cycle over the atoms of this cell
        j = lc.firstatom[icell]
        while j > 0
          # skip same particle and repeated
          if self && j >= i 
            j = lc.nextatom[j]
            continue
          end
          # Wrap particle j relative to particle xᵢ
          yⱼ = wrapone(y[j],sides,xᵢ)
          d2 = distance_sq(xᵢ,yⱼ)
          if d2 <= cutoff_sq
            output = f(xᵢ,yⱼ,i,j,d2,output)
          end
          j = lc.nextatom[j]
        end
      end
    end
  end
  return output
end

#
# Function that uses the naive algorithm, for testing
#
function map_naive!(f,output,x,box)
  @unpack sides, cutoff_sq = box
  for i in 1:length(x)-1
    xᵢ = x[i]
    for j in i+1:length(x)
      xⱼ = wrapone(x[j],sides,xᵢ)
      d2 = distance_sq(xᵢ,xⱼ) 
      if d2 <= cutoff_sq
        output = f(xᵢ,xⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

#
# Function that uses the naive algorithm, for testing
#
function map_naive_two!(f,output,x,y,box)
  @unpack sides, cutoff_sq = box
  for i in 1:length(x)
    xᵢ = x[i]
    for j in 1:length(y)
      yⱼ = wrapone(y[j],sides,xᵢ)
      d2 = distance_sq(xᵢ,yⱼ) 
      if d2 <= cutoff_sq
        output = f(xᵢ,yⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

#
# Serial version
#
function map_pairwise_serial!(
  f::Function, output, 
  x::AbstractVector,
  box::Box, cl::CellLists
)

  @unpack sides, nc, lcell, cutoff_sq = box

  # Iterator over all vicinal cells, not including the central one
  neighbour_cells = CartesianIndices((-lcell:lcell, -lcell:lcell, -lcell:lcell))

  # loop over cells that contain particles
  for i in 1:cl.ncwp
    ic = cl.cwp[i]
    ic_cartesian = cell_cartesian_indices(nc,ic)

    # loop over neighbouring cells of this cell
    for jc_cartesian in neighbour_cells
      jc_cartesian_wrapped = wrap_cell(nc,jc_cartesian+ic_cartesian)
      jc = cell_linear_index(nc,jc_cartesian_wrapped)

      # loop over list of particles of cell ic
      pᵢ = cl.fp[ic]
      for _ in 1:cl.ncp[ic]
        i = pᵢ.index
        xpᵢ = pᵢ.coordinates
        pⱼ = cl.fp[jc]
        # loop over particles of cell jc
        for _ in 1:cl.ncp[jc]
          j = pⱼ.index
          if i > j
            xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
            d2 = distance_sq(xpᵢ,xpⱼ)
            if d2 <= cutoff_sq
              output = f(xpᵢ,xpⱼ,i,j,d2,output)
            end
          end
          pⱼ = cl.np[pⱼ.index]
        end
        pᵢ = cl.np[pᵢ.index]
      end

    end # neighbouring cells
  end # cells that contain particles

  return output
end

function map_pairwise_serial_half!(
  f::Function, output, 
  x::AbstractVector,
  box::Box, cl::CellLists
)

  @unpack sides, nc, cutoff_sq = box

  # loop over cells that contain particles
  for i in 1:cl.ncwp
    ic = cl.cwp[i]
    ic_cartesian = cell_cartesian_indices(nc,ic)

    # loop over list of non-repeated particles of cell ic
    pᵢ = cl.fp[ic]
    for incp in 1:cl.ncp[ic]-1
      i = pᵢ.index
      xpᵢ = pᵢ.coordinates
      pⱼ = cl.np[pᵢ.index] 
      for jncp in incp+1:cl.ncp[ic]
        j = pⱼ.index
        xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
        d2 = distance_sq(xpᵢ,xpⱼ)
        if d2 <= cutoff_sq
          output = f(xpᵢ,xpⱼ,i,j,d2,output)
        end
        pⱼ = cl.np[pⱼ.index]
      end
      pᵢ = cl.np[pᵢ.index]
    end

    # cells that share faces
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1, 0, 0)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex(( 0,+1, 0)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex(( 0, 0,+1)))

    # Interactions of cells that share axes
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,+1, 0)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1, 0,+1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,-1, 0)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1, 0,-1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex(( 0,+1,+1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex(( 0,+1,-1)))

    # Interactions of cells that share vertices
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,+1,+1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,+1,-1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,-1,+1)))
    output = cell_output!(f,box,ic,cl,output,ic_cartesian+CartesianIndex((+1,-1,-1)))

  end # cells that contain particles

  return output
end

function cell_output!(f,box,ic,cl,output,jc_cartesian)

  @unpack sides, nc, cutoff_sq = box

  jc_cartesian_wrapped = wrap_cell(nc,jc_cartesian)
  jc = cell_linear_index(nc,jc_cartesian_wrapped)

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[ic]
  for _ in 1:cl.ncp[ic]
    i = pᵢ.index
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.fp[jc]
    for _ in 1:cl.ncp[jc]
      j = pⱼ.index
      xpⱼ = wrapone(pⱼ.coordinates,sides,xpᵢ)
      d2 = distance_sq(xpᵢ,xpⱼ)
      if d2 <= cutoff_sq
        output = f(xpᵢ,xpⱼ,i,j,d2,output)
      end
      pⱼ = cl.np[pⱼ.index]
    end
    pᵢ = cl.np[pᵢ.index]
  end

  return output
end

#
# Test examples
#
include("./examples.jl")

end # module


