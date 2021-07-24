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
Base.@kwdef struct Box{N,T,M}
  unit_cell::SMatrix{N,N,T,M}
  unit_cell_max::SVector{N,T}
  lcell::Int
  nc::SVector{N,Int}
  cutoff::T
  cutoff_sq::T
  ranges::SVector{N,UnitRange{Int}}
end

"""

```
Box(unit_cell::AbstractMatrix, cutoff; T::DataType=Float64, lcell::Int=1)
```

Construct box structure given the cell matrix of lattice vectors. 

### Example
```julia
julia> unit_cell = [ 100   50    0 
                       0  120    0
                       0    0  130 ];

julia> box = Box(unit_cell,10)
Box{3, Float64}
  unit cell: [100.0 50.0 0.0; 0.0 120.0 0.0; 0.0 0.0 130.0]
  cutoff: 10.0
  number of cells on each dimension: [16, 13, 14] (lcell: 1)
  Total number of cells: 2912

```

"""
function Box(unit_cell::AbstractMatrix, cutoff, T::DataType, lcell::Int=1)
  N = size(unit_cell)[1]
  @assert N == size(unit_cell)[2] "Unit cell matrix must be square."
  @assert count(unit_cell .< 0) == 0 "Unit cell lattice vectors must only contain non-negative coordinates."
  unit_cell_max = SVector{N,T}(sum(@view(unit_cell[:,i]) for i in 1:N))
  unit_cell = SMatrix{N,N,T}(unit_cell) 
  nc = SVector{N,Int}(
    ceil.(Int,max.(1,(unit_cell_max .+ 2*cutoff)/(cutoff/lcell)))
  ) 
  ranges = ranges_of_replicas(cutoff,nc,unit_cell,unit_cell_max)
  return Box{N,T,N*N}(
    unit_cell,
    unit_cell_max,
    lcell, nc,
    cutoff,
    cutoff^2,
    ranges
  )
end
Box(unit_cell::AbstractMatrix,cutoff;T::DataType=Float64,lcell::Int=1) =
  Box(unit_cell,cutoff,T,lcell)

function Base.show(io::IO,::MIME"text/plain",box::Box)
  println(typeof(box))
  println("  unit cell: ", box.unit_cell) 
  println("  unit cell maximum: ", box.unit_cell_max) 
  println("  cutoff: ", box.cutoff)
  println("  number of cells on each dimension: ",box.nc, " (lcell: ",box.lcell,")")
  print("  Total number of cells: ", prod(box.nc))
end

"""

```
Box(sides::AbstractVector, cutoff; T::DataType=Float64, lcell::Int=1)
```

For orthorhombic unit cells, `Box` can be initialized with a vector of the length of each side. 

### Example
```julia
julia> box = Box([120,150,100],10)
Box{3, Float64}
  unit cell: [120.0 0.0 0.0; 0.0 150.0 0.0; 0.0 0.0 100.0]
  cutoff: 10.0
  number of cells on each dimension: [13, 16, 11] (lcell: 1)
  Total number of cells: 2288

```

"""
function Box(sides::AbstractVector, cutoff, T::DataType, lcell::Int=1)
  N = length(sides)
  cart_idxs = CartesianIndices((1:N,1:N))
  # Build unit cell matrix from lengths
  unit_cell = SMatrix{N,N,T,N*N}( 
    ntuple(N*N) do i
      c = cart_idxs[i]
      if c[1] == c[2] 
        return sides[c[1]] 
      else
        return zero(T)
      end
    end
  )
  return Box(unit_cell,cutoff,T,lcell) 
end
Box(sides::AbstractVector,cutoff;T::DataType=Float64,lcell::Int=1) =
  Box(sides,cutoff,T,lcell)

"""

$(TYPEDEF)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort.

"""
struct AtomWithIndex{N,T}
  index::Int
  index_original::Int
  coordinates::SVector{N,T}
end
Base.zero(::Type{AtomWithIndex{N,T}}) where {N,T} =
  AtomWithIndex{N,T}(0,0,zeros(SVector{N,T}))

"""

$(TYPEDEF)

This structure contains the cell linear index and the information about if this cell
is in the border of the box (such that its neighbouring cells need to be wrapped) 

"""
struct Cell{N,T}
  icell::Int
  cartesian::CartesianIndex{N}
  center::SVector{N,T}
end
Base.zero(::Type{Cell{N,T}}) where {N,T} =
  Cell{N,T}(0,CartesianIndex{N}(0,0,0),zeros(SVector{N,T}))

"""

```
cell_center(c::CartesianIndex{N},cutoff::T) where {N,T}
```

Computes the geometric center of a cell, to be used in the projection
of points. Returns a `SVector{N,T}`

"""
cell_center(c::CartesianIndex{N},cutoff::T) where {N,T} =
  SVector{N,T}(ntuple(i-> cutoff*c[i]/2, N))

struct ProjectedParticle{N,T}
  index_original::Int
  xproj::T
  coordinates::SVector{N,T}
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef struct CellList{N,T}
  " *mutable* number of cells with particles "
  ncwp::Int
  " *mutable* number of particles in the computing box "
  ncp::Int
  " Indices of the unique cells with Particles "
  cwp::Vector{Cell{N,T}}
  " First particle of cell "
  fp::Vector{AtomWithIndex{N,T}}
  " Next particle of cell "
  np::Vector{AtomWithIndex{N,T}}
  " Number of particles in of cell "
  npcell::Vector{Int}
  " Auxiliar vector to store projected particles. "
  projected_particles::Vector{ProjectedParticle{N,T}}
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
  println(typeof(cl))
  println("  $(cl.ncwp) cells with real particles.")
  print("  $(cl.ncp) particles in computing box, including images.")
end

# Structure that will cointain the cell lists of two independent sets of
# particles for cross-computation of interactions
@with_kw struct CellListPair{V,N,T}
  small::V
  large::CellList{N,T}
  swap::Bool
end      
function Base.show(io::IO,::MIME"text/plain",cl::CellListPair)
  print(typeof(cl),"\n")
  print("   $(length(cl.small)) particles in the smallest vector.\n")
  print("   $(cl.large.ncwp) cells with particles.")
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
  number_of_cells = ceil(Int,prod(box.nc))
  # number_of_particles is a lower bound, will be resized when necessary to incorporate particle images
  number_of_particles = length(x)
  ncwp = 0
  ncp = 0
  cwp = Vector{Cell{N,T}}(undef,number_of_cells)
  fp = Vector{AtomWithIndex{N,T}}(undef,number_of_cells)
  np = Vector{AtomWithIndex{N,T}}(undef,number_of_particles)
  npcell = Vector{Int}(undef,number_of_cells)
  projected_particles = Vector{ProjectedParticle{N,T}}(undef,number_of_particles)
  cl = CellList{N,T}(ncwp,ncp,cwp,fp,np,npcell,projected_particles)
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
CellListMap.CellListPair{3, Float64}
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
    cl_pair = CellListPair(small=x,large=y_cl,swap=false)
  else
    x_cl = CellList(x,box,parallel=parallel)
    cl_pair = CellListPair(small=y,large=x_cl,swap=true)
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
  @unpack cwp, fp, np, npcell = cl

  number_of_cells = prod(box.nc)
  if number_of_cells > length(cwp) 
    number_of_cells = ceil(Int,1.1*number_of_cells) # some margin in case of box size variations
    resize!(cwp,number_of_cells)
    resize!(fp,number_of_cells)
    resize!(npcell,number_of_cells)
  end

  @set! cl.ncwp = 0
  if parallel
    @threads for i in eachindex(cwp)
      cwp[i] = zero(Cell{N,T})
      fp[i] = zero(AtomWithIndex{N,T})
      npcell[i] = 0
    end
    @threads for i in eachindex(np)
      np[i] = zero(AtomWithIndex{N,T})
    end
  else
    fill!(cwp,zero(Cell{N,T}))
    fill!(fp,zero(AtomWithIndex{N,T}))
    fill!(np,zero(AtomWithIndex{N,T}))
    fill!(npcell,0)
  end

  #
  # The following part cannot be *easily* paralelized, because 
  # there is concurrency on the construction of the cell lists
  #

  #
  # Add virtual particles to edge cells
  #
  for (ip,particle) in pairs(x)
    p = wrap_to_first(particle,box.unit_cell)
    cl = replicate_particle!(ip,p,box,cl)
  end
  #
  # Add true particles, such that the first particle of each cell is
  # always a true particle
  #
  for (ip,particle) in pairs(x)
    p = wrap_to_first(particle,box.unit_cell)
    cl = add_particle_to_celllist!(ip,p,box,cl) 
  end

  return cl
end

# Function that checks if the particule is outside the computation bounding box
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

Set one index of a cell list

"""
function add_particle_to_celllist!(ip,x::SVector{N,T},box,cl;real_particle::Bool=true) where {N,T}
  @unpack cutoff = box
  @unpack cwp, fp, np, npcell = cl
  @set! cl.ncp += 1
  icell_cartesian = particle_cell(x,box)
  icell = cell_linear_index(box.nc,icell_cartesian)
  # Cells starting with real particles are annotated to be run over
  if fp[icell].index == 0
    npcell[icell] = 1
    if real_particle 
      @set! cl.ncwp += 1
      cwp[cl.ncwp] = Cell{N,T}(icell,icell_cartesian,cell_center(icell_cartesian,cutoff))
    end
  else
    npcell[icell] += 1
  end
  if cl.ncp > length(np) 
    old_length = length(np)
    resize!(np,ceil(Int,1.2*old_length))
    for i in old_length+1:length(np)
      np[i] = zero(AtomWithIndex{N,T}) 
    end
  end
  np[cl.ncp] = fp[icell]
  fp[icell] = AtomWithIndex(cl.ncp,ip,x) 
  return cl
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

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellList; 
  show_progress::Bool=false
) where {F}
  show_progress && (p = Progress(cl.ncwp,dt=1))
  for icell in 1:cl.ncwp
    output = inner_loop!(f,box,icell,cl,output) 
    show_progress && next!(p)
  end
  return output
end

# Chunck splitter, we jump over cells such that possible 
# heterogeneities have greater changes of getting split into different chuncks
splitter(it,n) = it:nthreads():n

#
# Parallel version for self-pairwise computations
#
function map_pairwise_parallel!(f::F1, output, box::Box, cl::CellList;
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2}
  show_progress && (p = Progress(cl.ncwp,dt=1))
  @threads for it in 1:nthreads() 
    for icell in splitter(it,cl.ncwp)
      output_threaded[it] = inner_loop!(f,box,icell,cl,output_threaded[it]) 
      show_progress && next!(p)
    end
  end 
  output = reduce(output,output_threaded)
  return output
end

function inner_loop!(f,box,icell,cl::CellList,output)
  @unpack cutoff_sq = box
  cell = cl.cwp[icell]

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[cell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.np[i] 
    j = pⱼ.index
    while j > 0
      xpⱼ = pⱼ.coordinates
      d2 = distance_sq(xpᵢ,xpⱼ)
      if d2 <= cutoff_sq
        i_orig = pᵢ.index_original
        j_orig = pⱼ.index_original
        output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      end
      pⱼ = cl.np[pⱼ.index]
      j = pⱼ.index
    end
    pᵢ = cl.np[pᵢ.index]
    i = pᵢ.index
  end

  for jcell in neighbour_cells(box)
    output = cell_output!(f,box,cell,cl,output,cell.cartesian+jcell)
  end

  return output
end

"""

```
partialsort_cutoff!(x,cutoff)
```

Function that reorders x vector by putting in the first positions the
elements with values smaller than cutoff. Returns the number of elements
that satisfy the condition.

"""
function partialsort_cutoff!(x,cutoff)
  iswap = 1
  @inbounds for i in eachindex(x)
    if x[i].xproj <= cutoff
      if iswap != i
        x[iswap], x[i] = x[i], x[iswap]
      end
      iswap += 1
    end
  end
  return iswap - 1
end

#
# loops over the particles of a neighbour cell
#
function cell_output!(f,box,cell,cl,output,jc_cartesian)
  @unpack projected_particles = cl
  @unpack nc, cutoff, cutoff_sq = box
  jc = cell_linear_index(nc,jc_cartesian)

  # Copy coordinates of particles of cell jcell into continuous array
  pⱼ = cl.fp[jc]
  npcell = cl.npcell[jc]
  j = pⱼ.index
  for jp in 1:npcell
    j_orig = pⱼ.index_original
    xpⱼ = pⱼ.coordinates
    projected_particles[jp] = ProjectedParticle(j_orig,0.,xpⱼ) 
    pⱼ = cl.np[j]
    j = pⱼ.index
  end
  pp = @view(projected_particles[1:npcell])

  # Loop over particles of cell icell
  pᵢ = cl.fp[cell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates

    @turbo for j in 1:npcell 
      j_orig = pp[j].index_original
      xpⱼ = pp[j].coordinates
      d2 = distance_sq(xpᵢ,xpⱼ)
      projected_particles[j] = ProjectedParticle(j_orig,d2,xpⱼ) 
    end
    n = partialsort_cutoff!(pp, cutoff_sq)

    for j in 1:n
      xpⱼ = pp[j].coordinates
      d2 = pp[j].xproj
      i_orig = pᵢ.index_original
      j_orig = pp[j].index_original
      output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
    end

    pᵢ = cl.np[i]
    i = pᵢ.index
  end

  return output
end

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(f::F, output, box::Box, cl::CellListPair; 
  show_progress=show_progress
) where {F}
  show_progress && (p = Progress(length(cl.small),dt=1))
  for i in eachindex(cl.small)
    output = inner_loop!(f,output,i,box,cl)
    show_progress && next!(p)
  end
  return output
end

#
# Parallel version for cross-interaction computations
#
function map_pairwise_parallel!(f::F1, output, box::Box, cl::CellListPair;
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress=show_progress
) where {F1,F2}
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
function inner_loop!(f,output,i,box,cl::CellListPair)
  @unpack unit_cell, nc, cutoff_sq = box
  xpᵢ = wrap_to_first(cl.small[i],unit_cell)
  ic = particle_cell(xpᵢ,box)
  for neighbour_cell in neighbour_cells_all(box)
    jc = cell_linear_index(nc,neighbour_cell+ic)
    pⱼ = cl.large.fp[jc]
    j = pⱼ.index
    # loop over particles of cell jc
    while j > 0
      xpⱼ = pⱼ.coordinates
      d2 = distance_sq(xpᵢ,xpⱼ)
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



