module CellListMap

using Base.Threads
using Parameters
using StaticArrays
using DocStringExtensions
using ProgressMeter

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
  unit_cell::SMatrix{N,N,T}
  unit_cell_center::SVector{N,T}
  unit_cell_max::SVector{N,T}
  lcell::Int
  nc::SVector{N,Int}
  cutoff::T
  cutoff_sq::T
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
  unit_cell_center = unit_cell_max / 2 
  unit_cell = SMatrix{N,N,T}(unit_cell) 
  nc = SVector{N,Int}(ceil.(Int,max.(1,unit_cell_max/(cutoff/lcell)))) .+ 1
  return Box{N,T}(
    unit_cell,
    unit_cell_center,
    unit_cell_max,
    lcell,nc,
    cutoff,
    cutoff^2
  )
end
Box(unit_cell::AbstractMatrix,cutoff;T::DataType=Float64,lcell::Int=1) =
  Box(unit_cell,cutoff,T,lcell)

function Base.show(io::IO,::MIME"text/plain",box::Box)
  println(typeof(box))
  println("  unit cell: ", box.unit_cell) 
  println("  unit cell center: ", box.unit_cell_center) 
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
  unit_cell = SMatrix{N,N,T}( 
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
  coordinates::SVector{N,T}
end
AtomWithIndex{N,T}() where {N,T} = AtomWithIndex{N,T}(0,zero(SVector{N,T})) 

"""

$(TYPEDEF)

This structure contains the cell linear index and the information about if this cell
is in the border of the box (such that its neighbouring cells need to be wrapped) 

"""
struct Cell{N}
  icell::Int
  cartesian::CartesianIndex{N}
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef struct CellList{N,T}
  ncwp::Vector{Int} # One-element vector to contain the *mutable* number of cells with particles
  ncp::Vector{Int} # One-element vector to contain the *mutable* number of particles in the computing box
  cwp::Vector{Cell{N}} # Indices of the unique cells with Particles
  fp::Vector{AtomWithIndex{N,T}} # First particle of cell 
  np::Vector{AtomWithIndex{N,T}} # Next particle of cell
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
  println(typeof(cl))
  println("  $(cl.ncwp[1]) cells with particles.")
  println("  $(cl.ncp[1]) particles in computing box, including images.")
  print("  $(length(cl.np)/cl.ncwp[1]) particles per computing cell.")
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
  number_of_cells = ceil(Int,prod(box.nc))
  # next is a lower bound, will be resized when necessary to incorporate particle images
  number_of_particles = length(x)
  ncwp = zeros(Int,1)
  ncp = zeros(Int,1)
  cwp = Vector{Cell{N}}(undef,number_of_cells)
  fp = Vector{AtomWithIndex{N,T}}(undef,number_of_cells)
  np = Vector{AtomWithIndex{N,T}}(undef,number_of_particles)
  cl = CellList{N,T}(ncwp,ncp,cwp,fp,np)
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
      cwp[i] = Cell{N}(0,zero(CartesianIndex{N}))
      fp[i] = AtomWithIndex{N,T}()
    end
    @threads for i in eachindex(np)
      np[i] = AtomWithIndex{N,T}()
    end
  else
    fill!(cwp,Cell{N}(0,zero(CartesianIndex{N})))
    fill!(fp,AtomWithIndex{N,T}())
    fill!(np,AtomWithIndex{N,T}())
  end

  #
  # Not worth paralellizing probably (would need to take care of concurrency)
  #
  for particle in x
    add_images_to_celllist!(particle,box,cl)
  end

  return cl
end

# Function that checks if the particule is outside the computation bounding box
function out_of_bounding_box(x::SVector{N,T},box::Box{N,T}) where {N,T}
  for i in 1:N
    (x[i] >= box.nc[i]*box.cutoff) && return true
    (x[i] < 0) && return true
  end
  return false
end
out_of_bounding_box(p::AtomWithIndex,box::Box{N,T}) where {N,T} =
  out_of_bounding_box(p.coordinates,box)

function add_images_to_celllist!(particle,box::Box{N,T},cl) where {N,T}
  # Wrap to first image with positive coordinates
  p = wrap_to_first(particle,box.unit_cell)
  # current particle position
  add_particle_to_celllist!(p,box,cl) 
  # images that fall within the computing box
  steps = Iterators.filter(
    step -> count(isequal(0),step) != N,
    Iterators.product(
      ntuple(i -> -1:1, N)...
    ) 
  )
  for step in steps
    indices = ntuple( i -> step[i], N) 
    x = translation_image(p,box.unit_cell,indices)
    while ! out_of_bounding_box(x,box)
      add_particle_to_celllist!(x,box,cl) 
      indices = ntuple( i -> indices[i] + step[i], N) 
      x = translation_image(p,box.unit_cell,indices)
    end
  end
  return nothing
end

"""

```
view_celllist_particles(cl::CellList)
```

Auxiliary function to view the particles of a computing box, including images created
for computing purposes.

### Example
```julia
julia> box = Box([ 100 50; 50 100 ],10);

julia> p = [ box.unit_cell_max .* rand(SVector{2,Float64}) for i in 1:1000 ];

julia> cl = CellList(p,box);

julia> x, y = CellListMap.view_celllist_particles(cl);

julia> using Plots

julia> scatter(x,y,label=nothing,xlims=(-10,180),ylims=(-10,180));

```

"""
function view_celllist_particles(cl::CellList{N,T}) where {N,T}
  @unpack ncwp, cwp, ncp, fp, np = cl
  x = Vector{SVector{N,T}}(undef,ncp[1])
  ip = 0
  for i in 1:ncwp[1]
    ip += 1
    p = fp[cwp[i].icell]
    x[ip] = p.coordinates
    while np[p.index].index > 0
      ip += 1
      x[ip] = p.coordinates
      p = np[p.index]
    end
  end
  return ([x[i][j] for i in 1:ncp[1]] for j in 1:N)
end

"""

Set one index of a cell list

"""
function add_particle_to_celllist!(particle::SVector{N,T},box,cl) where {N,T}
  @unpack ncwp, ncp, cwp, fp, np = cl
  ncp[1] += 1
  p = AtomWithIndex(cl.ncp[1],particle)
  icell_cartesian = particle_cell(p.coordinates,box)
  icell = cell_linear_index(box.nc,icell_cartesian)
  if fp[icell].index == 0
    ncwp[1] += 1
    cwp[ncwp[1]] = Cell{N}(icell,icell_cartesian)
  end
  if ncp[1] > length(np) 
    old_length = length(np)
    resize!(np,ceil(Int,1.2*old_length))
    for i in old_length+1:length(np)
      np[i] = AtomWithIndex{N,T}() 
    end
  end
  np[ncp[1]] = fp[icell]
  fp[icell] = p
  return nothing
end

"""

wrap_cell_fraction(x::SVector{N,T}, cell)

"""
function wrap_cell_fraction(x::T, cell) where T
  p = MVector{length(x),eltype(x)}(rem.(cell\x,1))
  for i in eachindex(p)
    if p[i] < 0
      p[i] += 1
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
function wrap_to_first(x::T,cell) where T<:AbstractVector
  p = wrap_cell_fraction(x,cell)
  return T(cell*p)
end

"""

```
wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

"""
function wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
  return T(rem.(cell\(x-xref),1))
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

neighbour_cells_all(lcell) = 
  CartesianIndices((-lcell:lcell,-lcell:lcell,-lcell:lcell))


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
  CartesianIndex(ntuple(i -> floor(Int,x[i]/(box.cutoff/box.lcell) + 1), N))

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

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellList; 
  show_progress::Bool=false
) where {F}
  show_progress && (p = Progress(cl.ncwp[1],dt=1))
  for icell in 1:cl.ncwp[1]
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
  show_progress && (p = Progress(cl.ncwp[1],dt=1))
  @threads for it in 1:nthreads() 
    for icell in splitter(it,cl.ncwp[1])
      output_threaded[it] = inner_loop!(f,box,icell,cl,output_threaded[it]) 
      show_progress && next!(p)
    end
  end 
  output = reduce(output,output_threaded)
  return output
end

function inner_loop!(f,box,icell,cl::CellList,output)
  @unpack sides, cutoff_sq = box
  cell = cl.cwp[icell]

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[cell.icell]
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
    output = cell_output!(f,box,cell,cl,output,cell.cartesian+jcell)
  end

  return output
end

#
# loops over the particles of a neighbour cell
#
function cell_output!(f,box,cell,cl,output,jc_cartesian)
  @unpack sides, nc, cutoff_sq = box
  if cell.inborder
    jc_cartesian_wrapped = wrap_cell(nc,jc_cartesian)
  else
    jc_cartesian_wrapped = jc_cartesian
  end
  jc = cell_linear_index(nc,jc_cartesian_wrapped)

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[cell.icell]
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
function map_pairwise_serial!(f::F, output, box::Box, cl::CellListPair; 
  parallel::Bool=false,
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
  parallel::Bool=false,
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
  @unpack sides, nc, lcell, cutoff_sq = box
  xpᵢ = wrapone(cl.small[i],box.sides)
  ic = particle_cell(xpᵢ,box)
  inborder = cell_in_border(ic,box)
  for neighbour_cell in neighbour_cells_all(lcell)
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



