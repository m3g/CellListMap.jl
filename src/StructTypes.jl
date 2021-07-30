
#
# Type of systems, to dispatch methods accordingly
#
struct TinySystem end 
struct LowDensitySystem end
struct HighDensitySystem end

#
# We make difference between these two because wrapping in orthorhombic 
# cells is cheaper
#
struct TriclinicCell end
struct OrthorhombicCell end

struct UnitCell{UnitCellType,N,T,M}
  matrix::SMatrix{N,N,T,M}
end

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
Base.@kwdef struct Box{UnitCellType,N,T,M}
  unit_cell::UnitCell{UnitCellType,N,T,M}
  lcell::Int
  nc::SVector{N,Int}
  cutoff::T
  cutoff_sq::T
  ranges::SVector{N,UnitRange{Int}}
  scale_cutoff::T
  cell_size::T
  unit_cell_max::SVector{N,T}
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
function Box(
  unit_cell_matrix::AbstractMatrix, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  scale_cutoff=1.0,
  UnitCellType=TriclinicCell
)

  @assert lcell >= 1 "lcell must be greater or equal to 1"
  @assert scale_cutoff >= 1 "scale_cutoff must be greater or equal to 1.0"

  N = size(unit_cell_matrix)[1]
  @assert N == size(unit_cell_matrix)[2] "Unit cell matrix must be square."
  @assert count(unit_cell_matrix .< 0) == 0 "Unit cell lattice vectors must only contain non-negative coordinates."

  unit_cell = UnitCell{UnitCellType,N,T,N*N}(SMatrix{N,N,T,N*N}(unit_cell_matrix))
  cell_size = scale_cutoff*cutoff/lcell
  unit_cell_max = sum(@view(unit_cell_matrix[:,i]) for i in 1:N) 

### lcell*cell_size must be == unit_cell_max

  nc = SVector{N,Int}(ceil.(Int,unit_cell_max/cell_size) .+ 2)

  ranges = ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix)
  return Box{UnitCellType,N,T,N*N}(
    unit_cell,
    lcell, 
    nc,
    cutoff,
    cutoff^2,
    ranges,
    scale_cutoff,
    cell_size,
    unit_cell_max
  )
end
Box(
  unit_cell_matrix::AbstractMatrix,
  cutoff;
  T::DataType=Float64,
  lcell::Int=1,
  scale_cutoff=1.0,
  UnitCellType=TriclinicCell
) = Box(unit_cell_matrix,cutoff,T,lcell,scale_cutoff,UnitCellType)

function Base.show(io::IO,::MIME"text/plain",box::Box)
  println(typeof(box))
  println("  unit cell matrix: ", box.unit_cell.matrix) 
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
function Box(
  sides::AbstractVector, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  scale_cutoff=1.0,
  UnitCellType=OrthorhombicCell
)
  N = length(sides)
  cart_idxs = CartesianIndices((1:N,1:N))
  # Build unit cell matrix from lengths
  unit_cell_matrix = SMatrix{N,N,T,N*N}( 
    ntuple(N*N) do i
      c = cart_idxs[i]
      if c[1] == c[2] 
        return sides[c[1]] 
      else
        return zero(T)
      end
    end
  )
  return Box(
    unit_cell_matrix,
    cutoff,
    T,
    lcell,
    scale_cutoff,
    UnitCellType
  ) 
end
Box(
  sides::AbstractVector,
  cutoff;
  T::DataType=Float64,
  lcell::Int=1,
  scale_cutoff=1.0,
  UnitCellType=OrthorhombicCell
) = Box(sides,cutoff,T,lcell,scale_cutoff,UnitCellType)

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

This structure contains the cell linear index and the information 
about if this cell is in the border of the box (such that its 
neighbouring cells need to be wrapped) 

"""
struct Cell{N,T}
  icell::Int
  cartesian::CartesianIndex{N}
  center::SVector{N,T}
end
Base.zero(::Type{Cell{N,T}}) where {N,T} =
  Cell{N,T}(0,CartesianIndex{N}(ntuple(i->0,N)),zeros(SVector{N,T}))

"""

Auxiliary structure to contain projected particles in large and
and dense systems.

"""
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
Base.@kwdef struct CellList{SystemType,N,T}
  " *mutable* number of cells with real particles. "
  ncwp::Vector{Int}
  " *mutable* number of particles in the computing box "
  ncp::Vector{Int}
  " Auxiliary array to annotate if the cell contains real particles. "
  contains_real::Vector{Bool}
  " Indices of the unique cells with real particles. "
  cwp::Vector{Cell{N,T}}
  " First particle of cell. "
  fp::Vector{AtomWithIndex{N,T}}
  " Next particle of cell "
  np::Vector{AtomWithIndex{N,T}}
  " Number of particles of cell. "
  npcell::Vector{Int}
  " Auxiliar array to store projected particles. "
  projected_particles::Vector{Vector{ProjectedParticle{N,T}}}
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
  println(typeof(cl))
  println("  $(cl.ncwp[1]) cells with real particles.")
  println("  $(cl.ncp[1]) particles in computing box, including images.")
end

# Structure that will cointain the cell lists of two independent sets of
# particles for cross-computation of interactions
@with_kw struct CellListPair{SystemType,V,N,T}
  small::V
  large::CellList{SystemType,N,T}
  swap::Bool
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
function CellList(
  x::AbstractVector{SVector{N,T}},
  box::Box;
  parallel::Bool=true,
  SystemType=nothing
) where {N,T} 
  number_of_cells = prod(box.nc)
  # number_of_particles is a lower bound, will be resized when necessary to incorporate particle images
  number_of_particles = ceil(Int,1.2*length(x))
  ncwp = [0]
  ncp = [0]
  cwp = Vector{Cell{N,T}}(undef,number_of_cells)
  contains_real = Vector{Bool}(undef,number_of_cells)
  fp = Vector{AtomWithIndex{N,T}}(undef,number_of_cells)
  np = Vector{AtomWithIndex{N,T}}(undef,number_of_particles)
  npcell = Vector{Int}(undef,number_of_cells)
  projected_particles = [ Vector{ProjectedParticle{N,T}}(undef,0) for i in 1:nthreads() ]

  # Set automatically the system type, as a function of the 
  # total number of particles and the number of particles per cell
  if isnothing(SystemType) 
    SystemType = set_system_type(x,box)
  end

  cl = CellList{SystemType,N,T}(
    ncwp,
    ncp,
    contains_real,
    cwp,
    fp,
    np,
    npcell,
    projected_particles
  )
  return UpdateCellList!(x,box,cl,parallel=parallel)
end

function set_system_type(x,box)

  if length(x) <= 20 
    return TinySystem
  end

  particles_per_cell = length(x) / prod(box.nc)

  # what is low? what is high? to test
  if particles_per_cell < 32
    return LowDensitySystem
  end

  return HighDensitySystem
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

particles_per_cell(cl::CellList,box::Box) = cl.ncp[1] / prod(box.nc)

