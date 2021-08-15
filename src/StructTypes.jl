#
# This difference is important because for Orthorhombic cells it is
# possible to run over only half of the cells, and wrapping coordinates
# in Orthorhombic cells is slightly cheaper. 
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

## Examples

```julia-repl
julia> sides = [250,250,250];

julia> cutoff = 10;

julia> box = Box(sides,cutoff)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [250.0 0.0 0.0; 0.0 250.0 0.0; 0.0 0.0 250.0]
  cutoff: 10.0
  number of computing cells on each dimension: [27, 27, 27]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 19683

julia> box = Box([ 10  0  0 
                    0 10  5
                    0  0 10 ], 1)
Box{TriclinicCell, 3, Float64, 9}
unit cell matrix: [10.0 0.0 0.0; 0.0 10.0 5.0; 0.0 0.0 10.0]
cutoff: 1.0
number of computing cells on each dimension: [12, 17, 12]
computing cell sizes: [1.0, 1.0, 1.0] (lcell: 1)
Total number of cells: 2448

```

"""
Base.@kwdef struct Box{UnitCellType,N,T,M}
    unit_cell::UnitCell{UnitCellType,N,T,M}
    lcell::Int
    nc::SVector{N,Int}
    cutoff::T
    cutoff_sq::T
    ranges::SVector{N,UnitRange{Int}}
    cell_size::SVector{N,T}
    unit_cell_max::SVector{N,T}
end

"""

```
Box(
  unit_cell_matrix::AbstractMatrix, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  UnitCellType=TriclinicCell
)
```

Construct box structure given the cell matrix of lattice vectors. This 
constructor will always return a `TriclinicCell` box type, unless the
`UnitCellType` parameter is set manually to `OrthorhombicCell`

### Example
```julia
julia> unit_cell = [ 100   50    0 
                       0  120    0
                       0    0  130 ];

julia> box = Box(unit_cell,10)
Box{TriclinicCell, 3, Float64, 9}
  unit cell matrix: [100.0 50.0 0.0; 0.0 120.0 0.0; 0.0 0.0 130.0]
  cutoff: 10.0
  number of computing cells on each dimension: [17, 14, 15]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 3570

```

"""
function Box(
    unit_cell_matrix::AbstractMatrix, 
    cutoff, 
    T::DataType, 
    lcell::Int=1,
    UnitCellType=TriclinicCell
)

    s = size(unit_cell_matrix)
    unit_cell_matrix = SMatrix{s[1],s[2],Float64,s[1]*s[2]}(unit_cell_matrix)

    @assert lcell >= 1 "lcell must be greater or equal to 1"

    N = size(unit_cell_matrix)[1]
    @assert N == size(unit_cell_matrix)[2] "Unit cell matrix must be square."
    @assert check_unit_cell(unit_cell_matrix,cutoff) " Unit cell matrix does not satisfy required conditions."

    unit_cell = UnitCell{UnitCellType,N,T,N*N}(SMatrix{N,N,T,N*N}(unit_cell_matrix))
    unit_cell_max = sum(@view(unit_cell_matrix[:,i]) for i in 1:N) 
 
    nc_min = @. floor.(Int,unit_cell_max/cutoff) 
    cell_size_max = unit_cell_max ./ nc_min
    cell_size = cell_size_max/lcell

    nc = SVector{N,Int}(ceil.(Int,unit_cell_max ./ cell_size) .+ 2*lcell)

    #ranges = ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix)
    ranges = SVector{N,UnitRange{Int}}(ntuple(i->-1:1,N))
    return Box{UnitCellType,N,T,N*N}(
        unit_cell,
        lcell, 
        nc,
        cutoff,
        cutoff^2,
        ranges,
        cell_size,
        unit_cell_max
  )
end
Box(
    unit_cell_matrix::AbstractMatrix,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1,
    UnitCellType=TriclinicCell
) = Box(unit_cell_matrix,cutoff,T,lcell,UnitCellType)

function Base.show(io::IO,::MIME"text/plain",box::Box)
    println(typeof(box))
    println("  unit cell matrix: ", box.unit_cell.matrix) 
    println("  cutoff: ", box.cutoff)
    println("  number of computing cells on each dimension: ",box.nc)
    println("  computing cell sizes: ", box.cell_size, " (lcell: ",box.lcell,")")
    print("  Total number of cells: ", prod(box.nc))
end

"""

```
Box(
  sides::AbstractVector, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  UnitCellType=OrthorhombicCell
)
```

For orthorhombic unit cells, `Box` can be initialized with a vector of the 
length of each side. 

### Example
```julia
julia> box = Box([120,150,100],10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [120.0 0.0 0.0; 0.0 150.0 0.0; 0.0 0.0 100.0]
  cutoff: 10.0
  number of computing cells on each dimension: [14, 17, 12]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 2856

```

"""
function Box(
    sides::AbstractVector, 
    cutoff, 
    T::DataType, 
    lcell::Int=1,
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
        UnitCellType
    ) 
end
Box(
    sides::AbstractVector,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1,
    UnitCellType=OrthorhombicCell
) = Box(sides,cutoff,T,lcell,UnitCellType)

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort.

"""
struct ParticleWithIndex{N,T}
    index_original::Int
    coordinates::SVector{N,T}
    real::Bool
end
Base.zero(::Type{ParticleWithIndex{N,T}}) where {N,T} =
    ParticleWithIndex{N,T}(0,zeros(SVector{N,T}),false)

"""

$(TYPEDEF)

$(TYPEDFIELDS)

This structure contains the cell linear index and the information 
about if this cell is in the border of the box (such that its 
neighbouring cells need to be wrapped) 

"""
Base.@kwdef mutable struct Cell{N,T}
    linear_index::Int = 0
    cartesian_index::CartesianIndex{N} = CartesianIndex{N}(ntuple(i->0,N))
    center::SVector{N,T} = zeros{SVector{N,T}}
    contains_real_particles::Bool = false
    number_of_particles::Int = 0
    particles::Vector{ParticleWithIndex{N,T}} = Vector{ParticleWithIndex{N,T}}(undef,0)
end
Base.zero(::Type{Cell{N,T}}) where {N,T} = Cell{N,T}()
Cell(particles::Vector{ParticleWithIndex{N,T}}) where {N,T} = 
    Cell{N,T}(particles=particles)

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Auxiliary structure to contain projected particles in large and
and dense systems.

"""
Base.@kwdef struct ProjectedParticle{N,T}
    index_original::Int = 0
    xproj::T = zero(T)
    coordinates::SVector{N,T} = zeros(SVector{N,T})
    real::Bool = false
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef mutable struct CellList{N,T}
    " *mutable* number of cells with real particles. "
    number_of_cells_with_particles::Int
    " *mutable* number of particles in the computing box. "
    number_of_particles::Vector{Int}
    " Auxiliary array that contains the indexes in cwp of the cells with real particles. "
    cell_index_in_list::Vector{Int}
    " Vector containing cell lists of cells with particles. "
    list::Vector{Cell{N,T}}
    " Auxiliar array to store projected particles. "
    projected_particles::Vector{Vector{ProjectedParticle{N,T}}}
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
    println(typeof(cl))
    println("  $(cl.number_of_cells_with_particles) cells with real particles.")
    println("  $(cl.number_of_particles) particles in computing box, including images.")
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that will cointain the cell lists of two independent sets of
particles for cross-computation of interactions

"""
@with_kw struct CellListPair{V,N,T}
    small::V
    large::CellList{N,T}
    swap::Bool
end      
function Base.show(io::IO,::MIME"text/plain",cl::CellListPair)
    print(typeof(cl),"\n")
    print("   $(length(cl.small)) particles in the smallest vector.\n")
    print("   $(cl.large.number_of_cells_with_particles) cells with particles.")
end
  
"""

```
CellList(
    x::AbstractVector{AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
```

Function that will initialize a `CellList` structure from scracth, given a vector
or particle coordinates (a vector of vectors, typically of static vectors) 
and a `Box`, which contain the size ofthe system, cutoff, etc.  

### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:100000 ];

julia> cl = CellListMap.CellList(x,box)
CellList{3, Float64}
  1800 cells with real particles.
  161084 particles in computing box, including images.

```

"""
function CellList(
    x::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
    cl = init_cell_list(x,box)
    return UpdateCellList!(x,box,cl,parallel=parallel)
end

function init_cell_list(x,box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}
    number_of_cells = prod(box.nc)
    # number_of_particles is a lower bound, will be resized when necessary to incorporate 
    # particle images
    upper_bound_particles = ceil(Int,1.2*length(x))
    number_of_cells_with_particles = 0
    number_of_particles = 0
    particles_per_cell = upper_bound_particles ÷ number_of_cells
    list = 
        [ Cell(Vector{ParticleWithIndex{N,T}}(undef,particles_per_cell)) for i in 1:number_of_cells ] 
    projected_particles = 
        [ Vector{ProjectedParticle{N,T}}(undef,0) for _ in 1:nthreads() ]

  cl = CellList{N,T}(
      number_of_cells_with_particles,
      number_of_particles,
      list,
      projected_particles
  )
  return cl
end

function reset!(cwp::Vector{Cell{N,T}}) where {N,T}
    for i in eachindex(cwp)
        cwp[i] = Cell{N,T}(cwp[i].particles)
    end
    return cwp
end

function reset!(cl::CellList{N,T},box) where{N,T}
    number_of_cells = prod(box.nc)
    if number_of_cells > length(cwp) 
        number_of_cells = ceil(Int,1.2*number_of_cells) # some margin in case of box size variations
        resize!(cwp,number_of_cells)
    end
    cl.number_of_cells_with_particles = 0
    reset!(cl.cwp)
    return cl
end

"""

```
CellList(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
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
   1767 cells with particles.

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
particles_per_cell(cl::CellList,box::Box)
```

Returns the average number of particles per computing cell.

"""
particles_per_cell(cl::CellList,box::Box) = cl.ncp[1] / prod(box.nc)

"""

```
UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    box::Box,cl:CellList{N,T},
    parallel=true
) where {N,T}
```

Function that will update a previously allocated `CellList` structure, given new 
updated particle positions.

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
    x::AbstractVector{<:AbstractVector},
    box::Box,
    cl::CellList{N,T};
    parallel::Bool=true
) where {N,T}
    @unpack cwp, projected_particles = cl

    # Reset cell (resize if needed, and reset values)
    reset!(cl,box)
  
    #
    # Add particles to cell list
    #
    nt = min(length(x)÷250,nthreads())
    if !parallel || nt < 2
        cl = add_particles!(x,box,0,cl)
    else
        # Indices of the atoms that will be added by each thread
        idx_thread = Vector{UnitRange}(undef,nt)
        nperthread = length(x)÷nt
        nrem = length(x) - nt*nperthread
        first = 1
        for it in 1:nt
            nx = nperthread
            (it <= nrem) && (nx += 1)
            idx_thread[it] = first:(first-1)+nx
            first += nx
        end
        # Cell lists to be built by each thread
        clt = [ init_cell_list(@view(x[idx_thread[it]]),box) for it in 1:nt ]
        @threads for it in 1:nt
            reset!(clt[it],box)
            xt = @view(x[idx_thread[it]])  
            clt[it] = add_particles!(xt,box,idx_thread[it][1]-1,clt[it])
        end
        #
        # Merge threaded cell lists
        #
        @inbounds for it in 1:nt
            voltar
            cl.ncp[1] += clt[it].ncp[1]
            for (icwp,cell) in pairs(clt[it].cwp)
                if cell.number_of_particles > 0
                    if cell.contains_real
                        if !cl.cwp[icell].contains_real
                            cl.number_of_cells_with_particles += 1
                            @set! cl.cwp[icell].contains_real = true
                            @set! cl.cwp[icell].index = cell.index
                        end
                    end
                    nprev = cl.cwp[cl.number_of_cells_with_particles]
                    cl.cwp[cl.cnwp[1]] = npcell[icell] += clt[it].npcell[icell]
                    if cl.npcell[icell] > length(cl.list[icell]) 
                        resize!(cl.list[icell],cl.npcell[icell])
                    end
                    for ip in 1:clt[it].npcell[icell]
                        cl.list[icell][nprev + ip] = clt[it].list[icell][ip] 
                    end
                end
            end
        end
    end
  
    maximum_npcell = maximum(npcell)
    if maximum_npcell > length(projected_particles[1])
        for i in 1:nthreads()
            resize!(projected_particles[i],ceil(Int,1.2*maximum_npcell))
        end
    end
  
    return cl
end

function add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}
    #
    # Add virtual particles to edge cells
    #
    for ip in eachindex(x)
        xp = x[ip] 
        p = SVector{N,T}(ntuple(i->xp[i],N)) # in case the input was not static
        p = wrap_to_first(p,box)
        cl = replicate_particle!(ishift+ip,p,box,cl)
    end
    #
    # Add true particles, such that the first particle of each cell is
    # always a true particle
    #
    for ip in eachindex(x)
        xp = x[ip]
        p = SVector{N,T}(ntuple(i->xp[i],N))
        p = wrap_to_first(p,box)
        cl = add_particle_to_celllist!(ishift+ip,p,box,cl) 
    end
    return cl
end

"""

```
add_particle_to_celllist!(
    ip,
    x::SVector{N,T},
    box,
    cl::CellList{N,T};
    real_particle::Bool=true
) where {N,T}
```

Adds one particle to the cell lists, updating all necessary arrays.

"""
function add_particle_to_celllist!(
    ip,
    x::SVector{N,T},
    box,
    cl::CellList{N,T};
    real_particle::Bool=true
) where {N,T}
  @unpack cell_index_in_list,  
          number_of_particles,
          number_of_cells_with_real_particles
          list = cl

    number_of_particles += 1
    icell_cartesian = particle_cell(x,box)
    icell = cell_linear_index(box.nc,icell_cartesian)
    cell = list[icell]
    #
    # Cells with real particles are annotated to be run over
    #
    if real_particle && (!cell.contains_real)
        cell.contains_real = true
        number_of_cells_with_real_particles += 1
        cell_index_in_list[icell] = number_of_cells_with_real_particles
        cell.icell = icell
        cell.cartesian_index = icell_cartesian
        cell.center = cell_center(icell_cartesian,box)
    end

    cell.number_of_particles += 1
    if cell.number_of_particles > length(cell.particles)
        resize!(cell.number_of_particles,ceil(Int,2*length(cell.particles)))
    end
    cell.particles[cell.number_of_particles] = ParticleWithIndex(ip,x,real_particle) 

    return cl
end

"""

```
UpdateCellList!(
  x::AbstractVector{<:AbstractVector},
  y::AbstractVector{<:AbstractVector},
  box::Box{UnitCellType,N,T},cl:CellListPair,parallel=true
) where {UnitCellType,N,T}
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
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T},cl_pair::CellListPair;
    parallel::Bool=true
) where {UnitCellType,N,T}

    if length(x) <= length(y)
        UpdateCellList!(y,box,cl_pair.large,parallel=parallel)
    else
        UpdateCellList!(x,box,cl_pair.large,parallel=parallel)
    end

    return cl_pair
end



