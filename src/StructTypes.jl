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
    index::Int
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
Base.@kwdef struct Cell{N,T}
    linear_index::Int = 0
    cartesian_index::CartesianIndex{N} = CartesianIndex{N}(ntuple(i->0,N))
    center::SVector{N,T} = zeros{SVector{N,T}}
    contains_real::Bool = false
    n_particles::Int = 0
    particles::Vector{ParticleWithIndex{N,T}} = Vector{ParticleWithIndex{N,T}}(undef,0)
end
function Cell{N,T}(icell_cartesian::CartesianIndex,box::Box) where {N,T}
    return Cell{N,T}(
        linear_index=cell_linear_index(box.nc,icell_cartesian),
        cartesian_index=icell_cartesian,
        center=cell_center(icell_cartesian,box)
    )
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Auxiliary structure to contain projected particles in large and
and dense systems.

"""
Base.@kwdef struct ProjectedParticle{N,T}
    index::Int = 0
    xproj::T = zero(T)
    coordinates::SVector{N,T} = zeros(SVector{N,T})
    real::Bool = false
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef struct CellList{N,T}
    " Number of cells. "
    number_of_cells::Int = 0
    " *mutable* number of particles in the computing box. "
    n_particles::Int = 0
    " *mutable* number of cells with real particles. "
    n_cells_with_real_particles::Int = 0
    " *mutable* number of cells with particles, real or images. "
    n_cells_with_particles::Int = 0
    " Auxiliary array that contains the indexes in list of the cells with particles, real or images. "
    cell_indices::Vector{Int} = zeros(Int,number_of_cells)
    " Auxiliary array that contains the indexes in the cells with real particles. "
    cell_indices_real::Vector{Int} = zeros(Int,number_of_cells)
    " Vector containing cell lists of cells with particles. "
    cells::Vector{Cell{N,T}} = Cell{N,T}[]
    " Auxiliar array to store projected particles. "
    projected_particles::Vector{Vector{ProjectedParticle{N,T}}} = 
        [ Vector{ProjectedParticle{N,T}}(undef,0) for _ in 1:nthreads() ]
end
function Base.show(io::IO,::MIME"text/plain",cl::CellList)
    println(typeof(cl))
    println("  $(cl.n_cells_with_real_particles) cells with real particles.")
    println("  $(cl.n_particles) particles in computing box, including images.")
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
    print("   $(cl.large.n_cells_with_real_particles) cells with real particles.")
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
    cl = CellList{N,T}(number_of_cells=prod(box.nc))
    return UpdateCellList!(x,box,cl,parallel=parallel)
end

function reset!(cl::CellList{N,T},box) where{N,T}
    number_of_cells = ceil(Int,1.2*prod(box.nc)) # some margin in case of box size variations
    if number_of_cells > length(cl.cells) 
        resize!(cl.cells,number_of_cells)
        resize!(cl.cell_indices,number_of_cells)
        @. cl.cell_indices = 0
        resize!(cl.cell_indices_real,number_of_cells)
        @. cl.cell_indices_real = 0
    end
    for i in 1:cl.n_cells_with_particles
        index = cl.cell_indices[i]
        cl.cells[index] = Cell{N,T}(n_particles=0,particles=cl.cells[i].particles)
    end
    cl = CellList{N,T}(
        n_particles = 0,
        number_of_cells = number_of_cells,
        n_cells_with_real_particles = 0,
        n_cells_with_particles = 0,
        cell_indices = cl.cell_indices,
        cell_indices_real = cl.cell_indices_real,
        cells=cl.cells,
        projected_particles = cl.projected_particles
    )
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
    @unpack cells, projected_particles = cl

    # Reset cell (resize if needed, and reset values)
    cl = reset!(cl,box)

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
        clt = [ CellList{N,T}(number_of_cells=prod(box.nc)) for _ in 1:nt ]
        @threads for it in 1:nt
            xt = @view(x[idx_thread[it]])  
            clt[it] = add_particles!(xt,box,idx_thread[it][1]-1,clt[it])
        end
        #
        # Merge threaded cell lists
        #
        for it in 1:nt
            # Accumulate number of particles
            @set! cl.n_particles += clt[it].n_particles
            for icell in 1:clt[it].number_of_cells
                list_index = clt[it].cell_indices[icell]
                if list_index == 0 # empty cell
                    continue
                end
                cell = clt[it].cells[list_index]
                linear_index = cell.linear_index
                # If cell was yet not initialized in merge, copy 
                if cl.cell_indices[linear_index] == 0
                    @set! cl.n_cells_with_particles += 1
                    if length(cl.cells) >= cl.n_cells_with_particles
                        cl.cells[cl.n_cells_with_particles] = cell 
                    else
                        push!(cl.cells,cell)
                    end
                    cl.cell_indices[linear_index] = cl.n_cells_with_particles
                # Append particles to initialized cells
                else
                    cell_index = cl.cell_indices[linear_index]
                    prevcell = cl.cells[cell_index] 
                    @set! prevcell.n_particles += cell.n_particles
                    append!(prevcell.particles,cell.particles)
                    cl.cells[cell_index] = prevcell
                end
                # Check if this a new cell with real particles
                if (!cl.cells[cl.cell_indices[linear_index]].contains_real) && cell.contains_real 
                    @set! cl.n_cells_with_real_particles += 1
                    cl.cell_indices_real[cl.n_cells_with_real_particles] = cl.cell_indices[linear_index]
                end
            end
        end
    end
  
    maxnp = 0
    for i in 1:cl.n_cells_with_particles
        maxnp = max(maxnp,cl.cells[i].n_particles)
    end
    if maxnp > length(projected_particles[1])
        for i in 1:nthreads()
            resize!(projected_particles[i],ceil(Int,1.2*maxnp))
        end
    end
  
    return cl
end

function add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}
    #
    # Add particles to cell lists
    #
    for ip in eachindex(x)
        xp = x[ip] 
        p = SVector{N,T}(ntuple(i->xp[i],N)) # in case the input was not static
        p = wrap_to_first(p,box)
        cl = add_particle_to_celllist!(ishift+ip,p,box,cl) # add real particle
        cl = replicate_particle!(ishift+ip,p,box,cl) # add virtual particles to border cells
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
    @unpack n_particles,
            n_cells_with_real_particles,
            n_cells_with_particles,
            cell_indices,
            cell_indices_real,
            cells = cl

    # Cell of this particle
    n_particles += 1 
    icell_cartesian = particle_cell(x,box)
    icell_linear = cell_linear_index(box.nc,icell_cartesian)

    #
    # If cell was empty, initialize cell
    #
    if cell_indices[icell_linear] == 0
        n_cells_with_particles += 1
        cell_indices[icell_linear] = n_cells_with_particles
        cell = Cell{N,T}(icell_cartesian,box)
    else
        cell = cells[cell_indices[icell_linear]]
    end
    #
    # Cells with real particles are annotated to be run over
    #
    if real_particle && (!cell.contains_real)
        @set! cell.contains_real = true
        n_cells_with_real_particles += 1
        cell_indices_real[n_cells_with_real_particles] = cell_indices[icell_linear] 
    end

    #
    # Add particle to cell list
    #
    @set! cell.n_particles += 1
    if cell.n_particles > length(cell.particles)
        resize!(cell.particles,ceil(Int,2*cell.n_particles))
    end
    cell.particles[cell.n_particles] = ParticleWithIndex(ip,x,real_particle) 

    #
    # Update (imutable) cell in list
    #
    @set! cl.n_particles = n_particles
    @set! cl.cell_indices = cell_indices
    @set! cl.cell_indices_real = cell_indices_real
    @set! cl.n_cells_with_particles = n_cells_with_particles
    @set! cl.n_cells_with_real_particles = n_cells_with_real_particles
    if length(cl.cells) >= cell_indices[icell_linear]
        cl.cells[cell_indices[icell_linear]] = cell
    else
        push!(cl.cells,cell)
    end

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
        cl_pair = UpdateCellList!(y,box,cl_pair.large,parallel=parallel)
    else
        cl_pair = UpdateCellList!(x,box,cl_pair.large,parallel=parallel)
    end

    return cl_pair
end



