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

Structure that contains the maximum lengths on each direction,
to dispatch on the construction of boxes without periodic boundary
conditions.

"""
struct Limits{T<:AbstractVector} 
    limits::T
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains some data required to compute the linked cells. To
be initialized with the box size and cutoff. 

## Examples

```julia-repl
julia> using CellListMap

julia> sides = [250,250,250];

julia> cutoff = 10;

julia> box = Box(sides,cutoff)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [250.0 0.0 0.0; 0.0 250.0 0.0; 0.0 0.0 250.0]
  cutoff: 10.0
  number of computing cells on each dimension: [27, 27, 27]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 19683

```

```julia-repl
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
```julia-repl
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
    println(io,typeof(box))
    println(io,"  unit cell matrix: ", box.unit_cell.matrix) 
    println(io,"  cutoff: ", box.cutoff)
    println(io,"  number of computing cells on each dimension: ",box.nc)
    println(io,"  computing cell sizes: ", box.cell_size, " (lcell: ",box.lcell,")")
    print(io,"  Total number of cells: ", prod(box.nc))
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
```julia-repl
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

```
Box(
    limits::Limits,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1
)
```

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus the cutoff. This is used to 
emulate pairwise interactions in non-periodic boxes.

### Examples

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x),10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [109.99633932875878 0.0 0.0; 0.0 109.99780283179763 0.0; 0.0 0.0 109.99587254766517]
  cutoff: 10.0
  number of computing cells on each dimension: [12, 12, 12]
  computing cell sizes: [10.999633932875877, 10.999780283179764, 10.999587254766517] (lcell: 1)
  Total number of cells: 1728

julia> y = [ [150,150,50] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x,y),10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [159.99787690924168 0.0 0.0; 0.0 159.98878289444897 0.0; 0.0 0.0 109.99587254766517]
  cutoff: 10.0
  number of computing cells on each dimension: [18, 17, 12]
  computing cell sizes: [10.666525127282778, 10.665918859629931, 10.999587254766517] (lcell: 1)
  Total number of cells: 3672

```

"""
Box(
    limits::Limits,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1
) = Box(limits.limits .+ cutoff,cutoff,T,lcell) 

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
    center::SVector{N,T} = zeros(SVector{N,T})
    contains_real::Bool = false
    n_particles::Int = 0
    particles::Vector{ParticleWithIndex{N,T}} = Vector{ParticleWithIndex{N,T}}(undef,0)
end
function Cell{N,T}(cartesian_index::CartesianIndex,box::Box) where {N,T}
    return Cell{N,T}(
        linear_index=cell_linear_index(box.nc,cartesian_index),
        cartesian_index=cartesian_index,
        center=cell_center(cartesian_index,box)
    )
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Auxiliary structure to contain projected particles.

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
    cell_indices_real::Vector{Int} = zeros(Int,0)
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
    ref::V
    target::CellList{N,T}
    swap::Bool
end      
function Base.show(io::IO,::MIME"text/plain",cl::CellListPair)
    print(io,typeof(cl),"\n")
    print(io,"   $(length(cl.ref)) particles in the reference vector.\n")
    print(io,"   $(cl.target.n_cells_with_real_particles) cells with real particles of target vector.")
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
        resize!(cl.cell_indices,number_of_cells)
        @. cl.cell_indices = 0
        @. cl.cell_indices_real = 0
    end
    for i in 1:cl.n_cells_with_particles
        cl.cells[i] = Cell{N,T}(
            n_particles=0,
            contains_real=false,
            particles=cl.cells[i].particles
        )
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
    parallel::Bool=true,
    autoswap::Bool=true
) where {UnitCellType,N,T} 
```

Function that will initialize a `CellListPair` structure from scracth, given two vectors
of particle coordinates and a `Box`, which contain the size of the system, cutoff, etc.
By default, the cell list will be constructed for smallest vector, but this is not always
the optimal choice. Using `autoswap=false` the cell list is constructed for the second (`y`)

### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> y = [ 250*rand(SVector{3,Float64}) for i in 1:10000 ];

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   10000 particles in the reference vector.
   961 cells with real particles of target vector.

julia> cl = CellList(x,y,box,autoswap=false)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   1000 particles in the reference vector.
   7389 cells with real particles of target vector.

```

"""
function CellList(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true,
    autoswap=true
) where {UnitCellType,N,T} 
    if length(x) >= length(y) || !autoswap
        ref = [ SVector{N,T}(ntuple(i->el[i],N)...) for el in x ]
        target = CellList(y,box,parallel=parallel)
        swap = false
    else
        ref = [ SVector{N,T}(ntuple(i->el[i],N)...) for el in y ]
        target = CellList(x,box,parallel=parallel)
        swap = true
    end
    cl_pair = CellListPair(ref=ref,target=target,swap=swap)
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
    box::Box,
    cl:CellList{N,T},
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
    parallel::Bool=true,
    threaded_lists=nothing,
    threaded_idxs=nothing
) where {N,T}

    # Reset base cell list cell (resize if needed, and reset values)
    cl = reset!(cl,box)

    #
    # Add particles to cell list
    #
    nt = min(length(x)÷250,nthreads())
    if !parallel || nt < 2
        cl = add_particles!(x,box,0,cl)
    else
        if isnothing(threaded_lists)
            threaded_lists = [ deepcopy(cl) for _ in 1:nt-1] 
            for it in 1:nt-1
                reset!(threaded_lists[it],box)
            end
        end
        if isnothing(threaded_idxs)
            threaded_idxs = Vector{UnitRange}(undef,nt)
        end
        # Indices of the atoms that will be added by each thread
        nperthread = length(x)÷nt
        nrem = length(x) - nt*nperthread
        first = 1
        for it in 1:nt
            nx = nperthread
            (it <= nrem) && (nx += 1)
            threaded_idxs[it] = first:(first-1)+nx
            first += nx
        end
        # Cell lists to be built by each thread
        @threads for it in 1:nt
            if it < nt
                xt = @view(x[threaded_idxs[it+1]])  
                threaded_lists[it] = add_particles!(xt,box,threaded_idxs[it][1]-1,threaded_lists[it])
            else
                xt = @view(x[threaded_idxs[1]])  
                cl = add_particles!(xt,box,0,cl)
            end
        end
        #
        # Merge threaded cell lists
        #
        for it in 1:nt-1
            # Accumulate number of particles
            @set! cl.n_particles += threaded_lists[it].n_particles
            for icell in 1:threaded_lists[it].n_cells_with_particles
                cell = threaded_lists[it].cells[icell]
                linear_index = cell.linear_index
                # If cell was yet not initialized in merge, push it to the list
                if cl.cell_indices[linear_index] == 0
                    @set! cl.n_cells_with_particles += 1
                    if length(cl.cells) >= cl.n_cells_with_particles
                        cl.cells[cl.n_cells_with_particles] = cell 
                    else
                        push!(cl.cells,cell)
                    end
                    cl.cell_indices[linear_index] = cl.n_cells_with_particles
                    if cell.contains_real
                        @set! cl.n_cells_with_real_particles += 1
                        if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                            push!(cl.cell_indices_real,cl.cell_indices[linear_index])
                        else
                            cl.cell_indices_real[cl.n_cells_with_real_particles] = cl.cell_indices[linear_index] 
                        end
                    end
                # Append particles to initialized cells
                else
                    cell_index = cl.cell_indices[linear_index]
                    prevcell = cl.cells[cell_index] 
                    n_particles_old = prevcell.n_particles
                    @set! prevcell.n_particles += cell.n_particles
                    if prevcell.n_particles > length(prevcell.particles)
                        resize!(prevcell.particles,prevcell.n_particles)
                    end
                    for ip in 1:cell.n_particles
                        prevcell.particles[n_particles_old+ip] = cell.particles[ip]
                    end
                    cl.cells[cell_index] = prevcell
                    if (!cl.cells[cl.cell_indices[linear_index]].contains_real) && cell.contains_real 
                        cl_cell = cl.cells[cl.cell_indices[linear_index]]
                        @set! cl_cell.contains_real = true
                        cl.cells[cl.cell_indices[linear_index]] = cl_cell
                        @set! cl.n_cells_with_real_particles += 1
                        if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                            push!(cl.cell_indices_real,cl.cell_indices[linear_index])
                        else
                            cl.cell_indices_real[cl.n_cells_with_real_particles] = cl.cell_indices[linear_index]
                        end
                    end
                end
            end
        end
    end
  
    maxnp = 0
    for i in 1:cl.n_cells_with_particles
        maxnp = max(maxnp,cl.cells[i].n_particles)
    end
    if maxnp > length(cl.projected_particles[1])
        for i in 1:nthreads()
            resize!(cl.projected_particles[i],ceil(Int,1.2*maxnp))
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
    cartesian_index = particle_cell(x,box)
    linear_index = cell_linear_index(box.nc,cartesian_index)

    #
    # Check if this is the first particle of this cell
    #
    if cell_indices[linear_index] == 0
        n_cells_with_particles += 1
        cell_indices[linear_index] = n_cells_with_particles
        if n_cells_with_particles > length(cells)
            cell = Cell{N,T}(cartesian_index,box)
        else
            cell = cells[n_cells_with_particles]
            @set! cell.linear_index = linear_index
            @set! cell.cartesian_index = cartesian_index
            @set! cell.center = cell_center(cell.cartesian_index,box)
            @set! cell.contains_real = false
        end
        @set! cell.n_particles = 1
    else
        cell = cells[cell_indices[linear_index]]
        @set! cell.n_particles += 1
    end
    #
    # Cells with real particles are annotated to be run over
    #
    if real_particle && (!cell.contains_real)
        @set! cell.contains_real = true
        n_cells_with_real_particles += 1
        if n_cells_with_real_particles > length(cell_indices_real)
            push!(cell_indices_real,cell_indices[linear_index])
        else
            cell_indices_real[n_cells_with_real_particles] = cell_indices[linear_index] 
        end
    end

    #
    # Add particle to cell list
    #
    p = ParticleWithIndex(ip,x,real_particle) 
    if cell.n_particles > length(cell.particles)
        push!(cell.particles,p)
    else
        cell.particles[cell.n_particles] = p
    end

    #
    # Update (imutable) cell in list
    #
    @set! cl.n_particles = n_particles
    @set! cl.cell_indices = cell_indices
    @set! cl.cell_indices_real = cell_indices_real
    @set! cl.n_cells_with_particles = n_cells_with_particles
    @set! cl.n_cells_with_real_particles = n_cells_with_real_particles
    if n_cells_with_particles > length(cl.cells)
        push!(cl.cells,cell)
    else
        cl.cells[n_cells_with_particles] = cell
    end

    return cl
end

"""

```
UpdateCellList!(
  x::AbstractVector{<:AbstractVector},
  y::AbstractVector{<:AbstractVector},
  box::Box{UnitCellType,N,T},
  cl:CellListPair,
  parallel=true
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
        cl_pair = UpdateCellList!(y,box,cl_pair.target,parallel=parallel)
    else
        cl_pair = UpdateCellList!(x,box,cl_pair.target,parallel=parallel)
    end

    return cl_pair
end



