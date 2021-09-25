#
# This file contains all structre types and functions necessary for building
# the CellList and CellListPair structures.
#

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort. The index is a 32bit integer such 
that the complete struct has 32bytes.

"""
struct ParticleWithIndex{N,T}
    index::Int32
    real::Bool
    coordinates::SVector{N,T}
end
Base.zero(::Type{ParticleWithIndex{N,T}}) where {N,T} =
    ParticleWithIndex{N,T}(Int32(0),false,zeros(SVector{N,T}))

"""

```
set_nt(cl) = max(1,min(cl.n_real_particles÷500,nthreads()))
```

Don't use all threads to build cell lists if the number of particles
per thread is smaller than 500.

"""
set_nt(cl) = max(1,min(cl.n_real_particles÷500,nthreads()))

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

Auxiliary structure to contain projected particles. Types of 
scalars are chosen such that with a `SVector{3,Float64}` the
complete struct has 32bytes.

"""
Base.@kwdef struct ProjectedParticle{N,T}
    index::Int32
    xproj::Float32
    coordinates::SVector{N,T}
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that contains the cell lists information.

"""
Base.@kwdef struct CellList{N,T}
    " Number of real particles. "
    n_real_particles::Int = 0
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
    println(io,typeof(cl))
    println(io,"  $(cl.n_real_particles) real particles.")
    println(io,"  $(cl.n_cells_with_real_particles) cells with real particles.")
    print(io,"  $(cl.n_particles) particles in computing box, including images.")
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

$(TYPEDEF)

$(TYPEDEF)

Auxiliary structure to carry threaded lists and ranges of particles to 
be considered by each thread on parallel construction. 

"""
@with_kw struct AuxThreaded{N,T}
    idxs::Vector{UnitRange{Int}} = Vector{UnitRange{Int}}(undef,0)
    lists::Vector{CellList{N,T}} = Vector{CellList{N,T}}(undef,0)
end
function Base.show(io::IO,::MIME"text/plain",aux::AuxThreaded)
    println(io,typeof(aux))
    print(io," Auxiliary arrays for nthreads = ", length(aux.lists)) 
end

"""

```
AuxThreaded(cl::CellList{N,T}) where {N,T}
```

Constructor for the `AuxThreaded` type, to be passed to `UpdateCellList!` for in-place 
update of cell lists. 

## Example
```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(3) for _ in 1:100_000 ];

julia> cl = CellList(x,box);

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> cl = UpdateCellList!(x,box,cl,aux)
CellList{3, Float64}
  100000 real particles.
  31190 cells with real particles.
  1134378 particles in computing box, including images.

```
"""
function AuxThreaded(cl::CellList{N,T}) where {N,T}
    aux = AuxThreaded{N,T}()
    init_aux_threaded!(aux,cl)
    return aux
end

"""

```
AuxThreaded(cl::CellListPair{N,T}) where {N,T}
```

Constructor for the `AuxThreaded` type for lists of disjoint particle sets, 
to be passed to `UpdateCellList!` for in-place update of cell lists. 

## Example
```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(3) for i in 1:50_000 ];

julia> y = [ 250*rand(3) for i in 1:10_000 ];

julia> cl = CellList(x,y,box);

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> cl = UpdateCellList!(x,box,cl,aux)
CellList{3, Float64}
  100000 real particles.
  31190 cells with real particles.
  1134378 particles in computing box, including images.

```
"""
function AuxThreaded(cl_pair::CellListPair{V,N,T}) where {V,N,T}
    aux = AuxThreaded{N,T}()
    init_aux_threaded!(aux,cl_pair.target)
    return aux
end

"""

```
init_aux_threaded!(aux::AuxThreaded,cl::CellList)
```

Given an `AuxThreaded` object initialized with zero-length arrays,
push `ntrheads` copies of `cl` into `aux.lists` and resize `aux.idxs`
to the number of threads.  

"""
function init_aux_threaded!(aux::AuxThreaded,cl::CellList)
   nt = set_nt(cl)
   push!(aux.lists,cl)
   for it in 2:nt
       push!(aux.lists,deepcopy(cl))
   end
   # Indices of the atoms that will be added by each thread
   nrem = cl.n_real_particles%nt
   nperthread = (cl.n_real_particles-nrem)÷nt
   first = 1
   for it in 1:nt
       nx = nperthread
       if it <= nrem
           nx += 1
       end
       push!(aux.idxs,first:(first-1)+nx)
       first += nx
       cl_tmp = aux.lists[it]
       @set! cl_tmp.n_real_particles = length(aux.idxs[it]) 
       aux.lists[it] = cl_tmp
   end
   return aux
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

julia> cl = CellList(x,box)
CellList{3, Float64}
  100000 real particles.
  15600 cells with real particles.
  126276 particles in computing box, including images.

```

"""
function CellList(
    x::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
    cl = CellList{N,T}(
        n_real_particles=length(x),
        number_of_cells=prod(box.nc)
    )
    return UpdateCellList!(x,box,cl,parallel=parallel)
end

"""

```
function CellList(
    x::AbstractMatrix,
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
```

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

"""
function CellList(
    x::AbstractMatrix,
    box::Box{UnitCellType,N,T};
    parallel::Bool=true
) where {UnitCellType,N,T} 
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    cl = CellList{N,T}(
        n_real_particles=length(x_re),
        number_of_cells=prod(box.nc)
    )
    return UpdateCellList!(x_re,box,cl,parallel=parallel)
end


"""

```
reset!(cl::CellList{N,T},box) where{N,T}
```

Resets a cell list, by setting everything to zero, but retaining
the allocated `particles` and `projected_particles` vectors.

"""
function reset!(cl::CellList{N,T},box,n_real_particles) where{N,T}
    new_number_of_cells = prod(box.nc) 
    if new_number_of_cells > cl.number_of_cells
        resize!(cl.cell_indices,new_number_of_cells)
    end
    for i in 1:length(cl.cells)
        cl.cells[i] = Cell{N,T}(particles=cl.cells[i].particles)
    end
    @. cl.cell_indices = 0
    @. cl.cell_indices_real = 0
    cl = CellList{N,T}(
        n_real_particles = n_real_particles, 
        n_particles = 0,
        number_of_cells = new_number_of_cells,
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
function CellList(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T};
    parallel::Bool=true,
    autoswap=true
) where {UnitCellType,N,T} 
```

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function CellList(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T};
    parallel::Bool=true,
    autoswap=true
) where {UnitCellType,N,T} 
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    @assert size(y,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    y_re = reinterpret(reshape, SVector{N,T}, x)
    CellList(x_re,y_re,box,parallel=parallel,autoswap=autoswap)
end

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
updated particle positions. This function will allocate new threaded auxiliary
arrays in parallel calculations. To preallocate these auxiliary arrays, use
the `UpdateCellList!(x,box,cl,aux)` method instead. 

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
    aux = AuxThreaded{N,T}()
    if parallel && nthreads() > 1
        init_aux_threaded!(aux,cl)
    end
    return UpdateCellList!(x,box,cl,aux,parallel=parallel)
end

"""

```
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T};
    parallel::Bool=true
) where {N,T}
```

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T};
    parallel::Bool=true
) where {N,T}
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    return UpdateCellList!(x_re,box,cl,parallel=parallel)
end

"""

```
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    box::Box,
    cl::CellList{N,T},
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {N,T}
```

Function that updates the cell list `cl` new coordinates `x` and possibly a new
box `box`, and receives a preallocated `aux` structure of auxiliary vectors for
threaded cell list construction. Given a preallocated `aux` vector, allocations in
this function should be minimal, only associated with the spawning threads, or
to expansion of the cell lists if the number of cells or number of particles 
increased. 

### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:100000 ];

julia> cl = CellList(x,box);

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:100000 ];

julia> UpdateCellList!(x,box,cl,aux)
CellList{3, Float64}
  100000 real particles.
  15599 cells with real particles.
  125699 particles in computing box, including images.

```

To illustrate the expected ammount of allocations, which are a consequence
of thread spawning only:

```julia-repl
julia> using BenchmarkTools

julia> @btime UpdateCellList!(\$x,\$box,\$cl,\$aux)
  16.384 ms (41 allocations: 3.88 KiB)
CellList{3, Float64}
  100000 real particles.
  15599 cells with real particles.
  125699 particles in computing box, including images.

julia> @btime UpdateCellList!(\$x,\$box,\$cl,\$aux,parallel=false)
  20.882 ms (0 allocations: 0 bytes)
CellList{3, Float64}
  100000 real particles.
  15603 cells with real particles.
  125896 particles in computing box, including images.

```

"""
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    box::Box,
    cl::CellList{N,T},
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {N,T}


    # Add particles to cell list
    nt = set_nt(cl)
    if !parallel || nt < 2
        cl = reset!(cl,box,length(x))
        cl = add_particles!(x,box,0,cl)
    else
        # Cell lists to be built by each thread
        @threads for it in 1:nt
             aux.lists[it] = reset!(aux.lists[it],box,length(aux.idxs[it]))
             xt = @view(x[aux.idxs[it]])  
             aux.lists[it] = add_particles!(xt,box,aux.idxs[it][1]-1,aux.lists[it])
        end
        # Merge threaded cell lists
        cl = aux.lists[1]
        for it in 2:nt
            cl = merge_cell_lists!(cl,aux.lists[it])
        end
    end
  
    # allocate, or update the auxiliary projected_particles arrays
    maxnp = 0
    for i in 1:cl.n_cells_with_particles
        maxnp = max(maxnp,cl.cells[i].n_particles)
    end
    if maxnp > length(cl.projected_particles[1])
        for i in 1:nthreads()
            resize!(cl.projected_particles[i],maxnp)
        end
    end

    return cl
end

"""

```
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T},
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {N,T}
```

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T},
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {N,T}
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    return UpdateCellList!(x_re,box,cl,aux,parallel=parallel)
end

"""

```
add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}
```

Add all particles in vector `x` to the cell list `cl`. `ishift` is the shift in particle
index, meaning that particle `i` of vector `x` corresponds to the particle with original
index `i+ishift`. The shift is used to construct cell lists from fractions of the original
set of particles in parallel list construction.  

"""
function add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}
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
merge_cell_lists!(cl::CellList,aux::CellList)
```

Merges an auxiliary `aux` cell list to `cl`, and returns the modified `cl`. Used to
merge cell lists computed in parallel threads.

"""
function merge_cell_lists!(cl::CellList,aux::CellList)
    # One should never get here if the lists do not share the same
    # computing box
    @assert cl.number_of_cells == aux.number_of_cells
    # Accumulate number of particles
    @set! cl.n_particles += aux.n_particles
    @set! cl.n_real_particles += aux.n_real_particles
    for icell in 1:aux.n_cells_with_particles
        cell = aux.cells[icell]
        linear_index = cell.linear_index
        cell_index = cl.cell_indices[linear_index]
        # If cell was yet not initialized in merge, push it to the list
        if cell_index == 0
            @set! cl.n_cells_with_particles += 1
            cell_index = cl.n_cells_with_particles
            cl.cell_indices[linear_index] = cell_index
            if cell_index > length(cl.cells)
                push!(cl.cells,deepcopy(cell))
            else
                prevcell = cl.cells[cell_index]
                @set! prevcell.linear_index = cell.linear_index
                @set! prevcell.cartesian_index = cell.cartesian_index
                @set! prevcell.center = cell.center
                @set! prevcell.n_particles = cell.n_particles
                @set! prevcell.contains_real = cell.contains_real
                for ip in 1:cell.n_particles
                    p = cell.particles[ip]
                    if ip > length(prevcell.particles)
                        push!(prevcell.particles,p)
                    else
                        prevcell.particles[ip] = p
                    end
                end
                cl.cells[cell_index] = prevcell 
            end
            if cell.contains_real
                @set! cl.n_cells_with_real_particles += 1
                if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                    push!(cl.cell_indices_real,cell_index)
                else
                    cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
                end
            end
        # Append particles to initialized cells
        else
            prevcell = cl.cells[cell_index] 
            n_particles_old = prevcell.n_particles
            @set! prevcell.n_particles += cell.n_particles
            if prevcell.n_particles > length(prevcell.particles)
                resize!(prevcell.particles,prevcell.n_particles)
            end
            for ip in 1:cell.n_particles
                prevcell.particles[n_particles_old+ip] = cell.particles[ip]
            end
            # If the previous cell didn't contain real particles, but the current one
            # does, update
            if !prevcell.contains_real && cell.contains_real 
                @set! prevcell.contains_real = true
                @set! cl.n_cells_with_real_particles += 1
                if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                    push!(cl.cell_indices_real,cell_index)
                else
                    cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
                end
            end
            cl.cells[cell_index] = prevcell
        end
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
    @unpack n_cells_with_real_particles,
            n_cells_with_particles,
            cell_indices,
            cell_indices_real,
            cells = cl

    # Increase the counter of number of particles of this list
    @set! cl.n_particles += 1
    # Cell of this particle
    cartesian_index = particle_cell(x,box)
    linear_index = cell_linear_index(box.nc,cartesian_index)

    # Check if this is the first particle of this cell, if it is,
    # initialize a new cell, or reset a previously allocated one
    cell_index = cell_indices[linear_index]

    if cell_index == 0
        n_cells_with_particles += 1
        cell_index = n_cells_with_particles
        cell_indices[linear_index] = cell_index
        if cell_index > length(cells)
            push!(cells,Cell{N,T}(cartesian_index,box))
        else
            cell = cells[cell_index]
            @set! cell.linear_index = linear_index
            @set! cell.cartesian_index = cartesian_index
            @set! cell.center = cell_center(cell.cartesian_index,box)
            @set! cell.contains_real = false
            @set! cell.n_particles = 0
            cells[cell_index] = cell
        end
    end

    # Increase particle counter for this cell
    cell = cells[cell_index]
    @set! cell.n_particles += 1

    #
    # Cells with real particles are annotated to be run over
    #
    if real_particle && (!cell.contains_real)
        @set! cell.contains_real = true
        n_cells_with_real_particles += 1
        if n_cells_with_real_particles > length(cell_indices_real)
            push!(cell_indices_real,cell_index)
        else
            cell_indices_real[n_cells_with_real_particles] = cell_index 
        end
    end

    #
    # Add particle to cell list
    #
    p = ParticleWithIndex(Int32(ip),real_particle,x) 
    if cell.n_particles > length(cell.particles)
        push!(cell.particles,p)
    else
        cell.particles[cell.n_particles] = p
    end

    #
    # Update (imutable) cell in list
    #
    @set! cl.n_cells_with_particles = n_cells_with_particles
    @set! cl.n_cells_with_real_particles = n_cells_with_real_particles
    cells[cell_index] = cell

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

Function that will update a previously allocated `CellListPair` structure, given 
new updated particle positions, for example. This method will allocate new 
`aux` threaded auxiliary arrays. For a non-allocating version, see the 
`UpdateCellList!(x,y,box,cl,aux)` method.

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
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair;
    parallel::Bool=true
) where {UnitCellType,N,T}
    aux = AuxThreaded{N,T}()
    if parallel && nthreads() > 1
        init_aux_threaded!(aux,cl_pair.target)
    end
    return UpdateCellList!(x,y,box,cl_pair,aux,parallel=parallel)
end

"""

```
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair;
    parallel::Bool=true
) where {UnitCellType,N,T}
```

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair;
    parallel::Bool=true
) where {UnitCellType,N,T}
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    @assert size(y,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    y_re = reinterpret(reshape, SVector{N,T}, x)
    return UpdateCellList!(x_re,y_re,box,cl_pair,parallel=parallel)
end

"""

```
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair,
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {UnitCellType,N,T}
```

This function will update the `cl_pair` structure that contains the cell lists
for disjoint sets of particles. It receives the preallocated `aux` structure to
avoid reallocating auxiliary arrays necessary for the threaded construct of the
lists. 

### Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(3) for i in 1:50_000 ];

julia> y = [ 250*rand(3) for i in 1:10_000 ];

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   50000 particles in the reference vector.
   7381 cells with real particles of target vector.

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> x = [ 250*rand(3) for i in 1:50_000 ];

julia> y = [ 250*rand(3) for i in 1:10_000 ];

julia> cl = UpdateCellList!(x,y,box,cl,aux)
CellList{3, Float64}
  10000 real particles.
  7358 cells with real particles.
  12591 particles in computing box, including images.

```
To illustrate the expected ammount of allocations, which are a consequence
of thread spawning only:

```julia-repl
julia> using BenchmarkTools

julia> @btime UpdateCellList!(\$x,\$y,\$box,\$cl,\$aux)
  715.661 μs (41 allocations: 3.88 KiB)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   50000 particles in the reference vector.
   7414 cells with real particles of target vector.
   
julia> @btime UpdateCellList!(\$x,\$y,\$box,\$cl,\$aux,parallel=false)
   13.042 ms (0 allocations: 0 bytes)
 CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
    50000 particles in the reference vector.
    15031 cells with real particles of target vector.
 
```

"""
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair,
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {UnitCellType,N,T}
    if !cl_pair.swap 
        target = UpdateCellList!(x,box,cl_pair.target,aux,parallel=parallel)
    else
        target = UpdateCellList!(y,box,cl_pair.target,aux,parallel=parallel)
    end
    cl_pair = CellListPair(ref=cl_pair.ref,target=target,swap=cl_pair.swap)
    return cl_pair
end

"""

```
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair,
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {UnitCellType,N,T}
```

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N,T},
    cl_pair::CellListPair,
    aux::AuxThreaded{N,T};
    parallel::Bool=true
) where {UnitCellType,N,T}
    @assert size(x,1) == N "First dimension of input matrix must be $N"
    @assert size(y,1) == N "First dimension of input matrix must be $N"
    x_re = reinterpret(reshape, SVector{N,T}, x)
    y_re = reinterpret(reshape, SVector{N,T}, x)
    return UpdateCellList!(x_re,y_re,box,cl_pair,aux,parallel=parallel)
end

"""

```
particles_per_cell(cl::CellList,box::Box)
```

Returns the average number of particles per computing cell.

"""
particles_per_cell(cl::CellList,box::Box) = cl.ncp[1] / prod(box.nc)


