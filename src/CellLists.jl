#
# This file contains all structre types and functions necessary for building
# the CellList and CellListPair structures.
#

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort. 

"""
struct ParticleWithIndex{N,T}
    index::Int
    real::Bool
    coordinates::SVector{N,T}
end
Base.zero(::Type{ParticleWithIndex{N,T}}) where {N,T} =
    ParticleWithIndex{N,T}(0, false, zeros(SVector{N,T}))

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure to define the number of batches used in the parallel splitting of the calculations
of the cell list construction and of the `map_pairwise` computation. It is initialized with
a standard heuristic that returns at most the number of threads, but may return a smaller
number if the system is small. The two parameters can be tunned for optimal performance of each
step of the calculation (cell list construction and mapping of interactions). The construction
of the cell lists require a larger number of particles for threading to be effective, Thus
by default the system size that allows multi-threading is greater for this part of the calculation.  

"""
struct NumberOfBatches
    build_cell_lists::Int
    map_computation::Int
end
NumberOfBatches(t::Tuple{Int,Int}) = NumberOfBatches(t[1], t[2])
Base.zero(::Type{NumberOfBatches}) = NumberOfBatches(0, 0)
Base.iszero(x::NumberOfBatches) = (iszero(x.build_cell_lists) && iszero(x.map_computation))
function Base.show(io::IO, ::MIME"text/plain", nbatches::NumberOfBatches)
    _println(io, "  Number of batches for cell list construction: $(nbatches.build_cell_lists)")
    _print(io, "  Number of batches for function mapping: $(nbatches.map_computation)")
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

This structure contains the cell linear index and the information 
about if this cell is in the border of the box (such that its 
neighboring cells need to be wrapped) 

"""
Base.@kwdef struct Cell{N,T}
    linear_index::Int = 0
    cartesian_index::CartesianIndex{N} = CartesianIndex{N}(ntuple(i -> 0, N))
    center::SVector{N,T} = zeros(SVector{N,T})
    contains_real::Bool = false
    n_particles::Int = 0
    particles::Vector{ParticleWithIndex{N,T}} = Vector{ParticleWithIndex{N,T}}(undef, 0)
end
function Cell{N,T}(cartesian_index::CartesianIndex, box::Box; sizehint::Int=0) where {N,T}
    return Cell{N,T}(
        linear_index=cell_linear_index(box.nc, cartesian_index),
        cartesian_index=cartesian_index,
        center=cell_center(cartesian_index, box),
        particles=Vector{ParticleWithIndex{N,T}}(undef, sizehint)
    )
end

function copy_cell(cell::Cell{N,T}) where {N,T}
    return Cell{N,T}(
        linear_index=cell.linear_index,
        cartesian_index=cell.cartesian_index,
        center=cell.center,
        contains_real=cell.contains_real,
        n_particles=cell.n_particles,
        particles=ParticleWithIndex{N,T}[p for p in cell.particles]
    )
end


"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Auxiliary structure to contain projected particles. Types of 
scalars are chosen such that with a `SVector{3,Float64}` the
complete struct has 32bytes.

"""
Base.@kwdef struct ProjectedParticle{N,T}
    index::Int
    xproj::T
    coordinates::SVector{N,T}
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

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
    cell_indices::Vector{Int} = zeros(Int, number_of_cells)
    " Auxiliary array that contains the indexes in the cells with real particles. "
    cell_indices_real::Vector{Int} = zeros(Int, 0)
    " Vector containing cell lists of cells with particles. "
    cells::Vector{Cell{N,T}} = Cell{N,T}[]
    " Number of batches for the parallel calculations. "
    nbatches::NumberOfBatches = zero(NumberOfBatches)
    " Auxiliar array to store projected particles. "
    projected_particles::Vector{Vector{ProjectedParticle{N,T}}} =
        Vector{Vector{ProjectedParticle{N,T}}}(undef, 0)
end
function Base.show(io::IO, ::MIME"text/plain", cl::CellList)
    _println(io, typeof(cl))
    _println(io, "  $(cl.n_real_particles) real particles.")
    _println(io, "  $(cl.n_cells_with_real_particles) cells with real particles.")
    _print(io, "  $(cl.n_particles) particles in computing box, including images.")
end

"""

$(INTERNAL)

Structures to control dispatch on swapped vs. not swapped cell list pairs.

"""
struct Swapped end
struct NotSwapped end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that will cointain the cell lists of two independent sets of
particles for cross-computation of interactions

"""
struct CellListPair{V,N,T,Swap}
    ref::V
    target::CellList{N,T}
end
CellListPair(ref::V, target::CellList{N,T}, ::Swap) where {V,N,T,Swap} =
    CellListPair{V,N,T,Swap}(ref, target)

function Base.show(io::IO, ::MIME"text/plain", cl::CellListPair)
    _print(io, typeof(cl), "\n")
    _print(io, "   $(length(cl.ref)) particles in the reference vector.\n")
    _print(io, "   $(cl.target.n_cells_with_real_particles) cells with real particles of target vector.")
end

"""
    set_number_of_batches!(cl,nbatches::Tuple{Int,Int}=(0,0);parallel=true)  

$(INTERNAL)

# Extended help

Functions that set the default number of batches for the construction of the cell lists, 
and mapping computations. This is of course heuristic, and may not be the best choice for
every problem. See the parameter `nbatches` of the construction of the cell lists for 
tunning this.

"""
function set_number_of_batches!(cl::CellList{N,T}, nbatches::Tuple{Int,Int}=(0, 0); parallel=true) where {N,T}
    if parallel
        nbatches = NumberOfBatches(nbatches)
    else
        if nbatches != (0, 0) && nbatches != (1, 1)
            @warn "nbatches set to $nbatches, but parallel is set to false, implying nbatches == (1, 1)"
        end
        nbatches = NumberOfBatches((1, 1))
    end
    if nbatches.build_cell_lists < 1
        n1 = _nbatches_build_cell_lists(cl.n_real_particles)
    else
        n1 = nbatches.build_cell_lists
    end
    if nbatches.map_computation < 1
        n2 = _nbatches_map_computation(cl.n_real_particles)
    else
        n2 = nbatches.map_computation
    end
    nbatches = NumberOfBatches(n1, n2)
    @set! cl.nbatches = nbatches
    for _ in 1:cl.nbatches.map_computation
        push!(cl.projected_particles, Vector{ProjectedParticle{N,T}}(undef, 0))
    end
    return cl
end
# Heuristic choices for the number of batches, for an atomic system
_nbatches_build_cell_lists(n::Int) = max(1, min(n, min(8, nthreads())))
_nbatches_map_computation(n::Int) = max(1, min(n, min(floor(Int, 2^(log10(n) + 1)), nthreads())))

function set_number_of_batches!(
    cl::CellListPair{V,N,T,Swap},
    nbatches::Tuple{Int,Int}=(0, 0);
    parallel=true
) where {V,N,T,Swap}
    if parallel
        nbatches = NumberOfBatches(nbatches)
    else
        if nbatches != (0, 0) && nbatches != (1, 1)
            @warn "nbatches set to $nbatches, but parallel is set to false, implying nbatches == (1, 1)"
        end
        nbatches = NumberOfBatches((1, 1))
    end
    if nbatches.build_cell_lists < 1
        n1 = _nbatches_build_cell_lists(cl.target.n_real_particles)
    else
        n1 = nbatches.build_cell_lists
    end
    if nbatches.map_computation < 1
        n2 = _nbatches_map_computation(length(cl.ref))
    else
        n2 = nbatches.map_computation
    end
    nbatches = NumberOfBatches(n1, n2)
    target = cl.target
    @set! target.nbatches = nbatches
    return CellListPair{V,N,T,Swap}(cl.ref, target)
end

"""
    nbatches(cl)

Returns the number of batches for parallel processing that will be used in the pairwise function mappings associated to cell list `cl`. 
It returns the `cl.nbatches.map_computation` value. This function is important because it must be used to set the number of copies
of custom preallocated output arrays.

A second argument can be provided, which may be `:map` or `:build`, in which case the function returns either the number of batches used 
for pairwise mapping or for the construction of the cell lists. Since this second value is internal and does not affect the interface, 
it can be usually ignored. 

## Example

```julia-repl
julia> x = rand(3,1000); box = Box([1,1,1],0.1);

julia> cl = CellList(x,box,nbatches=(2,16));

julia> nbatches(cl)
16

julia> nbatches(cl,:map)
16

julia> nbatches(cl,:build)
2
```

"""
nbatches(cl::CellList) = cl.nbatches.map_computation
function nbatches(cl::CellList, s::Symbol)
    s == :map_computation || s == :map && return cl.nbatches.map_computation
    s == :build_cell_lists || s == :build && return cl.nbatches.build_cell_lists
end
nbatches(cl::CellListPair) = nbatches(cl.target)
nbatches(cl::CellListPair, s::Symbol) = nbatches(cl.target, s)

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Auxiliary structure to carry threaded lists and ranges of particles to 
be considered by each thread on parallel construction. 

"""
@with_kw struct AuxThreaded{N,T}
    particles_per_batch::Int
    idxs::Vector{UnitRange{Int}} = Vector{UnitRange{Int}}(undef, 0)
    lists::Vector{CellList{N,T}} = Vector{CellList{N,T}}(undef, 0)
end
function Base.show(io::IO, ::MIME"text/plain", aux::AuxThreaded)
    _println(io, typeof(aux))
    _print(io, " Auxiliary arrays for nbatches = ", length(aux.lists))
end

"""
    AuxThreaded(cl::CellList{N,T}) where {N,T}

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
function AuxThreaded(cl::CellList{N,T}; particles_per_batch=10_000) where {N,T}
    nbatches = cl.nbatches.build_cell_lists
    aux = AuxThreaded{N,T}(
        particles_per_batch=particles_per_batch,
        idxs=Vector{UnitRange{Int}}(undef, nbatches),
        lists=Vector{CellList{N,T}}(undef, nbatches)
    )
    # If the calculation is not parallel, no need to initialize this
    nbatches == 1 && return aux
    for ibatch in eachindex(aux.lists)
        cl_batch = CellList{N,T}(
            n_real_particles=particles_per_batch, # this is reset before filling, in UpdateCellList!
            number_of_cells=cl.number_of_cells,
        )
        aux.lists[ibatch] = cl_batch
    end
    # Set indices of the atoms that will be considered by each thread
    # these indices may be updated by an update of cell lists, if the number
    # of particles change.
    set_idxs!(aux.idxs, cl.n_real_particles, nbatches)
    return aux
end


"""
    set_idxs!(idxs, n_particles, nbatches)

$(INTERNAL)

# Extended help

Sets the indexes of the particles that will be considered for each batch in parallel runs.
Modifies the `idxs` array of ranges, which is usually the `aux.idxs` array of the the 
corresponding `AuxThreaded` structure.

"""
function set_idxs!(idxs, n_particles, nbatches)
    length(idxs) == nbatches || throw(ArgumentError("Modifying `nbatches` requires an explicit update of the AuxThreaded auxiliary array."))
    nrem = n_particles % nbatches
    nperthread = (n_particles - nrem) ÷ nbatches
    first = 1
    for ibatch in eachindex(idxs)
        nx = nperthread
        if ibatch <= nrem
            nx += 1
        end
        idxs[ibatch] = first:(first-1)+nx
        first += nx
    end
    return nothing
end

"""
    AuxThreaded(cl::CellListPair{N,T}) where {N,T}

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
AuxThreaded(cl_pair::CellListPair; particles_per_batch=10_000) =
    AuxThreaded(cl_pair.target, particles_per_batch=particles_per_batch)

"""
    CellList(
        x::AbstractVector{AbstractVector},
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0)
    ) where {UnitCellType,N,T} 

Function that will initialize a `CellList` structure from scratch, given a vector
or particle coordinates (a vector of vectors, typically of static vectors) 
and a `Box`, which contain the size ofthe system, cutoff, etc. Except for small
systems, the number of parallel batches is equal to the number of threads, but it can
be tunned for optimal performance in some cases.

## Example

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
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0)
) where {UnitCellType,N,T}
    cl = CellList{N,T}(
        n_real_particles=length(x),
        number_of_cells=prod(box.nc),
    )
    cl = set_number_of_batches!(cl, nbatches, parallel=parallel)
    return UpdateCellList!(x, box, cl, parallel=parallel)
end

"""
    CellList(x::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T} 

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

"""
function CellList(x::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return CellList(x_re, box; kargs...)
end

"""
    reset!(cl::CellList{N,T},box,n_real_particles) where{N,T}

$(INTERNAL)

# Extended help

Resets a cell list, by setting everything to zero, but retaining
the allocated `particles` and `projected_particles` vectors.

"""
function reset!(cl::CellList{N,T}, box, n_real_particles) where {N,T}
    new_number_of_cells = prod(box.nc)
    if new_number_of_cells > cl.number_of_cells
        resize!(cl.cell_indices, new_number_of_cells)
    end
    for i in eachindex(cl.cells)
        cl.cells[i] = Cell{N,T}(particles=cl.cells[i].particles)
    end
    @. cl.cell_indices = 0
    @. cl.cell_indices_real = 0
    cl = CellList{N,T}(
        n_real_particles=n_real_particles,
        n_particles=0,
        number_of_cells=new_number_of_cells,
        n_cells_with_real_particles=0,
        n_cells_with_particles=0,
        cell_indices=cl.cell_indices,
        cell_indices_real=cl.cell_indices_real,
        cells=cl.cells,
        nbatches=cl.nbatches,
        projected_particles=cl.projected_particles
    )
    return cl
end

"""
    CellList(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0),
        autoswap::Bool=true
    ) where {UnitCellType,N,T} 

Function that will initialize a `CellListPair` structure from scratch, given two vectors
of particle coordinates and a `Box`, which contain the size of the system, cutoff, etc.
By default, the cell list will be constructed for smallest vector, but this is not always
the optimal choice. Using `autoswap=false` the cell list is constructed for the second (`y`)

## Example

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
    nbatches::Tuple{Int,Int}=(0, 0),
    autoswap=true
) where {UnitCellType,N,T}
    if !autoswap || length(x) >= length(y)
        ref = [SVector{N,T}(ntuple(i -> el[i], N)) for el in x]
        target = CellList(y, box, parallel=parallel)
        swap = NotSwapped()
    else
        ref = [SVector{N,T}(ntuple(i -> el[i], N)) for el in y]
        target = CellList(x, box, parallel=parallel)
        swap = Swapped()
    end
    cl_pair = CellListPair(ref, target, swap)
    cl_pair = set_number_of_batches!(cl_pair, nbatches, parallel=parallel)
    return cl_pair
end

"""
    CellList(x::AbstractMatrix, y::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T} 

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function CellList(x::AbstractMatrix, y::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    CellList(x_re, y_re, box; kargs...)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        box::Box,
        cl:CellList,
        parallel=true
    ) 

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
    cl::CellList;
    parallel::Bool=true
)
    if parallel
        aux = AuxThreaded(cl)
        return UpdateCellList!(x, box, cl, aux, parallel=parallel)
    else
        return UpdateCellList!(x, box, cl, nothing, parallel=parallel)
    end
end

"""
    function UpdateCellList!(
        x::AbstractMatrix,
        box::Box,
        cl::CellList{N,T};
        parallel::Bool=true
    ) where {N,T}

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
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return UpdateCellList!(x_re, box, cl, parallel=parallel)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        box::Box,
        cl::CellList{N,T},
        aux::Union{Nothing,AuxThreaded{N,T}};
        parallel::Bool=true
    ) where {N,T}

Function that updates the cell list `cl` new coordinates `x` and possibly a new
box `box`, and receives a preallocated `aux` structure of auxiliary vectors for
threaded cell list construction. Given a preallocated `aux` vector, allocations in
this function should be minimal, only associated with the spawning threads, or
to expansion of the cell lists if the number of cells or number of particles 
increased. 

## Example

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
    aux::Union{Nothing,AuxThreaded{N,T}};
    parallel::Bool=true
) where {N,T}

    # Provide a better error message if the unit cell dimension does not match the dimension of the positions.
    if length(x) > 0 && (length(x[begin]) != size(box.input_unit_cell.matrix, 1))
        n1 = length(x[begin])
        n2 = size(box.input_unit_cell.matrix, 1)
        throw(DimensionMismatch("Positions have dimension $n1, but the unit cell has dimension $n2."))
    end

    # Add particles to cell list
    nbatches = cl.nbatches.build_cell_lists
    if !parallel || nbatches == 1
        cl = reset!(cl, box, length(x))
        cl = add_particles!(x, box, 0, cl)
    else
        # Reset cell list
        cl = reset!(cl, box, 0)
        # Update the aux.idxs ranges, for if the number of particles changed
        set_idxs!(aux.idxs, length(x), nbatches)
        @sync for ibatch in eachindex(aux.idxs, aux.lists)
            Threads.@spawn begin
                prange = aux.idxs[ibatch]
                aux.lists[ibatch] = reset!(aux.lists[ibatch], box, length(prange))
                xt = @view(x[prange])
                aux.lists[ibatch] = add_particles!(xt, box, prange[begin] - 1, aux.lists[ibatch])
            end
        end
        for ibatch in eachindex(aux.lists)
            cl = merge_cell_lists!(cl, aux.lists[ibatch])
        end
    end

    # allocate, or update the auxiliary projected_particles arrays
    maxnp = 0
    for i in 1:cl.n_cells_with_particles
        maxnp = max(maxnp, cl.cells[i].n_particles)
    end
    for i in eachindex(cl.projected_particles)
        if maxnp > length(cl.projected_particles[i])
            resize!(cl.projected_particles[i], maxnp)
        end
    end

    return cl
end

"""
    UpdateCellList!(
        x::AbstractMatrix,
        box::Box,
        cl::CellList{N,T},
        aux::Union{Nothing,AuxThreaded{N,T}};
        parallel::Bool=true
    ) where {N,T}

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T},
    aux::Union{Nothing,AuxThreaded{N,T}};
    parallel::Bool=true
) where {N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return UpdateCellList!(x_re, box, cl, aux, parallel=parallel)
end

"""
    add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}

$(INTERNAL)

# Extended help

Add all particles in vector `x` to the cell list `cl`. `ishift` is the shift in particle
index, meaning that particle `i` of vector `x` corresponds to the particle with original
index `i+ishift`. The shift is used to construct cell lists from fractions of the original
set of particles in parallel list construction.  

"""
function add_particles!(x, box, ishift, cl::CellList{N,T}) where {N,T}
    for ip in eachindex(x)
        xp = x[ip]
        # This converts the coordinates to static arrays, if necessary
        p = SVector{N,T}(ntuple(i -> xp[i], N))
        p = box.rotation * wrap_to_first(p, box.input_unit_cell.matrix)
        cl = add_particle_to_celllist!(ishift + ip, p, box, cl) # add real particle
        cl = replicate_particle!(ishift + ip, p, box, cl) # add virtual particles to border cells
    end
    return cl
end

# define method for the ParticleWithIndex type
out_of_computing_box(p::ParticleWithIndex, box::Box) = out_of_computing_box(p.coordinates, box)

"""
    copydata!(cell1::Cell,cell2::Cell)

$(INTERNAL)

# Extended help

Copies the data from `cell2` to `cell1`, meaning that particles are
copied element-wise from `cell2` to `cell1`, with the `particles` array
of `cell1` being resized (increased) if necessary.

"""
function copydata!(cell1::Cell, cell2::Cell)
    @set! cell1.linear_index = cell2.linear_index
    @set! cell1.cartesian_index = cell2.cartesian_index
    @set! cell1.center = cell2.center
    @set! cell1.n_particles = cell2.n_particles
    @set! cell1.contains_real = cell2.contains_real
    for ip in 1:cell2.n_particles
        p = cell2.particles[ip]
        if ip > length(cell1.particles)
            push!(cell1.particles, p)
        else
            cell1.particles[ip] = p
        end
    end
    return cell1
end

"""
    append_particles!(cell1::Cell,cell2::Cell)

$(INTERNAL)

# Extended help

Add the particles of `cell2` to `cell1`, updating the cell data and, if necessary,
resizing (increasing) the `particles` array of `cell1`

"""
function append_particles!(cell1::Cell, cell2::Cell)
    if cell2.contains_real
        @set! cell1.contains_real = true
    end
    n_particles_old = cell1.n_particles
    @set! cell1.n_particles += cell2.n_particles
    if cell1.n_particles > length(cell1.particles)
        resize!(cell1.particles, cell1.n_particles)
    end
    for ip in 1:cell2.n_particles
        cell1.particles[n_particles_old+ip] = cell2.particles[ip]
    end
    return cell1
end

"""
    merge_cell_lists!(cl::CellList,aux::CellList)

$(INTERNAL)

# Extended help

Merges an auxiliary `aux` cell list to `cl`, and returns the modified `cl`. Used to
merge cell lists computed in parallel threads.

"""
function merge_cell_lists!(cl::CellList, aux::CellList)
    # One should never get here if the lists do not share the same # computing box
    if cl.number_of_cells != aux.number_of_cells
        error("cell lists must have the same number of cells to be merged.")
    end
    @set! cl.n_particles += aux.n_particles
    @set! cl.n_real_particles += aux.n_real_particles
    for icell in 1:aux.n_cells_with_particles
        aux_cell = aux.cells[icell]
        linear_index = aux_cell.linear_index
        cell_index = cl.cell_indices[linear_index]
        # If cell was yet not initialized in merge, push it to the list
        if cell_index == 0
            @set! cl.n_cells_with_particles += 1
            cell_index = cl.n_cells_with_particles
            cl.cell_indices[linear_index] = cell_index
            if cell_index > length(cl.cells)
                push!(cl.cells, copy_cell(aux_cell))
            else
                cl.cells[cell_index] = copydata!(cl.cells[cell_index], aux_cell)
            end
            if aux_cell.contains_real
                @set! cl.n_cells_with_real_particles += 1
                if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                    push!(cl.cell_indices_real, cell_index)
                else
                    cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
                end
            end
            # Append particles to initialized cells
        else
            # If the previous cell didn't contain real particles, but the current one
            # does, update the list information 
            if !cl.cells[cell_index].contains_real && aux_cell.contains_real
                @set! cl.n_cells_with_real_particles += 1
                if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
                    push!(cl.cell_indices_real, cell_index)
                else
                    cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
                end
            end
            cl.cells[cell_index] = append_particles!(cl.cells[cell_index], aux_cell)
        end
    end
    return cl
end

"""
    add_particle_to_celllist!(
        ip,
        x::SVector{N,T},
        box,
        cl::CellList{N,T};
        real_particle::Bool=true
    ) where {N,T}

$(INTERNAL)

# Extended help

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
    cartesian_index = particle_cell(x, box)
    linear_index = cell_linear_index(box.nc, cartesian_index)

    # Check if this is the first particle of this cell, if it is,
    # initialize a new cell, or reset a previously allocated one
    cell_index = cell_indices[linear_index]

    if cell_index == 0
        n_cells_with_particles += 1
        cell_index = n_cells_with_particles
        cell_indices[linear_index] = cell_index
        if cell_index > length(cells)
            particles_sizehint = cl.n_real_particles ÷ prod(box.nc)
            push!(cells, Cell{N,T}(cartesian_index, box, sizehint=particles_sizehint))
        else
            cell = cells[cell_index]
            @set! cell.linear_index = linear_index
            @set! cell.cartesian_index = cartesian_index
            @set! cell.center = cell_center(cell.cartesian_index, box)
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
            push!(cell_indices_real, cell_index)
        else
            cell_indices_real[n_cells_with_real_particles] = cell_index
        end
    end

    #
    # Add particle to cell list
    #
    p = ParticleWithIndex(ip, real_particle, x)
    if cell.n_particles > length(cell.particles)
        push!(cell.particles, p)
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
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box,
        cl:CellListPair,
        parallel=true
    )

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
    box::Box,
    cl_pair::CellListPair;
    parallel::Bool=true
)
    if parallel
        aux = AuxThreaded(cl_pair)
        return UpdateCellList!(x, y, box, cl_pair, aux, parallel=parallel)
    else
        return UpdateCellList!(x, y, box, cl_pair, nothing, parallel=parallel)
    end
end

"""
    UpdateCellList!(
        x::AbstractMatrix,
        y::AbstractMatrix,
        box::Box{UnitCellType,N},
        cl_pair::CellListPair;
        parallel::Bool=true
    ) where {UnitCellType,N}

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N},
    cl_pair::CellListPair;
    parallel::Bool=true
) where {UnitCellType,N}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return UpdateCellList!(x_re, y_re, box, cl_pair, parallel=parallel)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box,
        cl_pair::CellListPair,
        aux::Union{Nothing,AuxThreaded};
        parallel::Bool=true
    )

This function will update the `cl_pair` structure that contains the cell lists
for disjoint sets of particles. It receives the preallocated `aux` structure to
avoid reallocating auxiliary arrays necessary for the threaded construct of the
lists. 

## Example

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
    box::Box,
    cl_pair::CellListPair{V,N,T,Swap},
    aux::Union{Nothing,AuxThreaded};
    parallel::Bool=true
) where {V,N,T,Swap<:NotSwapped}
    ref = x
    target = UpdateCellList!(y, box, cl_pair.target, aux, parallel=parallel)
    cl_pair = _update_CellListPair!(ref, target, cl_pair)
    return cl_pair
end
# Swapped vectors version
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box,
    cl_pair::CellListPair{V,N,T,Swap},
    aux::Union{Nothing,AuxThreaded};
    parallel::Bool=true
) where {V,N,T,Swap<:Swapped}
    ref = y
    target = UpdateCellList!(x, box, cl_pair.target, aux, parallel=parallel)
    cl_pair = _update_CellListPair!(ref, target, cl_pair)
    return cl_pair
end

# Function barrier that was required to avoid the `Swap` type to cause some instability
function _update_CellListPair!(ref, target, cl_pair::CellListPair{V,N,T,Swap}) where {V,N,T,Swap}
    # This resizing will fail if the data was input as a (N,M) matrix, because resizing
    # is not implemented for reinterpreted arrays. 
    if length(ref) != length(cl_pair.ref)
        resize!(cl_pair.ref, length(ref))
    end
    for i in eachindex(ref, cl_pair.ref)
        @inbounds cl_pair.ref[i] = SVector{N,T}(ntuple(j -> ref[i][j], N))
    end
    cl_pair = CellListPair{V,N,T,Swap}(cl_pair.ref, target)
    return cl_pair
end

"""
    UpdateCellList!(
        x::AbstractMatrix,
        y::AbstractMatrix,
        box::Box,
        cl_pair::CellListPair,
        aux::Union{Nothing,AuxThreaded};
        parallel::Bool=true
    ) where {UnitCellType,N}

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

"""
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N},
    cl_pair::CellListPair,
    aux::Union{Nothing,AuxThreaded};
    parallel::Bool=true
) where {UnitCellType,N}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return UpdateCellList!(x_re, y_re, box, cl_pair, aux, parallel=parallel)
end

"""
    particles_per_cell(cl)

$(INTERNAL)

Returns the average number of real particles per computing cell.

"""
particles_per_cell(cl::CellList) = cl.n_real_particles / cl.number_of_cells
particles_per_cell(cl::CellListPair) = particles_per_cell(cl.target)


