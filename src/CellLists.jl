#
# This file contains all structure types and functions necessary for building
# the CellList and CellListPair structures.
#

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Copies particle coordinates and associated index, to build contiguous particle lists
in memory when building the cell lists. This strategy duplicates the particle coordinates
data, but is probably worth the effort. 

=#
struct ParticleWithIndex{N,T}
    index::Int
    real::Bool
    coordinates::SVector{N,T}
end
Base.zero(::Type{ParticleWithIndex{N,T}}) where {N,T} =
    ParticleWithIndex{N,T}(0, false, zeros(SVector{N,T}))

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Structure to define the number of batches used in the parallel splitting of the calculations
of the cell list construction and of the `map_pairwise` computation. It is initialized with
a standard heuristic that returns at most the number of threads, but may return a smaller
number if the system is small. The two parameters can be tunned for optimal performance of each
step of the calculation (cell list construction and mapping of interactions). The construction
of the cell lists require a larger number of particles for threading to be effective, Thus
by default the system size that allows multi-threading is greater for this part of the calculation.  

=#
struct NumberOfBatches
    build_cell_lists::Tuple{Bool,Int}
    map_computation::Tuple{Bool,Int}
end
NumberOfBatches(auto::Tuple{Bool,Bool}, n::Tuple{Int,Int}) = 
    NumberOfBatches((first(auto), first(n)), (last(auto), last(n)))
Base.zero(::Type{NumberOfBatches}) = NumberOfBatches((true,0), (true,0))
function Base.show(io::IO, ::MIME"text/plain", nbatches::NumberOfBatches)
    _println(io,"  Number of batches for cell list construction: $(last(nbatches.build_cell_lists)) (auto-update: $(first(nbatches.build_cell_lists)))")
      _print(io,"  Number of batches for function mapping: $(last(nbatches.map_computation)) (auto-update: $(first(nbatches.map_computation)))")
end

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

This structure contains the cell linear index and the information 
about if this cell is in the border of the box (such that its 
neighboring cells need to be wrapped) 

=#
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

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Auxiliary structure to contain projected particles. Types of 
scalars are chosen such that with a `SVector{3,Float64}` the
complete struct has 32bytes.

=#
Base.@kwdef struct ProjectedParticle{N,T}
    index::Int
    xproj::T
    coordinates::SVector{N,T}
end

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Structure that contains the cell lists information.

=#
Base.@kwdef mutable struct CellList{N,T}
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

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Structure that will contains the cell lists of two independent sets of
particles for cross-computation of interactions

=#
struct CellListPair{N,T}
    small_set::CellList{N,T}
    large_set::CellList{N,T}
    swap::Bool
end

function Base.show(io::IO, ::MIME"text/plain", cl::CellListPair)
    _print(io, typeof(cl), "\n")
    _print(io, "   $(cl.small_set.n_cells_with_real_particles) cells with real particles of the smallest set.\n")
    _print(io, "   $(cl.large_set.n_cells_with_real_particles) cells with real particles of the largest set.")
end

#=
    update_number_of_batches!(cl, nbatches::NumberOfBatches; parallel=true)  

Set the default number of batches for the construction of the cell lists, 
and mapping computations. This is of course heuristic, and may not be the best choice for
every problem. See the parameter `nbatches` of the construction of the cell lists for 
tunning this.

=#
function update_number_of_batches!(cl::CellList{N,T}, _nbatches=cl.nbatches; parallel=true) where {N,T}
    auto = (first(_nbatches.build_cell_lists), first(_nbatches.map_computation))
    n1 = last(_nbatches.build_cell_lists)
    n2 = last(_nbatches.map_computation)
    if !parallel
        if !all(auto) && (n1, n2) != (1, 1)
            @warn begin
                """\n
                WARNING: nbatches set to ($n1,$n2), but parallel is set to false, implying nbatches == (1, 1)
    
                """
            end _file=nothing _line=nothing
        end
        nbatches = NumberOfBatches(auto, (1, 1))
    else # Heuristic choices
        if first(auto)
            n1 = _nbatches_build_cell_lists(cl.n_real_particles)
        end
        if last(auto)
            n2 = _nbatches_map_computation(cl.n_real_particles)
        end
        nbatches = NumberOfBatches(auto, (n1, n2))
    end
    for _ in 1:n2
        push!(cl.projected_particles, Vector{ProjectedParticle{N,T}}(undef, 0))
    end
    cl.nbatches = nbatches
    return cl
end

# Heuristic choices for the number of batches, for an atomic system
_nbatches_build_cell_lists(n::Int) = max(1, min(n, min(8, nthreads())))
_nbatches_map_computation(n::Int) = max(1, min(n, min(floor(Int, 2^(log10(n) + 1)), nthreads())))

function update_number_of_batches!(cl::CellListPair{N,T}, _nbatches::NumberOfBatches; parallel=true) where {N,T}
    large_set = update_number_of_batches!(cl.large_set, _nbatches; parallel)
    return CellListPair{N,T}(
        update_number_of_batches!(cl.small_set, large_set.nbatches; parallel),
        large_set,
        cl.swap,
    )
end

#
# Functions for initialization of the batches, called from the API functions. Receives a 
# tuple with the number of batches for the construction of the cell lists and the mapping.
# If the number of batches are smaller than 1, the function will set the number of batches
# in automatic mode.
#
function set_number_of_batches!(cl::Union{CellList,CellListPair}, _nbatches::Tuple{Int,Int}; parallel=true)
    auto = _nbatches .<= 0
    nbatches = NumberOfBatches((first(auto), first(_nbatches)), (last(auto), last(_nbatches)))
    return update_number_of_batches!(cl, nbatches; parallel)
end

"""
    nbatches(cl::CellList)
    nbatches(cl::CellListPair)
    nbatches(system::AbstractParticleSystem)

Returns the number of batches for parallel processing that will be used. Returns a tuple, where 
the first element is the number of batches for the construction of the cell lists, and the second
element is the number of batches for the pairwise mapping.

A second argument can be provided, which may be `:map` or `:build`, in which case the function returns either the number of batches used 
for pairwise mapping or for the construction of the cell lists. Since this second value is internal and does not affect the interface, 
it can be usually ignored. 

## Example

```jldoctest
julia> using CellListMap

julia> x = rand(3,1000); box = Box([1,1,1],0.1);

julia> cl = CellList(x, box, nbatches=(2,16));

julia> nbatches(cl)
(2, 16)

julia> nbatches(cl,:build)
2

julia> nbatches(cl,:map)
16

```

"""
function nbatches(cl::CellList, s::Symbol)
    s == :map_computation || s == :map && return last(cl.nbatches.map_computation)
    s == :build_cell_lists || s == :build && return last(cl.nbatches.build_cell_lists)
end
nbatches(cl::CellList) = (nbatches(cl, :build), nbatches(cl, :map))
nbatches(cl::CellListPair) = nbatches(cl.large_set)
nbatches(cl::CellListPair, s::Symbol) = nbatches(cl.large_set, s)

@testitem "output of nbatches" begin
    using CellListMap, StaticArrays
    x = rand(3,100)
    y = rand(3,10)
    box = Box([1,1,1],0.1)
    cl = CellList(x,box; nbatches=(2,4))
    @test nbatches(cl) == (2, 4)
    @test nbatches(cl, :build) == 2
    @test nbatches(cl, :map) == 4
    cl = CellList(x,y,box; nbatches=(2,4))
    @test nbatches(cl) == (2, 4)
    @test nbatches(cl, :build) == 2
    @test nbatches(cl, :map) == 4
    cl = CellList(x,y,box)
    @test nbatches(cl.small_set) == nbatches(cl.large_set)
    cl = CellList(x,box; nbatches=(2,4), parallel=false)
    @test nbatches(cl) == (1,1)
    # The automatic set of number of batches for this small system:
    if Threads.nthreads() == 10
        x = rand(SVector{3,Float64},10)
        sys = ParticleSystem(positions=x, cutoff=0.1, unitcell=[1,1,1], output=0.0)
        @test nbatches(sys) == (8,4)
        x = rand(SVector{3,Float64},10000)
        resize!(sys.xpositions, length(x))
        sys.xpositions .= x
        map_pairwise((_, _, _, _, d2, out) -> out += d2, sys)
        @test nbatches(sys) == (8,10)
    else
        @warn "Test not run because it is invalid for this number of threads"
    end
end

@testitem "automatic nbataches update" begin
    using CellListMap, StaticArrays
    x, box = CellListMap.xatomic(5000)

end

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Auxiliary structure to carry threaded lists and ranges of particles to 
be considered by each thread on parallel construction. 

=#
@with_kw struct AuxThreaded{N,T}
    idxs::Vector{UnitRange{Int}} = Vector{UnitRange{Int}}(undef, 0)
    lists::Vector{CellList{N,T}} = Vector{CellList{N,T}}(undef, 0)
end
function Base.show(io::IO, ::MIME"text/plain", aux::AuxThreaded)
    _println(io, typeof(aux))
    _print(io, " Auxiliary arrays for nbatches = ", length(aux.lists))
end

@with_kw struct AuxThreadedPair{N,T}
    small_set::AuxThreaded{N,T}
    large_set::AuxThreaded{N,T}
end
function Base.show(io::IO, ::MIME"text/plain", aux::AuxThreadedPair)
    _println(io, typeof(aux))
    _print(io, " Auxiliary arrays for nbatches = ", length(aux.small_set.lists))
end

"""
    AuxThreaded(cl::CellList{N,T}) where {N,T}

Constructor for the `AuxThreaded` type, to be passed to `UpdateCellList!` for in-place 
update of cell lists. 

## Example
```julia-repl
julia> using CellListMap

julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(3) for i in 1:10000 ];

julia> cl = CellList(x,box);

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> UpdateCellList!(x,box,cl,aux)
CellList{3, Float64}
  100000 real particles.
  31190 cells with real particles.
  1134378 particles in computing box, including images.

```
"""
function AuxThreaded(cl::CellList{N,T}) where {N,T}
    _nbatches = nbatches(cl, :build)
    aux = AuxThreaded{N,T}(
        idxs=Vector{UnitRange{Int}}(undef, _nbatches),
        lists=Vector{CellList{N,T}}(undef, _nbatches)
    )
    # If the calculation is not parallel, no need to initialize this
    _nbatches == 1 && return aux
    @sync for ibatch in eachindex(aux.lists)
        @spawn begin
            cl_batch = CellList{N,T}(number_of_cells=cl.number_of_cells)
            aux.lists[ibatch] = cl_batch
        end
    end
    # Set indices of the atoms that will be considered by each thread
    # these indices may be updated by an update of cell lists, if the number
    # of particles change.
    set_idxs!(aux.idxs, cl.n_real_particles, _nbatches)
    return aux
end


#=
    set_idxs!(idxs, n_particles, nbatches)

# Extended help

Sets the indexes of the particles that will be considered for each batch in parallel runs.
Modifies the `idxs` array of ranges, which is usually the `aux.idxs` array of the the 
corresponding `AuxThreaded` structure.

=#
function set_idxs!(idxs, n_particles, nbatches)
    if length(idxs) != nbatches 
        throw(ArgumentError("""\n 
            Modifying `nbatches` requires an explicit update of the AuxThreaded auxiliary array.

        """))
    end
    nperthread, nrem = divrem(n_particles, nbatches)
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
julia> using CellListMap

julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(3) for i in 1:50_000 ];

julia> y = [ 250*rand(3) for i in 1:10_000 ];

julia> cl = CellList(x,y,box);

julia> aux = CellListMap.AuxThreaded(cl)
CellListMap.AuxThreaded{3, Float64}
 Auxiliary arrays for nthreads = 8

julia> UpdateCellList!(x,y,box,cl,aux)
CellList{3, Float64}
  100000 real particles.
  31190 cells with real particles.
  1134378 particles in computing box, including images.

```
"""
AuxThreaded(cl_pair::CellListPair) = 
    AuxThreadedPair(AuxThreaded(cl_pair.small_set), AuxThreaded(cl_pair.large_set))

"""
    CellList(
        x::AbstractVector{AbstractVector},
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0),
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
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
    nbatches::Tuple{Int,Int}=(0, 0),
    validate_coordinates::Union{Function,Nothing}=_validate_coordinates,
) where {UnitCellType,N,T}
    cl = CellList{N,T}(n_real_particles=length(x), number_of_cells=prod(box.nc))
    set_number_of_batches!(cl, nbatches; parallel)
    return UpdateCellList!(x, box, cl; parallel, validate_coordinates)
end

#=
    CellList(x::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T} 

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

=#
function CellList(
    x::AbstractMatrix, 
    box::Box{UnitCellType,N,T}; 
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0),
    validate_coordinates::Union{Function,Nothing}=_validate_coordinates,
) where {UnitCellType,N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return CellList(x_re, box; parallel, nbatches, validate_coordinates)
end

#=
    reset!(cl::CellList{N,T},box,n_real_particles) where{N,T}

# Extended help

Resets a cell list, by setting everything to zero, but retaining
the allocated `particles` and `projected_particles` vectors.

=#
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
    cl.n_real_particles = n_real_particles
    cl.n_particles = 0
    cl.number_of_cells = new_number_of_cells
    cl.n_cells_with_real_particles = 0
    cl.n_cells_with_particles = 0
    return cl
end

"""
    CellList(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0),
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
    ) where {UnitCellType,N,T} 

Function that will initialize a `CellListPair` structure from scratch, given two vectors
of particle coordinates and a `Box`, which contain the size of the system, cutoff, etc.

## Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> y = [ 250*rand(SVector{3,Float64}) for i in 1:10000 ];

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   10000 particles in the reference vector.
   961 cells with real particles of target vector.

```

"""
function CellList(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box{UnitCellType,N,T};
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0),
    validate_coordinates::Union{Function,Nothing}=_validate_coordinates,
) where {UnitCellType,N,T}
    xsmall, xlarge, swap = length(x) <= length(y) ? (x, y, false) : (y, x, true)
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    small_set = CellList(xsmall, box; parallel, validate_coordinates) 
    large_set = CellList(xlarge, box; parallel, validate_coordinates)
    cl_pair = CellListPair{N,T}(small_set, large_set, swap)
    cl_pair = set_number_of_batches!(cl_pair, nbatches; parallel)
    return cl_pair
end

#=
    CellList(x::AbstractMatrix, y::AbstractMatrix, box::Box{UnitCellType,N,T}; kargs...) where {UnitCellType,N,T} 

Reinterprets the matrices `x` and `y` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrices must be the dimension of the points (`2` or `3`).

=#
function CellList(
    x::AbstractMatrix, 
    y::AbstractMatrix, 
    box::Box{UnitCellType,N,T}; 
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0),
    validate_coordinates::Union{Function,Nothing}=_validate_coordinates,
) where {UnitCellType,N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return CellList(x_re, y_re, box; parallel, nbatches, validate_coordinates)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        box::Box,
        cl:CellList;
        parallel=true,
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
    ) 

Function that will update a previously allocated `CellList` structure, given new 
updated particle positions. This function will allocate new threaded auxiliary
arrays in parallel calculations. To preallocate these auxiliary arrays, use
the `UpdateCellList!(x,box,cl,aux)` method instead. 

The `validate_coordinates` function is called before the update of the cell list, and
should throw an error if the coordinates are invalid. By default, this function 
throws an error if some coordinates are missing or are NaN. Set to `nothing` to disable
this check, or provide a custom function.

## Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> cl = CellList(x,box);

julia> box = Box([260,260,260],10);

julia> x = [ 260*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> UpdateCellList!(x,box,cl); # update lists

```

"""
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    box::Box,
    cl::CellList;
    parallel::Bool=true,
    validate_coordinates=_validate_coordinates,
)
    cl = if parallel
        aux = AuxThreaded(cl)
        return UpdateCellList!(x, box, cl, aux; parallel, validate_coordinates)
    else
        return UpdateCellList!(x, box, cl, nothing; parallel, validate_coordinates)
    end
    cl = update_number_of_batches!(cl; parallel)
    return cl
end

#=
    function UpdateCellList!(
        x::AbstractMatrix,
        box::Box,
        cl::CellList{N,T};
        parallel::Bool=true,
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates,
    ) where {N,T}

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

=#
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T};
    kargs...
) where {N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return UpdateCellList!(x_re, box, cl; kargs...)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        box::Box,
        cl::CellList{N,T},
        aux::Union{Nothing,AuxThreaded{N,T}};
        parallel::Bool=true,
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
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
    parallel::Bool=true,
    validate_coordinates=_validate_coordinates,
) where {N,T}

    # validate coordinates 
    isnothing(validate_coordinates) || validate_coordinates(x)

    # Provide a better error message if the unit cell dimension does not match the dimension of the positions.
    if length(x) > 0 && (length(x[begin]) != size(box.input_unit_cell.matrix, 1))
        n1 = length(x[begin])
        n2 = size(box.input_unit_cell.matrix, 1)
        throw(DimensionMismatch("""\n 
            Positions have dimension $n1, but the unit cell has dimension $n2.

        """))
    end

    # Add particles to cell list
    _nbatches = nbatches(cl, :build)
    if !parallel || _nbatches == 1
        reset!(cl, box, length(x))
        add_particles!(x, box, 0, cl)
    else
        # Reset cell list
        reset!(cl, box, 0)
        # Update the aux.idxs ranges, for if the number of particles changed
        set_idxs!(aux.idxs, length(x), _nbatches)
        merge_lock = ReentrantLock()
        @sync for ibatch in eachindex(aux.idxs, aux.lists)
            @spawn let cl = $cl
                prange = aux.idxs[ibatch]
                aux.lists[ibatch] = reset!(aux.lists[ibatch], box, length(prange))
                xt = @view(x[prange])
                aux.lists[ibatch] = add_particles!(xt, box, prange[begin] - 1, aux.lists[ibatch])
                @lock merge_lock begin
                    merge_cell_lists!(cl, aux.lists[ibatch])
                end
            end
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

#=
    UpdateCellList!(
        x::AbstractMatrix,
        box::Box,
        cl::CellList{N,T},
        aux::Union{Nothing,AuxThreaded{N,T}};
        parallel::Bool=true,
        validate_coordinates=_validate_coordinates,
    ) where {N,T}

Reinterprets the matrix `x` as vectors of static vectors and calls the
equivalent function with the reinterprted input. The first dimension of the 
matrix must be the dimension of the points (`2` or `3`).

=#
function UpdateCellList!(
    x::AbstractMatrix,
    box::Box,
    cl::CellList{N,T},
    aux::Union{Nothing,AuxThreaded{N,T}};
    kargs...
) where {N,T}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return UpdateCellList!(x_re, box, cl, aux; kargs...)
end

#=
    add_particles!(x,box,ishift,cl::CellList{N,T}) where {N,T}

# Extended help

Add all particles in vector `x` to the cell list `cl`. `ishift` is the shift in particle
index, meaning that particle `i` of vector `x` corresponds to the particle with original
index `i+ishift`. The shift is used to construct cell lists from fractions of the original
set of particles in parallel list construction.  

=#
function add_particles!(x, box, ishift, cl::CellList{N,T}) where {N,T}
    for ip in eachindex(x)
        xp = x[ip]
        # This converts the coordinates to static arrays, if necessary
        p = SVector{N,T}(ntuple(i -> xp[i], N))
        p = box.rotation * wrap_to_first(p, box.input_unit_cell.matrix)
        add_particle_to_celllist!(ishift + ip, p, box, cl) # add real particle
        replicate_particle!(ishift + ip, p, box, cl) # add virtual particles to border cells
    end
    return cl
end

# define method for the ParticleWithIndex type
out_of_computing_box(p::ParticleWithIndex, box::Box) = out_of_computing_box(p.coordinates, box)

#=
    copydata!(cell1::Cell,cell2::Cell)

# Extended help

Copies the data from `cell2` to `cell1`, meaning that particles are
copied element-wise from `cell2` to `cell1`, with the `particles` array
of `cell1` being resized (increased) if necessary.

=#
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

#=
    append_particles!(cell1::Cell,cell2::Cell)

# Extended help

Add the particles of `cell2` to `cell1`, updating the cell data and, if necessary,
resizing (increasing) the `particles` array of `cell1`

=#
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

#=
    merge_cell_lists!(cl::CellList,aux::CellList)

# Extended help

Merges an auxiliary `aux` cell list to `cl`, and returns the modified `cl`. Used to
merge cell lists computed in parallel threads.

=#
function merge_cell_lists!(cl::CellList, aux::CellList)
    # One should never get here if the lists do not share the same # computing box
    if cl.number_of_cells != aux.number_of_cells
        throw(ArgumentError("""\n
            Cell lists must have the same number of cells to be merged.
            Got inconsistent number of cells: $(cl.number_of_cells) and $(aux.number_of_cells).

        """))
    end
    cl.n_particles += aux.n_particles
    cl.n_real_particles += aux.n_real_particles
    for icell in 1:aux.n_cells_with_particles
        aux_cell = aux.cells[icell]
        linear_index = aux_cell.linear_index
        cell_index = cl.cell_indices[linear_index]
        # If cell was yet not initialized in merge, push it to the list
        if cell_index == 0
            cl.n_cells_with_particles += 1
            cell_index = cl.n_cells_with_particles
            cl.cell_indices[linear_index] = cell_index
            if cell_index > length(cl.cells)
                push!(cl.cells, copy_cell(aux_cell))
            else
                cl.cells[cell_index] = copydata!(cl.cells[cell_index], aux_cell)
            end
            if aux_cell.contains_real
                cl.n_cells_with_real_particles += 1
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
                cl.n_cells_with_real_particles += 1
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

# Deal with corner cases where a real particle is found in the exact bundary of real box.
# This cannot happen because then running over the neighboring boxes can cause an 
# invalid access to an index of a cell.
function real_particle_border_case(cartesian_index::CartesianIndex{N}, box) where {N}
    cidxs = ntuple(i -> cartesian_index[i], N)
    for i in 1:N
        if cidxs[i] == box.lcell
            @set! cidxs[i] += 1
        end
        if cidxs[i] == box.nc[i] - box.lcell + 1
            @set! cidxs[i] -= 1
        end
    end
    return CartesianIndex{N}(cidxs)
end

#=
    add_particle_to_celllist!(
        ip,
        x::SVector{N,T},
        box,
        cl::CellList{N,T};
        real_particle::Bool=true
    ) where {N,T}

# Extended help

Adds one particle to the cell lists, updating all necessary arrays.

=#
function add_particle_to_celllist!(
    ip,
    x::SVector{N,T},
    box,
    cl::CellList{N,T};
    real_particle::Bool=true
) where {N,T}

    # Increase the counter of number of particles of this list
    cl.n_particles += 1
    # Cell of this particle
    cartesian_index = particle_cell(x, box)
    if real_particle
        cartesian_index = real_particle_border_case(cartesian_index, box)
    end

    # Linear index of the cell
    linear_index = cell_linear_index(box.nc, cartesian_index)

    # Check if this is the first particle of this cell, if it is,
    # initialize a new cell, or reset a previously allocated one
    cell_index = cl.cell_indices[linear_index]

    if cell_index == 0
        cl.n_cells_with_particles += 1
        cell_index = cl.n_cells_with_particles
        cl.cell_indices[linear_index] = cell_index
        if cell_index > length(cl.cells)
            particles_sizehint = cl.n_real_particles ÷ prod(box.nc)
            push!(cl.cells, Cell{N,T}(cartesian_index, box, sizehint=particles_sizehint))
        else
            cell = cl.cells[cell_index]
            @set! cell.linear_index = linear_index
            @set! cell.cartesian_index = cartesian_index
            @set! cell.center = cell_center(cell.cartesian_index, box)
            @set! cell.contains_real = false
            @set! cell.n_particles = 0
            cl.cells[cell_index] = cell
        end
    end

    # Increase particle counter for this cell
    cell = cl.cells[cell_index]
    @set! cell.n_particles += 1

    #
    # Cells with real particles are annotated to be run over
    #
    if real_particle && (!cell.contains_real)
        @set! cell.contains_real = true
        cl.n_cells_with_real_particles += 1
        if cl.n_cells_with_real_particles > length(cl.cell_indices_real)
            push!(cl.cell_indices_real, cell_index)
        else
            cl.cell_indices_real[cl.n_cells_with_real_particles] = cell_index
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
    cl.cells[cell_index] = cell

    return cl
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box,
        cl:CellListPair,
        parallel=true,
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
    )

Function that will update a previously allocated `CellListPair` structure, given 
updated particle positions, for example. This method will allocate new 
`aux` threaded auxiliary arrays. For a non-allocating version, see the 
`UpdateCellList!(x,y,box,cl,aux)` method.

The `validate_coordinates` function is called before the update of the cell list, and
should throw an error if the coordinates are invalid. By default, this function
throws an error if some coordinates are missing or are NaN. Set to `nothing` to disable
this check, or provide a custom function.

```jlcodetest
julia> using CellListMap, StaticArrays

julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> y = [ 250*rand(SVector{3,Float64}) for i in 1:10000 ];

julia> cl = CellList(x,y,box);

julia> UpdateCellList!(x,y,box,cl); # update lists
```

"""
function UpdateCellList!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector},
    box::Box,
    cl_pair::CellListPair;
    parallel::Bool=true,
    kargs...
)
    cl_pair = if parallel
        aux = AuxThreaded(cl_pair)
        UpdateCellList!(x, y, box, cl_pair, aux; kargs...)
    else
        UpdateCellList!(x, y, box, cl_pair, nothing; kargs...)
    end
    cl_pair = update_number_of_batches!(cl_pair; parallel)
    return cl_pair
end

#=
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

=#
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N},
    cl_pair::CellListPair;
    kargs...
) where {UnitCellType,N}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return UpdateCellList!(x_re, y_re, box, cl_pair; kargs...)
end

"""
    UpdateCellList!(
        x::AbstractVector{<:AbstractVector},
        y::AbstractVector{<:AbstractVector},
        box::Box,
        cl_pair::CellListPair,
        aux::Union{Nothing,AuxThreaded};
        parallel::Bool=true,
        validate_coordinates::Union{Function,Nothing}=_validate_coordinates
    )

This function will update the `cl_pair` structure that contains the cell lists
for disjoint sets of particles. It receives the preallocated `aux` structure to
avoid reallocating auxiliary arrays necessary for the threaded construct of the
lists. 

## Example

```julia-repl
julia> using CellListMap

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

julia> UpdateCellList!(x,y,box,cl,aux)
CellList{3, Float64}
  10000 real particles.
  7358 cells with real particles.
  12591 particles in computing box, including images.

```

The allocations in the above calls are a consequence of the thread spawning only:

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
    cl_pair::CellListPair{N,T},
    aux::Union{Nothing,AuxThreadedPair};
    parallel::Bool=true,
    validate_coordinates::Union{Nothing,Function}=_validate_coordinates,
) where {N,T}
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    xsmall, xlarge, swap = length(x) <= length(y) ? (x, y, false) : (y, x, true)
    small_set = UpdateCellList!(xsmall, box, cl_pair.small_set, aux.small_set; parallel, validate_coordinates)
    large_set = UpdateCellList!(xlarge, box, cl_pair.large_set, aux.large_set; parallel, validate_coordinates)
    return CellListPair{N,T}(small_set, large_set, swap)
end

#=
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

=#
function UpdateCellList!(
    x::AbstractMatrix,
    y::AbstractMatrix,
    box::Box{UnitCellType,N},
    cl_pair::CellListPair,
    aux::Union{Nothing,AuxThreadedPair};
    kargs...
) where {UnitCellType,N}
    size(x, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    size(y, 1) == N || throw(DimensionMismatch("First dimension of input matrix must be $N"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return UpdateCellList!(x_re, y_re, box, cl_pair, aux; kargs...)
end

#=
    particles_per_cell(cl)

Returns the average number of real particles per computing cell.

=#
particles_per_cell(cl::CellList) = cl.n_real_particles / cl.number_of_cells
particles_per_cell(cl::CellListPair) = particles_per_cell(cl.small_set) + particle_cell(cl.large_set)

@testitem "celllists - validate coordinates" begin
    using CellListMap
    using StaticArrays
    x = rand(SVector{3,Float64}, 100)
    x[50] = SVector(1.0, NaN, 1.0)
    box = Box([1.0, 1.0, 1.0], 0.1)
    y = rand(SVector{3,Float64}, 100)
    @test_throws ArgumentError CellList(x, box)
    cl = CellList(y, box)
    @test_throws ArgumentError UpdateCellList!(x, box, cl)
    @test_throws ArgumentError CellList(x, y, box)
    @test_throws ArgumentError CellList(y, x, box)
    cl = CellList(y, y, box)
    @test_throws ArgumentError UpdateCellList!(x, y, box, cl)
    @test_throws ArgumentError UpdateCellList!(y, x, box, cl)
    x = rand(3, 100)
    x[2, 50] = NaN
    box = Box([1.0, 1.0, 1.0], 0.1)
    y = rand(3, 100)
    @test_throws ArgumentError CellList(x, box)
    cl = CellList(y, box)
    @test_throws ArgumentError UpdateCellList!(x, box, cl)
    @test_throws ArgumentError CellList(x, y, box)
    @test_throws ArgumentError CellList(y, x, box)
    cl = CellList(y, y, box)
    @test_throws ArgumentError UpdateCellList!(x, y, box, cl)
    @test_throws ArgumentError UpdateCellList!(y, x, box, cl)
end
