#
# This file contains the structures and functions for the construction of the
# cell lists (`CellList` and `CellListPair`).
#
# The cell lists are built on a grid defined in fractional coordinates of the unit
# cell (see `Box`). Each particle is stored once, wrapped into the unit cell, together
# with the lattice translation (`shift`) applied by the wrapping. The particles are
# sorted by cell (counting sort), such that the particles of each cell are contiguous
# in memory. Periodic images are not replicated: they are reached in the pairwise
# computation by wrapping the index of the neighboring cells.
#

#=
Coordinates of a particle (wrapped into the unit cell) and its index in the
original array of coordinates.

=#
struct ParticleWithIndex{N, T}
    index::Int
    coordinates::SVector{N, T}
end

#=
Structure to define the number of batches used in the parallel splitting of the calculations
of the cell list construction and of the `pairwise` computation. It is initialized with
a standard heuristic that returns at most the number of threads, but may return a smaller
number if the system is small.

=#
struct NumberOfBatches
    build_cell_lists::Tuple{Bool, Int}
    map_computation::Tuple{Bool, Int}
end
NumberOfBatches(auto::Tuple{Bool, Bool}, n::Tuple{Int, Int}) =
    NumberOfBatches((first(auto), first(n)), (last(auto), last(n)))
Base.zero(::Type{NumberOfBatches}) = NumberOfBatches((true, 0), (true, 0))
function Base.show(io::IO, ::MIME"text/plain", nbatches::NumberOfBatches)
    _println(io, "  Number of batches for cell list construction: $(last(nbatches.build_cell_lists)) (auto-update: $(first(nbatches.build_cell_lists)))")
    return _print(io, "  Number of batches for function mapping: $(last(nbatches.map_computation)) (auto-update: $(first(nbatches.map_computation)))")
end

#=
Structure that contains the cell lists of one set of particles.

=#
Base.@kwdef mutable struct CellList{N, T}
    " Number of particles. "
    n_real_particles::Int = 0
    " Number of computing cells. "
    number_of_cells::Int = 0
    " Particles (wrapped coordinates and original index), sorted by cell. "
    particles::Vector{ParticleWithIndex{N, T}} = ParticleWithIndex{N, T}[]
    " Lattice translation applied to each particle by the wrapping (sorted as `particles`): `particles[k].coordinates = x[index] + unitcell * shifts[k]`. "
    shifts::Vector{SVector{N, Int32}} = SVector{N, Int32}[]
    " Particles of cell `c` are `particles[cell_start[c]:cell_start[c+1]-1]`. "
    cell_start::Vector{Int} = [1]
    " Linear indices of the cells that contain particles. "
    nonempty_cells::Vector{Int} = Int[]
    " Scratch buffers used in the construction of the lists. "
    cell_of::Vector{Int} = Int[]
    wrapped::Vector{SVector{N, T}} = SVector{N, T}[]
    wrap_shift::Vector{SVector{N, Int32}} = SVector{N, Int32}[]
    write_pos::Vector{Int} = Int[]
    " Offsets of neighboring cells: one of each {v,-v} (self-computations) and all of them (cross-computations). "
    offsets_half::Vector{SVector{N, Int}} = SVector{N, Int}[]
    offsets_full::Vector{SVector{N, Int}} = SVector{N, Int}[]
    " Key of the box for which the offsets were computed. "
    offsets_key::UInt = UInt(0)
    " Number of batches for the parallel calculations. "
    nbatches::NumberOfBatches = zero(NumberOfBatches)
end

function Base.show(io::IO, ::MIME"text/plain", cl::CellList)
    _println(io, typeof(cl))
    _println(io, "  $(cl.n_real_particles) particles.")
    return _print(io, "  $(length(cl.nonempty_cells)) cells with particles.")
end

#=
Structure that will contains the cell lists of two independent sets of
particles for cross-computation of interactions

=#
struct CellListPair{N, T}
    ref_list::CellList{N, T}  # Reference list (iterated in outer loop)
    target_list::CellList{N, T}  # Target list (searched for neighbors)
end

function Base.show(io::IO, ::MIME"text/plain", cl::CellListPair)
    _print(io, typeof(cl), "\n")
    _print(io, "   $(length(cl.ref_list.nonempty_cells)) cells with particles in reference set.\n")
    return _print(io, "   $(length(cl.target_list.nonempty_cells)) cells with particles in target set.")
end

#=
    update_number_of_batches!(cl, nbatches::NumberOfBatches; parallel=true)

Set the default number of batches for the construction of the cell lists,
and mapping computations. This is of course heuristic, and may not be the best choice for
every problem. See the parameter `nbatches` of the construction of the cell lists for
tuning this.

=#
function update_number_of_batches!(cl::CellList{N, T}, _nbatches = cl.nbatches; parallel = true) where {N, T}
    auto = (first(_nbatches.build_cell_lists), first(_nbatches.map_computation))
    n1 = last(_nbatches.build_cell_lists)
    n2 = last(_nbatches.map_computation)
    if !parallel
        if !all(auto) && (n1, n2) != (1, 1)
            @warn begin
                """\n
                WARNING: nbatches set to ($n1,$n2), but parallel is set to false, implying nbatches == (1, 1)

                """
            end _file = nothing _line = nothing
        end
        nbatches = NumberOfBatches((false, false), (1, 1))
    else # Heuristic choices
        if first(auto)
            n1 = _nbatches_build_cell_lists(cl.n_real_particles)
        end
        if last(auto)
            n2 = _nbatches_map_computation(cl.n_real_particles)
        end
        nbatches = NumberOfBatches(auto, (n1, n2))
    end
    cl.nbatches = nbatches
    return cl
end

# Heuristic choices for the number of batches, for an atomic system
_nbatches_default(n::Int) = max(1, min(n, min(floor(Int, 2^(log10(n) + 1)), nthreads())))
_nbatches_build_cell_lists(n::Int) = _nbatches_default(n)
_nbatches_map_computation(n::Int) = _nbatches_default(n)

function update_number_of_batches!(
        cl::CellListPair{N, T},
        _nbatches::Union{NumberOfBatches, Nothing} = nothing;
        parallel = true
    ) where {N, T}
    if isnothing(_nbatches)
        # Default call (e.g. when particle count changed): reconstruct auto flags from stored
        # state. build_cell_lists.auto is correctly preserved across updates. map_computation.auto
        # is always stored as false in individual lists (to prevent per-list recalculation based
        # on individual particle counts), so we use build_auto as the pair-level map auto flag.
        build_auto = first(cl.ref_list.nbatches.build_cell_lists)
        _nbatches = NumberOfBatches(
            (build_auto, build_auto),
            (last(cl.ref_list.nbatches.build_cell_lists), last(cl.ref_list.nbatches.map_computation))
        )
    end
    auto = (first(_nbatches.build_cell_lists), first(_nbatches.map_computation))
    n_build_ref = last(_nbatches.build_cell_lists)
    n_build_target = last(_nbatches.build_cell_lists)
    n_map = last(_nbatches.map_computation)

    if !parallel
        nbatches_ref = NumberOfBatches((false, false), (1, 1))
        nbatches_target = NumberOfBatches((false, false), (1, 1))
    else
        # Build: independent for each set
        if first(auto)
            n_build_ref = _nbatches_build_cell_lists(cl.ref_list.n_real_particles)
            n_build_target = _nbatches_build_cell_lists(cl.target_list.n_real_particles)
        end
        # Map: based on product of particle counts (cross-interaction complexity)
        if last(auto)
            n_total = cl.ref_list.n_real_particles * cl.target_list.n_real_particles
            n_map = _nbatches_map_computation(n_total)
        end
        nbatches_ref = NumberOfBatches((first(auto), false), (n_build_ref, n_map))
        nbatches_target = NumberOfBatches((first(auto), false), (n_build_target, n_map))
    end

    ref_list = update_number_of_batches!(cl.ref_list, nbatches_ref; parallel)
    target_list = update_number_of_batches!(cl.target_list, nbatches_target; parallel)
    return CellListPair{N, T}(ref_list, target_list)
end

#
# Functions for initialization of the batches, called from the API functions. Receives a
# tuple with the number of batches for the construction of the cell lists and the mapping.
# If the number of batches are smaller than 1, the function will set the number of batches
# in automatic mode.
#
function set_number_of_batches!(cl::Union{CellList, CellListPair}, _nbatches::Tuple{Int, Int}; parallel = true)
    auto = _nbatches .<= 0
    nbatches = NumberOfBatches((first(auto), first(_nbatches)), (last(auto), last(_nbatches)))
    return update_number_of_batches!(cl, nbatches; parallel)
end

#=
    nbatches(cl, [:map|:build])

Returns the number of batches for parallel processing that will be used. Returns a tuple, where
the first element is the number of batches for the construction of the cell lists, and the second
element is the number of batches for the pairwise mapping.

=#
function nbatches(cl::CellList, s::Symbol)
    if s == :map_computation || s == :map
        return last(cl.nbatches.map_computation)
    elseif s == :build_cell_lists || s == :build
        return last(cl.nbatches.build_cell_lists)
    end
end
nbatches(cl::CellList) = (nbatches(cl, :build), nbatches(cl, :map))
nbatches(cl::CellListPair) = nbatches(cl.ref_list)
nbatches(cl::CellListPair, s::Symbol) = nbatches(cl.ref_list, s)

#=
Auxiliary structures for the parallel construction of the cell lists. The construction
does not require auxiliary data in the current implementation (all buffers are part of the
`CellList` structure), but the types are kept for compatibility of the internal interface.

=#
struct AuxThreaded{N, T} end
struct AuxThreadedPair{N, T}
    ref_list::AuxThreaded{N, T}
    target_list::AuxThreaded{N, T}
end
AuxThreaded(::CellList{N, T}) where {N, T} = AuxThreaded{N, T}()
AuxThreaded(::CellListPair{N, T}) where {N, T} = AuxThreadedPair{N, T}(AuxThreaded{N, T}(), AuxThreaded{N, T}())
function Base.show(io::IO, ::MIME"text/plain", aux::Union{AuxThreaded, AuxThreadedPair})
    return _print(io, typeof(aux))
end
_create_aux(::Box, cl::Union{CellList, CellListPair}) = AuxThreaded(cl)

#=
    set_idxs!(idxs, n_particles, nbatches)

# Extended help

Sets the indexes of the particles that will be considered for each batch in parallel runs.

=#
function set_idxs!(idxs, n_particles, nbatches)
    if length(idxs) != nbatches
        throw(
            ArgumentError(
                """\n
                    Modifying `nbatches` requires an explicit update of the AuxThreaded auxiliary array.

                """
            )
        )
    end
    nperthread, nrem = divrem(n_particles, nbatches)
    first = 1
    for ibatch in eachindex(idxs)
        nx = nperthread
        if ibatch <= nrem
            nx += 1
        end
        idxs[ibatch] = first:((first - 1) + nx)
        first += nx
    end
    return nothing
end

#=
    CellList(
        x::ParticleSystemPositions,
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0),
        validate_coordinates::Function=_validate_coordinates
    ) where {UnitCellType,N,T}

Function that will initialize a `CellList` structure from scratch, given a vector
or particle coordinates and a `Box`, which contain the size of the system, cutoff, etc.

=#
function CellList(
        x::ParticleSystemPositions,
        box::Box{UnitCellType, N, T};
        parallel::Bool = true,
        nbatches::Tuple{Int, Int} = (0, 0),
        validate_coordinates::F = _validate_coordinates,
    ) where {UnitCellType, N, T, F <: Function}
    cl = CellList{N, T}(n_real_particles = length(x), number_of_cells = prod(box.nc))
    set_number_of_batches!(cl, nbatches; parallel)
    return UpdateCellList!(x, box, cl, nothing; parallel, validate_coordinates)
end

#=
    CellList(
        x::ParticleSystemPositions,
        y::ParticleSystemPositions,
        box::Box{UnitCellType,N,T};
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0,0),
        validate_coordinates::Function=_validate_coordinates
    ) where {UnitCellType,N,T}

Function that will initialize a `CellListPair` structure from scratch, given two vectors
of particle coordinates and a `Box`, which contain the size of the system, cutoff, etc.

=#
function CellList(
        x::ParticleSystemPositions,
        y::ParticleSystemPositions,
        box::Box{UnitCellType, N, T};
        parallel::Bool = true,
        nbatches::Tuple{Int, Int} = (0, 0),
        validate_coordinates::F = _validate_coordinates,
    ) where {UnitCellType, N, T, F <: Function}
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    ref_list = CellList(x, box; parallel, validate_coordinates)
    target_list = CellList(y, box; parallel, validate_coordinates)
    cl_pair = CellListPair{N, T}(ref_list, target_list)
    cl_pair = set_number_of_batches!(cl_pair, nbatches; parallel)
    return cl_pair
end

#=
    UpdateCellList!(
        x::ParticleSystemPositions,
        box::Box,
        cl::CellList,
        [aux::Union{Nothing,AuxThreaded}];
        parallel=true,
        validate_coordinates::Function=_validate_coordinates,
    )

Function that will update a previously allocated `CellList` structure, given new
updated particle positions. The method without the `aux` argument only updates
the cell lists if the positions were flagged as updated, and updates the number of
batches for parallel computations (which depends on the number of particles).

The `validate_coordinates` function is called before the update of the cell list, and
should throw an error if the coordinates are invalid. By default, this function
throws an error if some coordinates are missing or are NaN. Set to `nothing` to disable
this check, or provide a custom function.

=#
function UpdateCellList!(
        x::ParticleSystemPositions,
        box::Box,
        cl::CellList;
        parallel::Bool = true,
        validate_coordinates::F = _validate_coordinates,
    ) where {F <: Function}
    if x.updated[]
        cl = UpdateCellList!(x, box, cl, nothing; parallel, validate_coordinates)
        cl = update_number_of_batches!(cl; parallel)
    end
    return cl
end

function UpdateCellList!(
        x::ParticleSystemPositions,
        box::Box{UnitCellType, N},
        cl::CellList{N, T},
        ::Union{Nothing, AuxThreaded};
        parallel::Bool = true,
        validate_coordinates::F = _validate_coordinates,
    ) where {UnitCellType, N, T, F <: Function}
    isnothing(validate_coordinates) || validate_coordinates(x)
    # Provide a better error message if the unit cell dimension does not match the dimension of the positions.
    if length(x) > 0 && (length(x[begin]) != size(box.input_unit_cell.matrix, 1))
        n1 = length(x[begin])
        n2 = size(box.input_unit_cell.matrix, 1)
        throw(
            DimensionMismatch(
                """\n
                    Positions have dimension $n1, but the unit cell has dimension $n2.

                """
            )
        )
    end
    cl.n_real_particles = length(x)
    cl.number_of_cells = prod(box.nc)
    _update_offsets!(cl, box)
    _build!(cl, x, box; parallel)
    # set updated to false, to indicate that the cell list was constructed
    # for these coordinates.
    x.updated[] = false
    return cl
end

#
# Update the offsets of neighboring cells, if the box changed since they were computed.
#
function _update_offsets!(cl::CellList, box::Box)
    key = hash((box.input_unit_cell.matrix, box.cutoff, box.nc, box.n_search))
    if key != cl.offsets_key
        _cell_offsets!(cl.offsets_half, cl.offsets_full, box)
        cl.offsets_key = key
    end
    return cl
end

#
# Fractional-coordinate cell of a particle, its coordinates wrapped into the unit cell,
# and the lattice translation applied by the wrapping, such that
# wrapped == x + unitcell * shift.
#
@inline function _locate(x::SVector{N, T}, box::Box{<:PeriodicCellType, N, T}) where {N, T}
    (; nc) = box
    m = box.input_unit_cell.matrix
    f = box.inv_unit_cell * ((x - box.origin) ./ oneunit(T))
    fl = floor.(f)
    fw = f - fl
    # Boundary coordinates belong to the lower boundary
    at_upper = fw .== one(eltype(f))
    fw = ifelse.(at_upper, zero(eltype(f)), fw)
    fl = ifelse.(at_upper, fl .+ one(eltype(f)), fl)
    c = min.(floor.(Int, fw .* nc), nc .- 1)
    shift = SVector{N, Int32}(Int32.(-fl))
    wrapped = x + m * shift
    return _cell_linear_index(c, nc), wrapped, shift
end

# Non-periodic systems: particles are inside the box, no wrapping
@inline function _locate(x::SVector{N, T}, box::Box{NonPeriodicCell, N, T}) where {N, T}
    (; nc) = box
    f = box.inv_unit_cell * ((x - box.origin) ./ oneunit(T))
    c = clamp.(floor.(Int, f .* nc), 0, nc .- 1)
    return _cell_linear_index(c, nc), x, zero(SVector{N, Int32})
end

# Linear index (1-based) of a cell given its 0-based cartesian indices
@inline function _cell_linear_index(c::SVector{N, Int}, nc::SVector{N, Int}) where {N}
    li = c[N]
    for d in (N - 1):-1:1
        li = li * nc[d] + c[d]
    end
    return li + 1
end

# 0-based cartesian indices of a cell given its (1-based) linear index
@inline function _cell_cartesian_index(li::Int, nc::SVector{N, Int}) where {N}
    r = li - 1
    return SVector{N, Int}(
        ntuple(Val(N)) do d
            q, m = divrem(r, nc[d])
            r = q
            m
        end
    )
end

#=
    _build!(cl, x, box; parallel)

# Extended help

Builds the cell lists: (1) computes the cell, wrapped coordinates, and wrap shift of
each particle (in parallel, if requested); (2) counts the particles in each cell;
(3) scatters the particles into the cell-sorted arrays.

=#
function _build!(cl::CellList{N, T}, x::ParticleSystemPositions, box::Box; parallel::Bool) where {N, T}
    n = length(x)
    ncell = cl.number_of_cells
    resize!(cl.cell_of, n)
    resize!(cl.wrapped, n)
    resize!(cl.wrap_shift, n)
    resize!(cl.particles, n)
    resize!(cl.shifts, n)
    resize!(cl.cell_start, ncell + 1)
    resize!(cl.write_pos, ncell)

    # 1. Locate particles
    _nbatches = nbatches(cl, :build)
    if !parallel || _nbatches == 1 || n == 0
        _locate_range!(cl, x, box, 1:n)
    else
        n_chunks = _n_workqueue_chunks(_nbatches, n)
        chunks = collect(index_chunks(1:n; n = n_chunks, split = Consecutive()))
        next_chunk = Atomic{Int}(1)
        tasks = Vector{Task}(undef, _nbatches)
        for ibatch in 1:_nbatches
            tasks[ibatch] = @spawn begin
                while true
                    ic = atomic_add!($next_chunk, 1)
                    ic > length($chunks) && break
                    _locate_range!($cl, $x, $box, $chunks[ic])
                end
            end
        end
        foreach(wait, tasks)
    end

    # 2. Count particles per cell, and compute the first index of each cell
    cell_start = cl.cell_start
    fill!(cell_start, 0)
    @inbounds for i in 1:n
        cell_start[cl.cell_of[i] + 1] += 1
    end
    empty!(cl.nonempty_cells)
    cell_start[1] = 1
    @inbounds for c in 1:ncell
        cell_start[c + 1] > 0 && push!(cl.nonempty_cells, c)
        cell_start[c + 1] += cell_start[c]
    end

    # 3. Scatter particles into the cell-sorted arrays
    copyto!(cl.write_pos, 1, cell_start, 1, ncell)
    @inbounds for i in 1:n
        c = cl.cell_of[i]
        p = cl.write_pos[c]
        cl.write_pos[c] = p + 1
        cl.particles[p] = ParticleWithIndex(i, cl.wrapped[i])
        cl.shifts[p] = cl.wrap_shift[i]
    end
    return cl
end

function _locate_range!(cl::CellList{N, T}, x::ParticleSystemPositions, box::Box, range) where {N, T}
    @inbounds for i in range
        xi = x[i]
        # This converts the coordinates to static arrays of the type of the cell list, if necessary
        p = SVector{N, T}(ntuple(k -> xi[k], Val(N)))
        cl.cell_of[i], cl.wrapped[i], cl.wrap_shift[i] = _locate(p, box)
    end
    return cl
end

#=
    UpdateCellList!(
        x::ParticleSystemPositions,
        y::ParticleSystemPositions,
        box::Box,
        cl_pair::CellListPair,
        [aux::Union{Nothing,AuxThreadedPair}];
        parallel=true,
        validate_coordinates::Function=_validate_coordinates,
    )

Function that will update a previously allocated `CellListPair` structure, given
updated particle positions. The method without `aux` only updates the lists of the
sets for which the positions were flagged as updated, and updates the number of batches
for parallel computations.

=#
function UpdateCellList!(
        x::ParticleSystemPositions,
        y::ParticleSystemPositions,
        box::Box,
        cl_pair::CellListPair;
        parallel::Bool = true,
        validate_coordinates::F = _validate_coordinates,
    ) where {F <: Function}
    x_updated, y_updated = x.updated[], y.updated[]
    ref_list = x_updated ? UpdateCellList!(x, box, cl_pair.ref_list, nothing; parallel, validate_coordinates) : cl_pair.ref_list
    target_list = y_updated ? UpdateCellList!(y, box, cl_pair.target_list, nothing; parallel, validate_coordinates) : cl_pair.target_list
    cl_pair = CellListPair(ref_list, target_list)
    if x_updated || y_updated
        cl_pair = update_number_of_batches!(cl_pair; parallel)
    end
    return cl_pair
end

function UpdateCellList!(
        x::ParticleSystemPositions,
        y::ParticleSystemPositions,
        box::Box,
        cl_pair::CellListPair{N, T},
        aux::Union{Nothing, AuxThreadedPair};
        parallel::Bool = true,
        validate_coordinates::F = _validate_coordinates,
    ) where {N, T, F <: Function}
    ref_aux = isnothing(aux) ? nothing : aux.ref_list
    target_aux = isnothing(aux) ? nothing : aux.target_list
    ref_list = UpdateCellList!(x, box, cl_pair.ref_list, ref_aux; parallel, validate_coordinates)
    target_list = UpdateCellList!(y, box, cl_pair.target_list, target_aux; parallel, validate_coordinates)
    return CellListPair{N, T}(ref_list, target_list)
end
