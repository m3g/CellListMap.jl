#=

$(TYPEDEF)

Structure that containst the system information for neighborlist computations. All fields are internal.

## Extended help

$(TYPEDFIELDS)

=#
mutable struct InPlaceNeighborList{P}
    sys::P
    show_progress::Bool
end

# getters
get_unitcell(nb::InPlaceNeighborList) = nb.sys.unitcell 
get_cutoff(nb::InPlaceNeighborList) = nb.sys.cutoff

"""
    InPlaceNeighborList(;
        xpositions::AbstractVecOrMat,
        #or
        positions::AbstractVecOrMat,
        ypositions::Union{AbstractVecOrMat,Nothing}=nothing,
        cutoff::T,
        unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
        parallel::Bool=true,
        show_progress::Bool=false,
    ) where {T}

Function that initializes the `InPlaceNeighborList` structure, to be used for in-place
computation of neighbor lists.

- `xpositions` (or `positions` as a shortcut for single-set computations) provides the
  coordinates of the particles.
- If `ypositions` is also provided, the neighbor list between the two sets is computed.
- If `unitcell` is provided, periodic boundary conditions will be used. The `unitcell` can
  be a vector of Orthorhombic box sides, or an actual unitcell matrix for general cells.
- If `unitcell` is not provided (value `nothing`), no periodic boundary conditions will
  be considered.

## Examples

Here the neighborlist structure is constructed for the first time, and used
to compute the neighbor lists with the mutating `neighborlist!` function:

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 10^4);

julia> system = InPlaceNeighborList(xpositions=x, cutoff=0.1, unitcell=[1,1,1])
InPlaceNeighborList with types:
CellList{3, Float64}
Box{OrthorhombicCell, 3, Float64, Float64, 9}
Current list buffer length: 0

julia> neighborlist!(system)
210034-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 357, 0.09922225615002134)
 (1, 488, 0.043487074695938925)
 (1, 2209, 0.017779967072139684)
 ⋮
 (9596, 1653, 0.0897570322108541)
 (9596, 7927, 0.0898266280344037)
```

The coordinates of the system, its unitcell, or the cutoff can be changed with
the `update!` function. If the number of pairs of the list does not change
significantly, the new calculation is minimally allocating, or non-allocating
at all, in particular if the computation is run without parallelization:

!!! note
    The order of the pairs in the output of `neighborlist!` is not guaranteed,
    and may change, in particular, in parallel runs.

If the structure is used repeatedly for similar systems, the allocations will
vanish, except for minor allocations used in the threading computation (if a
non-parallel computation is executed, the allocations will vanish completely):

"""
function InPlaceNeighborList(; show_progress::Bool = false, cutoff, unitcell = nothing, kargs...)
    T = _promote_types(unitcell, cutoff)
    sys = ParticleSystem(;
        cutoff=cutoff,
        unitcell=unitcell,
        output=NeighborList{T}(0, Vector{Tuple{Int, Int, T}}[]),
        output_name=:nb,
        kargs...
    )
    return InPlaceNeighborList(sys, show_progress)
end

"""
    update!(
        system::InPlaceNeighborList;
        xpositions = nothing,
        #or
        positions = nothing,
        ypositions = nothing,
        cutoff = nothing,
        unitcell = nothing,
    )

Updates a `InPlaceNeighborList` system. Only the keyword arguments that are provided
(i.e. not `nothing`) are updated. `positions` is accepted as an alias for `xpositions`.

## Examples

### For self-pairs computations

```julia-repl
julia> x = rand(SVector{3,Float64}, 10^3);

julia> system = InPlaceNeighborList(xpositions=x, cutoff=0.1)
InPlaceNeighborList with types:
CellList{3, Float64}
Box{NonPeriodicCell, 3, Float64, Float64, 9}
Current list buffer length: 0

julia> neighborlist!(system);

julia> new_x = rand(SVector{3,Float64}, 10^3);

julia> update!(system; xpositions=new_x, cutoff=0.05)
InPlaceNeighborList with types:
CellList{3, Float64}
Box{NonPeriodicCell, 3, Float64, Float64, 9}
Current list buffer length: 1826

julia> neighborlist!(system)
224-element Vector{Tuple{Int64, Int64, Float64}}:
 (25, 486, 0.03897345036790646)
 ⋮
 (723, 533, 0.04795768478723409)
 (868, 920, 0.042087156715720137)
```

"""
function update!(system::InPlaceNeighborList; kargs...) 
    update!(system.sys; kargs...)
    return system
end

function Base.show(io::IO, ::MIME"text/plain", system::InPlaceNeighborList)
    _print(io, "InPlaceNeighborList with types: \n")
    _print(io, typeof(system.sys._cell_list), "\n")
    _print(io, typeof(system.sys._box), "\n")
    return _print(io, "Current list buffer length: $(length(system.sys.nb.list))")
end

"""
    neighborlist(system::InPlaceNeighborList)

Computes the neighbor list in-place, given a `InPlaceNeighborList` system.

## Example

In the following example, we compute the neighbor list of a set of random
particles, and then we change the coordinates and recompute the neighbor list
without reallocations (or with minimal allocations if run in parallel):

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 10^4);

julia> system = InPlaceNeighborList(xpositions=x, cutoff=0.1, unitcell=[1,1,1], parallel=false);

julia> neighborlist!(system)
210034-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 2669, 0.04444346517920411)
 (1, 8475, 0.02554075837438248)
 ⋮
 (9463, 5955, 0.08698158178214915)
 (9463, 2308, 0.09482635540291776)

julia> x .= rand(SVector{3,Float64}, 10^4); # change coordinates

julia> @time neighborlist!(system; parallel=false)
  0.007978 seconds
209418-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1253, 0.09420839394144173)
 (1, 2048, 0.01448095691254145)
 ⋮
 (9728, 6367, 0.08204145034963985)
 (9728, 2594, 0.06536277710826768)
```

"""
function neighborlist!(system::InPlaceNeighborList)
    (; sys) = system
    sys.output = _reset_all_output!(sys.output, sys._output_threaded; reset = true)
    UpdateParticleSystem!(sys)
    _sizehint_neighbor_lists!(sys)
    sys.output = _pairwise!(
        push_pair!,
        sys.output, sys._box, sys._cell_list;
        output_threaded = sys._output_threaded,
        parallel = sys.parallel,
        show_progress = system.show_progress,
    )
    resize!(sys.output.list, sys.output.n)
    return sys.output.list
end

"""
    neighborlist(;
        xpositions,
        #or
        positions,
        ypositions=nothing,
        cutoff,
        unitcell=nothing,
        parallel=true,
        show_progress=false,
        nbatches=(0,0),
    )

Computes the list of pairs of particles closer to each other than `cutoff`.

- `xpositions` (or `positions` as a shortcut for single-set computations) provides the
  coordinates of the particles.
- If `ypositions` is also provided, returns the cross-pair neighbor list between the two sets.
- If `unitcell` is provided (as a vector of sides or a unit cell matrix), periodic boundary
  conditions are used.

!!! note
    The order of the pairs in the output of `neighborlist` is not guaranteed,
    and may change, in particular, in parallel runs.

## Example

Compute the neighborlist within a set of Argon atoms, non-periodic:

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file));

julia> neighborlist(xpositions=x, cutoff=8.0, parallel=false)
857-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 20, 3.163779526466901)
 (1, 61, 4.08865164675529)
 (1, 67, 5.939772435456664)
 ⋮
 (78, 88, 7.0061163797598445)
 (88, 54, 7.933654063435482)
```

And now, considering the system periodic:

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file));

julia> neighborlist(xpositions=x, cutoff=8.0, unitcell=[21.0, 21.0, 21.0], parallel=false)
1143-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 7, 3.3638756414119397)
 (1, 20, 3.163779526466901)
 (1, 47, 6.243868666689442)
 ⋮
 (68, 38, 7.409287768713663)
 (68, 90, 7.8758006026725464)
```

For cross-pair neighbor lists between two sets:

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file, "index <= 50"));

julia> y = coor(read_pdb(CellListMap.argon_pdb_file, "index > 50"));

julia> CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=8.0, parallel=false)
439-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 11, 4.08865164675529)
 (1, 17, 5.939772435456664)
 (1, 30, 2.4572288423012236)
 ⋮
 (46, 48, 4.9269093987894745)
 (46, 1, 7.99947286297016)
```

"""
neighborlist(; kargs...) = neighborlist!(InPlaceNeighborList(; kargs...))
