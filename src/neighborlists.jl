export InPlaceNeighborList
export update!
export neighborlist, neighborlist!

#
# Wrapper of the list of neighbors, that allows in-place updating of the lists
#
mutable struct NeighborList{T}
    n::Int
    list::Vector{Tuple{Int,Int,T}}
end

import Base: push!, empty!, resize!, copy
empty!(x::NeighborList) = x.n = 0
function push!(x::NeighborList, pair)
    x.n += 1
    if x.n > length(x.list)
        push!(x.list, pair)
    else
        x.list[x.n] = pair
    end
    return x
end
function resize!(x::NeighborList, n::Int)
    x.n = n
    resize!(x.list, n)
    return x
end
copy(x::NeighborList{T}) where {T} = NeighborList{T}(x.n, copy(x.list))

@testitem "NeighborList operations" begin
    using CellListMap
    nb = CellListMap.NeighborList(0, Tuple{Int,Int,Float64}[])
    @test length(nb.list) == 0
    push!(nb, (0, 0, 0.0))
    @test (nb.n, length(nb.list)) == (1, 1)
    empty!(nb)
    @test (nb.n, length(nb.list)) == (0, 1)
    resize!(nb, 5)
    @test (nb.n, length(nb.list), nb.n) == (5, 5, 5)
    nb.list = [(i, i, Float64(i)) for i in 1:5]
    nb2 = copy(nb)
    @test (nb.n, nb.list) == (nb2.n, nb2.list)
end

# Function adds pair to the list
function push_pair!(i, j, d2, list::NeighborList)
    d = sqrt(d2)
    push!(list, (i, j, d))
    return list
end

# We have to define our own reduce function here (for the parallel version)
# this reduction can be dum assynchronously on a preallocated array
function reduce_lists(list::NeighborList{T}, list_threaded::Vector{<:NeighborList{T}}) where {T}
    ranges = cumsum(nb.n for nb in list_threaded)
    npairs = ranges[end]
    # need to resize here for the case where length(list) < npairs
    list = resize!(list, npairs)
    @sync for it in eachindex(list_threaded)
        lt = list_threaded[it]
        range = ranges[it]-lt.n+1:ranges[it]
        @spawn list.list[range] .= @view(lt.list[1:lt.n])
    end
    return list
end

@testitem "Neighborlist push/reduce" begin
    using CellListMap
    nb1 = CellListMap.NeighborList(2, [(0, 0, 0.0), (1, 1, 1.0)])
    CellListMap.push_pair!(3, 3, 9.0, nb1)
    @test (nb1.n, nb1.list[3]) == (3, (3, 3, 3.0))
    nb2 = [copy(nb1), CellListMap.NeighborList(1, [(4, 4, 4.0)])]
    CellListMap.reduce_lists(nb1, nb2)
    @test nb1.n == 4
    @test nb1.list == [(0, 0, 0.0), (1, 1, 1.0), (3, 3, 3.0), (4, 4, 4.0)]
end

#=

$(TYPEDEF)

Structure that containst the system information for neighborlist computations. All fields are internal.

## Extended help

$(TYPEDFIELDS)

=#
mutable struct InPlaceNeighborList{B,C,A,NB<:NeighborList}
    box::B
    cl::C
    aux::A
    nb::NB
    nb_threaded::Vector{NB}
    parallel::Bool
    show_progress::Bool
end

"""
    InPlaceNeighborList(;
        x::AbstractVecOrMat,
        y::Union{AbstractVecOrMat,Nothing}=nothing,
        cutoff::T,
        unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
        parallel::Bool=true,
        show_progress::Bool=false,
    ) where {T}

Function that initializes the `InPlaceNeighborList` structure, to be used for in-place
computation of neighbor lists.

- If only `x` is provided, the neighbor list of the set is computed. 
- If `x` and `y` are provided, the neighbor list between the sets is computed.
- If `unitcell` is provided, periodic boundary conditions will be used. The `unitcell` can
  be a vector of Orthorhombic box sides, or an actual unitcell matrix for general cells. 
- If `unicell` is not provide (value `nothing`), no periodic boundary conditions will
  be considered. 

## Examples

Here the neighborlist structure is constructed for the first time, and used
to compute the neighbor lists with the mutating `neighborlist!` function:

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 10^4);

julia> system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1]) 
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{OrthorhombicCell, 3, Float64, Float64, 9}
Current list buffer size: 0

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

```julia-repl
julia> x = rand(SVector{3,Float64}, 10^4);

julia> system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1]);

julia> @time neighborlist!(system);
  0.008004 seconds (228 allocations: 16.728 MiB)

julia> update!(system, rand(SVector{3,Float64}, 10^4); cutoff = 0.1, unitcell = [1,1,1]);

julia> @time neighborlist!(system);
  0.024811 seconds (167 allocations: 7.887 MiB)

julia> update!(system, rand(SVector{3,Float64}, 10^4); cutoff = 0.1, unitcell = [1,1,1]);

julia> @time neighborlist!(system);
  0.005213 seconds (164 allocations: 1.439 MiB)

julia> update!(system, rand(SVector{3,Float64}, 10^4); cutoff = 0.1, unitcell = [1,1,1]);

julia> @time neighborlist!(system);
  0.005276 seconds (162 allocations: 15.359 KiB)

```

"""
function InPlaceNeighborList(;
    x::AbstractVecOrMat,
    y::Union{AbstractVecOrMat,Nothing}=nothing,
    cutoff,
    unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
    parallel::Bool=true,
    show_progress::Bool=false,
    nbatches=(0, 0)
)
    T = cutoff isa Integer ? Float64 : eltype(cutoff)
    if isnothing(y)
        if isnothing(unitcell)
            unitcell = limits(x)
        end
        box = Box(unitcell, cutoff)
        cl = CellList(x, box, parallel=parallel, nbatches=nbatches)
        aux = AuxThreaded(cl)
    else
        if isnothing(unitcell)
            unitcell = limits(x, y)
        end
        box = Box(unitcell, cutoff)
        cl = CellList(x, y, box, parallel=parallel, nbatches=nbatches)
        aux = AuxThreaded(cl)
    end
    nb = NeighborList{T}(0, Vector{Tuple{Int,Int,T}}[])
    nb_threaded = [copy(nb) for _ in 1:CellListMap.nbatches(cl, :map)]
    return InPlaceNeighborList(box, cl, aux, nb, nb_threaded, parallel, show_progress)
end

"""
    update!(system::InPlaceNeighborList, x::AbstractVecOrMat; cutoff=nothing, unitcell=nothing)
    update!(system::InPlaceNeighborList, x::AbstractVecOrMat, y::AbstractVecOrMat; cutoff=nothing, unitcell=nothing)

Updates a `InPlaceNeighborList` system, by updating the coordinates, cutoff, and unitcell.

## Examples

### For self-pairs computations

```julia-repl
julia> x = rand(SVector{3,Float64}, 10^3);

julia> system = InPlaceNeighborList(x=x; cutoff=0.1)
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{NonPeriodicCell, 3, Float64, Float64, 9}
Current list buffer size: 0

julia> neighborlist!(system);

julia> new_x = rand(SVector{3,Float64}, 10^3);

julia> update!(system, new_x; cutoff = 0.05)
InPlaceNeighborList with types: 
CellList{3, Float64}
Box{NonPeriodicCell, 3, Float64, Float64, 9}
Current list buffer size: 1826

julia> neighborlist!(system)
224-element Vector{Tuple{Int64, Int64, Float64}}:
 (25, 486, 0.03897345036790646)
 ⋮
 (723, 533, 0.04795768478723409)
 (868, 920, 0.042087156715720137)
```

"""
function update!(
    system::InPlaceNeighborList{<:Box{UnitCellType},C},
    x::AbstractVecOrMat;
    cutoff=nothing, unitcell=nothing
) where {UnitCellType,C<:CellList}
    if UnitCellType == NonPeriodicCell
        isnothing(unitcell) || throw(ArgumentError("Cannot set unitcell for NonPeriodicCell."))
        system.box = update_box(system.box; unitcell=limits(x), cutoff=cutoff)
    else
        system.box = update_box(system.box; unitcell=unitcell, cutoff=cutoff)
    end
    system.cl = UpdateCellList!(x, system.box, system.cl, system.aux, parallel=system.parallel)
    return system
end

#
# update system for cross-computations
#
function update!(
    system::InPlaceNeighborList{<:Box{UnitCellType},C},
    x::AbstractVecOrMat,
    y::AbstractVecOrMat;
    cutoff=nothing, unitcell=nothing
) where {UnitCellType,C<:CellListPair}
    if UnitCellType == NonPeriodicCell
        isnothing(unitcell) || throw(ArgumentError("Cannot set unitcell for NonPeriodicCell."))
        system.box = update_box(system.box; unitcell=limits(x, y), cutoff=cutoff)
    else
        system.box = update_box(system.box; unitcell=unitcell, cutoff=cutoff)
    end
    system.cl = UpdateCellList!(x, y, system.box, system.cl, system.aux; parallel=system.parallel)
    return system
end

@testitem "InPlaceNeighborLists Updates" begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra: diag
    import CellListMap: _sides_from_limits

    # Non-periodic systems
    x = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, cutoff=0.1)
    @test diag(system.box.input_unit_cell.matrix) == _sides_from_limits(limits(x), 0.1)
    x = rand(SVector{3,Float64}, 10^3)
    update!(system, x)
    @test system.box.cutoff == 0.1
    update!(system, x; cutoff=0.05)
    @test system.box.cutoff == 0.05
    @test diag(system.box.input_unit_cell.matrix) == _sides_from_limits(limits(x),0.05)

    x = rand(SVector{3,Float64}, 10^3)
    y = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, y=y, cutoff=0.1)
    @test diag(system.box.input_unit_cell.matrix) ≈ _sides_from_limits(limits(x, y), 0.1)
    x = rand(SVector{3,Float64}, 10^3)
    y = rand(SVector{3,Float64}, 10^3)
    update!(system, x, y)
    @test system.box.cutoff == 0.1
    update!(system, x, y; cutoff=0.05)
    @test system.box.cutoff == 0.05
    @test diag(system.box.input_unit_cell.matrix) ≈ _sides_from_limits(limits(x, y), 0.05)

    # Orthorhombic systems
    x = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1, 1, 1])
    update!(system, x)
    @test system.box.cutoff == 0.1
    update!(system, x; cutoff=0.05)
    @test system.box.cutoff == 0.05
    update!(system, x; cutoff=0.05, unitcell=[2, 2, 2])
    @test (system.box.cutoff, system.box.input_unit_cell.matrix) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    system = InPlaceNeighborList(x=x, y=y, cutoff=0.1, unitcell=[1, 1, 1])
    update!(system, x, y)
    @test system.box.cutoff == 0.1
    update!(system, x, y; cutoff=0.05)
    @test system.box.cutoff == 0.05
    update!(system, x, y; cutoff=0.05, unitcell=[2, 2, 2])
    @test (system.box.cutoff, system.box.input_unit_cell.matrix) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    # Triclinic systems
    x = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1 0 0; 0 1 0; 0 0 1])
    update!(system, x)
    @test system.box.cutoff == 0.1
    update!(system, x; cutoff=0.05)
    @test system.box.cutoff == 0.05
    update!(system, x; cutoff=0.05, unitcell=[2 0 0; 0 2 0; 0 0 2])
    @test (system.box.cutoff, system.box.input_unit_cell.matrix) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    system = InPlaceNeighborList(x=x, y=y, cutoff=0.1, unitcell=[1 0 0; 0 1 0; 0 0 1])
    update!(system, x, y)
    @test system.box.cutoff == 0.1
    update!(system, x, y; cutoff=0.05)
    @test system.box.cutoff == 0.05
    update!(system, x, y; cutoff=0.05, unitcell=[2 0 0; 0 2 0; 0 0 2])
    @test (system.box.cutoff, system.box.input_unit_cell.matrix) == (0.05, [2 0 0; 0 2 0; 0 0 2])

end

function Base.show(io::IO, ::MIME"text/plain", system::InPlaceNeighborList)
    _print(io, "InPlaceNeighborList with types: \n")
    _print(io, typeof(system.cl), "\n")
    _print(io, typeof(system.box), "\n")
    _print(io, "Current list buffer size: $(length(system.nb.list))")
end

function neighborlist!(system::InPlaceNeighborList)
    # Empty lists and auxiliary threaded arrays
    empty!(system.nb)
    for i in eachindex(system.nb_threaded)
        empty!(system.nb_threaded[i])
    end
    # Compute the neighbor lists
    map_pairwise!(
        (x, y, i, j, d2, nb) -> push_pair!(i, j, d2, nb),
        system.nb, system.box, system.cl,
        reduce=reduce_lists,
        parallel=system.parallel,
        output_threaded=system.nb_threaded,
        show_progress=system.show_progress
    )
    # need to resize here to return the correct number of pairs for serial runs
    # (this resizing is redundant for parallel runs, since it occurs at the reduction function
    # before updating)
    !system.parallel && resize!(system.nb, system.nb.n)
    return system.nb.list
end

@testitem "InPlaceNeighborList vs. NearestNeighbors" begin

    using CellListMap
    using CellListMap.TestingNeighborLists
    using NearestNeighbors

    for N in [2, 3]

        x = rand(N, 500)
        r = 0.1
        nb = nl_NN(BallTree, inrange, x, x, r)
        system = InPlaceNeighborList(x=x, cutoff=r)
        cl = neighborlist!(system)
        @test is_unique(cl; self=true)
        @test compare_nb_lists(cl, nb, x, r)[1]
        # Test system updating for self-lists
        r = 0.05
        new_x = rand(N, 450)
        nb = nl_NN(BallTree, inrange, new_x, new_x, r)
        update!(system, new_x; cutoff=r)
        cl = neighborlist!(system)
        @test is_unique(cl; self=true)
        @test compare_nb_lists(cl, nb, x, r)[1]

        # Test system updating for cross-lists
        x = rand(N, 500)
        y = rand(N, 1000)
        r = 0.1
        nb = nl_NN(BallTree, inrange, x, y, r)
        system = InPlaceNeighborList(x=x, y=y, cutoff=r)
        cl = neighborlist!(system)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        r = 0.05
        new_x = rand(N, 500)
        new_y = rand(N, 831)
        nb = nl_NN(BallTree, inrange, new_x, new_y, r)
        update!(system, new_x, new_y; cutoff=r)
        cl = neighborlist!(system)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]

    end

end

@testitem "Allocations" setup=[AllocTest] begin
    using CellListMap
    using StaticArrays
    using BenchmarkTools
    using .AllocTest: Allocs

    #
    # Single set of particles
    #

    # Periodic systems
    x = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1, 1, 1], parallel=false)
    neighborlist!(system)
    x = rand(SVector{3,Float64}, 10^3)
    allocs = @ballocated update!($system, $x) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x; cutoff=0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)

    # Non-Periodic systems
    x = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, cutoff=0.1, parallel=false)
    neighborlist!(system)
    x = rand(SVector{3,Float64}, 10^3)
    allocs = @ballocated update!($system, $x) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x; cutoff=0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)

    #
    # Two sets of particles
    #

    # Periodic systems
    y = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, y=y, cutoff=0.1, unitcell=[1, 1, 1], parallel=false)
    neighborlist!(system)
    x = rand(SVector{3,Float64}, 10^3)
    y = rand(SVector{3,Float64}, 10^3)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x, $y) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x, $y; cutoff=0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)

    # Non-Periodic systems
    y = rand(SVector{3,Float64}, 10^3)
    system = InPlaceNeighborList(x=x, y=y, cutoff=0.1, parallel=false)
    neighborlist!(system)
    x = rand(SVector{3,Float64}, 10^3)
    y = rand(SVector{3,Float64}, 10^3)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x, $y) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system, $x, $y; cutoff=0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)

end

"""
    neighborlist(x, cutoff; unitcell=nothing, parallel=true, show_progress=false)

Computes the list of pairs of particles in `x` which are closer to each other than `cutoff`.
If the keyword parameter `unitcell` is provided (as a vector of sides or a general unit cell
matrix, periodic boundary conditions are considered). 

!!! note
    The order of the pairs in the output of `neighborlist` is not guaranteed,
    and may change, in particular, in parallel runs.

## Example

Compute the neighborlist between within a set Argon atoms, considering the system
non-periodic (do not provide a `unitcell`):

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file));

julia> neighborlist(x, 8.0; parallel=false)
857-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 20, 3.163779543520692)
 (1, 61, 4.0886518560523095)
 (1, 67, 5.939772807102978)
 (1, 80, 2.457228927063981)
 (1, 94, 5.394713986857875)
 (13, 15, 2.678764267344178)
 (13, 41, 4.408015539900014)
 (13, 44, 6.960112211739117)
 (13, 61, 5.939197673086828)
 (13, 64, 4.560755858407684)
 ⋮
 (46, 18, 6.114385414741209)
 (46, 51, 7.999472795128439)
 (51, 68, 2.200357470957844)
 (51, 22, 6.638020940009152)
 (54, 45, 4.423308377221737)
 (73, 78, 2.853611220891874)
 (73, 88, 6.078711047582372)
 (78, 88, 7.006116541993863)
 (88, 54, 7.933654076149277)
```

And now, considering the system periodic:

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file));

julia> neighborlist(x, 8.0; unitcell = [21.0, 21.0, 21.0], parallel=false)
1143-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 7, 3.36387559222989)
 (1, 20, 3.163779543520693)
 (1, 47, 6.243868272153088)
 (1, 63, 7.017762962654125)
 (1, 74, 7.976895636774997)
 (1, 79, 3.177028328485598)
 (1, 94, 5.394713986857875)
 (1, 95, 5.4248765884580274)
 (7, 20, 3.3995637955478935)
 (7, 26, 7.96292025578556)
 ⋮
 (57, 34, 6.536566147450816)
 (57, 84, 7.225401442134547)
 (57, 88, 7.971591246420004)
 (68, 14, 5.2021891545771375)
 (68, 34, 3.955899012866733)
 (68, 84, 5.650943284089833)
 (68, 88, 7.254140403934848)
 (68, 38, 7.4092885623384905)
 (68, 90, 7.875801229081395)
```

"""
function neighborlist(
    x, cutoff;
    unitcell=nothing,
    parallel=true,
    show_progress=false,
    nbatches=(0, 0)
)
    system = InPlaceNeighborList(;
        x=x,
        cutoff=cutoff,
        unitcell=unitcell,
        parallel=parallel,
        show_progress=show_progress,
        nbatches=nbatches
    )
    return neighborlist!(system)
end

"""
    neighborlist(
        x, y, cutoff; 
        unitcell=nothing, 
        parallel=true, 
        show_progress=false, 
        nbatches=(0,0)
    )

Computes the list of pairs of particles of `x` which are closer than `r` to the particles of `y`. 

!!! note
    The order of the pairs in the output of `neighborlist!` is not guaranteed,
    and may change, in particular, in parallel runs.
    
## Examples

Compute the neighborlist between two sets of Argon atoms, considering the system
non-periodic (do not provide a `unitcell`):

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file, "index <= 50"));

julia> y = coor(read_pdb(CellListMap.argon_pdb_file, "index > 50"));

julia> CellListMap.neighborlist(x, y, 8.0; parallel=false)
439-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 11, 4.088651646755291)
 (1, 17, 5.939772435456664)
 (1, 30, 2.4572288423012236)
 (1, 44, 5.394714484195586)
 (13, 11, 5.9391977223435495)
 (13, 14, 4.560755938642345)
 (13, 17, 5.323270872969311)
 (13, 31, 4.201549872989818)
 (15, 11, 6.710523644838785)
 (15, 14, 6.58933933286106)
 ⋮
 (46, 29, 7.402029970260964)
 (46, 50, 4.926250116154994)
 (46, 5, 6.738722573577668)
 (46, 12, 6.363161177968381)
 (46, 22, 5.082701606032681)
 (46, 41, 4.531261514830008)
 (46, 45, 4.251787182648767)
 (46, 48, 4.9269093987894745)
 (46, 1, 7.99947286297016)
```

Now, considering the system periodic:

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file, "index <= 50"));

julia> y = coor(read_pdb(CellListMap.argon_pdb_file, "index > 50"));

julia> CellListMap.neighborlist(x, y, 8.0; unitcell = [21.0, 21.0, 21.0], parallel=false)
584-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 13, 7.0177634180502215)
 (1, 24, 7.97689645513632)
 (1, 29, 3.177029085967527)
 (1, 44, 5.3947144841955845)
 (1, 45, 5.424876392996784)
 (7, 13, 5.245861876160856)
 (7, 24, 7.56131349268393)
 (7, 29, 2.2629708276336706)
 (7, 44, 7.7484227088218764)
 (7, 45, 3.3104152017215167)
 ⋮
 (45, 36, 7.530346972366768)
 (18, 5, 6.8561687188298315)
 (18, 12, 4.461142949385965)
 (18, 41, 5.01835758100573)
 (45, 16, 7.493858435501202)
 (45, 25, 6.0791278577215815)
 (45, 24, 6.97677144471301)
 (18, 10, 6.9654396670725)
 (18, 37, 6.222988130894417)
```

"""
function neighborlist(
    x, y, cutoff;
    unitcell=nothing,
    parallel=true,
    show_progress=false,
    nbatches=(0, 0)
)
    system = InPlaceNeighborList(
        x=x,
        y=y,
        cutoff=cutoff,
        unitcell=unitcell,
        parallel=parallel,
        show_progress=show_progress,
        nbatches=nbatches
    )
    return neighborlist!(system)
end

@testitem "Neighborlist - pathological" begin
    using CellListMap
    using CellListMap.TestingNeighborLists
    using StaticArrays

    @test neighborlist([[0.0, 0.0, 1.0], [0.0, 0.0, 10.0], [0.0, 0.0, 7.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0, 1.0], [0.0, 0.0, 10.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 1.0], [0.0, 10.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 1.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0, 0.0]], 2.0) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0]], 1.0; unitcell=[2.0, 2.0] .+ nextfloat(1.0)) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0], [0.0, 1.0]], 1.0; unitcell=[2.0, 2.0] .+ nextfloat(1.0)) in ([(1, 2, 1.0)],[(2, 1, 1.0)])
    @test neighborlist([[0.0, 0.0], [0.0, 1.0]], prevfloat(1.0); unitcell=[2.0, 2.0]) == Tuple{Int64,Int64,Float64}[]
    @test neighborlist([[0.0, 0.0], [0.0, 1.0] .+ nextfloat(1.0)], prevfloat(1.0); unitcell=[2.0, 2.0]) in ([(1, 2, 0.9999999999999998)],[(2, 1, 0.9999999999999998)])

    # Some pathological cases related to bug 84
    l = SVector{3, Float32}[[0.0, 0.0, 0.0], [0.154, 1.136, -1.827], [-1.16, 1.868, 4.519], [-0.089, 2.07, 4.463],  [0.462, -0.512, 5.473]]
    nl = neighborlist(l, 7.0) 
    @test is_unique(nl; self=true)
    lr = Ref(x_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)
    lr = Ref(y_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)
    lr = Ref(z_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)
    lr = Ref(z_rotation(π/2) * y_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)
    lr = Ref(z_rotation(π/2) * x_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)
    lr = Ref(y_rotation(π/2) * x_rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)

    # in 2D
    rotation(x) = @SMatrix[ cos(x) sin(x); -sin(x) cos(x)]

    l = SVector{2, Float32}[[0.0, 0.0], [0.0, -2.0], [-0.1, 5.0],  [0.0, 5.5]]
    nl = neighborlist(l, 7.0) 
    @test is_unique(nl; self=true)
    lr = Ref(rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)

    l = SVector{2, Float32}[[0.0, 0.0], [-0.1, 5.0]]
    nl = neighborlist(l, 7.0; unitcell=[14.01, 14.51])
    @test length(nl) == 1
    l = Ref(rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)

    l = SVector{2, Float64}[[0.0, 0.0], [-1, 0.0]]
    unitcell = [14.01, 14.02]
    nl = neighborlist(l, 5.0; unitcell=unitcell)
    @test length(nl) == 1
    l = Ref(rotation(π/2)) .* l
    nr = neighborlist(l, 7.0) 
    @test is_unique(nr; self=true)

    unitcell=[1.0,1.0]
    for x in [nextfloat(0.1),prevfloat(0.9)]
        local l, nl, lr
        l = [[0.0,0.0],[x,0.0]] 
        nl = neighborlist(l, 0.1; unitcell=unitcell)
        @test length(nl) == 0
        lr = Ref(rotation(π/2)) .* l
        nl = neighborlist(l, 0.1; unitcell=unitcell)
        @test length(nl) == 0
    end
    for x in [-0.1,0.1,0.9]
        local l, nl, lr
        l = [[0.0,0.0],[x,0.0]] 
        nl = neighborlist(l, 0.1; unitcell=unitcell)
        @test length(nl) == 1
        lr = Ref(rotation(π/2)) .* l
        nl = neighborlist(l, 0.1; unitcell=unitcell)
        @test length(nl) == 1
    end

    # allow cutoff as an integer, promoting it to Float64
    x = [[1,2], [3,4]]
    nb = neighborlist(x, 3)
    @test length(nb) == 1
    @test nb isa Vector{Tuple{Int64, Int64, Float64}}

end

@testitem "Neighborlist with units" begin
    using CellListMap
    using Unitful
    using StaticArrays

    positions = [SVector(0.1, 0.0, 0.0), SVector(0.11, 0.01, 0.01) ]u"nm"
    cutoff = 0.1u"nm"
    nb = neighborlist(positions, cutoff)
    @test unit(nb[1][3]) == u"nm"

    # and with boundary coordinates (to test the fix for upper boundary shifts)
    l = [SVector(0.0, 0.0)u"nm", SVector(-1, 0.0)u"nm"]
    unitcell = [14.01, 14.02]u"nm"
    nl = neighborlist(l, 7.0u"nm")
    @test length(nl) == 1
    @test nl[1][3] ≈ 1.0u"nm"

end

@testitem "Compare with NearestNeighbors" begin

    using CellListMap
    using CellListMap.TestingNeighborLists
    using StaticArrays
    using NearestNeighbors

    r = 0.1

    for N in [2, 3]

        #
        # Using vectors as input
        #

        # With y smaller than x
        x = [rand(SVector{N,Float64}) for _ in 1:500]
        y = [rand(SVector{N,Float64}) for _ in 1:250]

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(x, r)
        @test is_unique(cl; self=true)
        @test compare_nb_lists(cl, nb, x, r)[1]

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(y, x, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # with x smaller than y
        x = [rand(SVector{N,Float64}) for _ in 1:500]
        y = [rand(SVector{N,Float64}) for _ in 1:1000]
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(y, x, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # Using matrices as input
        x = rand(N, 1000)
        y = rand(N, 500)

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(x, r)
        @test is_unique(cl; self=true)
        @test compare_nb_lists(cl, nb, x, r)[1]

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(y, x, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # with x smaller than y
        x = rand(N, 500)
        y = rand(N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(y, x, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # Check random coordinates to test the limits more thoroughly
        check_random_NN = true
        for i in 1:500
            x = rand(SVector{N,Float64}, 100)
            y = rand(SVector{N,Float64}, 50)
            nb = nl_NN(BallTree, inrange, x, y, r)
            cl = CellListMap.neighborlist(x, y, r)
            @test is_unique(cl; self=false)
            check_random_NN = compare_nb_lists(cl, nb, x, y, r)[1]
        end
        @test check_random_NN

        # with different types
        x = rand(Float32, N, 500)
        y = rand(Float32, N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(y, x, r)
        @test is_unique(cl; self=false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

    end

end

@testitem "list buffer reduction" begin
    using CellListMap
    using StaticArrays
    x = [ SVector{3,Float64}(0,0,0), SVector{3,Float64}(0,0,0.05) ];
    system = InPlaceNeighborList(x=x, cutoff=0.1, unitcell=[1,1,1], parallel=false)
    list0 = neighborlist!(system) # correct
    @test length(list0) == 1
    xnew = [ SVector{3,Float64}(0,0,0), SVector{3,Float64}(0,0,0.2) ];
    update!(system, xnew)
    list1 = neighborlist!(system)
    @test length(list1) == 0
end


#
# some auxiliary functions for testing neighbor lists
#
include("../test/TestingNeighborLists.jl")