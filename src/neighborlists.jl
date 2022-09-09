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
copy(x::Tuple{Int,Int,T}) where {T} = (x[1], x[2], x[3])

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
    list = resize!(list, npairs)
    @sync for it in eachindex(list_threaded)
        lt = list_threaded[it]
        range = ranges[it]-lt.n+1:ranges[it]
        Threads.@spawn list.list[range] .= @view(lt.list[1:lt.n])
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

"""

$(TYPEDEF)

$(INTERNAL)

Structure that containst the system information for neighborlist computations. All fields are internal.

## Extended help

$(TYPEDFIELDS)

"""
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

```
function InPlaceNeighborList(;
    x::AbstractVecOrMat,
    y::Union{AbstractVecOrMat,Nothing}=nothing,
    cutoff::T,
    unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
    parallel::Bool=true,
    show_progress::Bool=false,
) where {T<:Real}
```

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
    cutoff::T,
    unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
    parallel::Bool=true,
    show_progress::Bool=false,
    autoswap=true,
    nbatches=(0, 0)
) where {T<:Real}
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
        cl = CellList(x, y, box, autoswap=autoswap, parallel=parallel, nbatches=nbatches)
        aux = AuxThreaded(cl)
    end
    nb = NeighborList{T}(0, Vector{Tuple{Int,Int,T}}[])
    nb_threaded = [copy(nb) for _ in 1:CellListMap.nbatches(cl)]
    return InPlaceNeighborList(box, cl, aux, nb, nb_threaded, parallel, show_progress)
end

#
# Function that updates the box if the cutoff or unitcell has changed
#
# System without perioidc boundary conditions: update limits
function update_box!(system::InPlaceNeighborList{<:Box{NonPeriodicCell}}, x, cutoff::Real, unitcell::Nothing)
    system.box = Box(limits(x), cutoff)
    return system
end
function update_box!(system::InPlaceNeighborList{<:Box{NonPeriodicCell}}, x, y, cutoff::Real, unitcell::Nothing)
    system.box = Box(limits(x,y), cutoff)
    return system
end
# Nothing change for systems with periodic boundary conditions: just return
update_box!(system::InPlaceNeighborList{<:Box{<:PeriodicCellType}}, x, cutoff::Nothing, unitcell::Nothing) = system
update_box!(system::InPlaceNeighborList{<:Box{<:PeriodicCellType}}, x, y, cutoff::Nothing, unitcell::Nothing) = system
# Systems with periodic boundary conditions
function update_box!(system::InPlaceNeighborList{<:Box{UnitCellType}}, cutoff::Real, unitcell) where {UnitCellType<:PeriodicCellType}
    if isnothing(cutoff)
        cutoff = system.box.cutoff
    end
    if isnothing(unitcell)
        unitcell = system.box.unit_cell.matrix
        box = Box(unitcell, cutoff, UnitCellType=UnitCellType)
    else
        box = Box(unitcell, cutoff)
    end
    # Error if the type of box was changed: this is not supported
    B = typeof(system.box)
    if typeof(box) != B
        throw(ArgumentError(" The unitcell must be of the same type (Orthorhombic **or** general) as of the initial constructor"))
    end
    system.box = box
    return system
end

#
# update system for self computations
#
function update!(
    system::InPlaceNeighborList{<:Box{UnitCellType},C},
    x::AbstractVecOrMat;
    cutoff=nothing, unitcell=nothing
) where {UnitCellType,C<:CellList}
    update_box!(system, x, cutoff, unitcell)
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
    update_box!(system, x, y, cutoff, unitcell)
    system.cl = UpdateCellList!(x, y, system.box, system.cl, system.aux; parallel=system.parallel)
    return system
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
    empty!.(system.nb_threaded)
    # Compute the neighbor lists
    map_pairwise!(
        (x, y, i, j, d2, nb) -> push_pair!(i, j, d2, nb),
        system.nb, system.box, system.cl,
        reduce=reduce_lists,
        parallel=system.parallel,
        output_threaded=system.nb_threaded,
        show_progress=system.show_progress
    )
    return system.nb.list
end


@testitem "InPlaceNeighborList vs. NearestNeighbors" begin

    using CellListMap
    using CellListMap.TestingNeighborLists
    using NearestNeighbors

    for N in [2,3]
    
        x = rand(N, 500)
        cutoff = 0.1
        nb = nl_NN(BallTree, inrange, x, x, cutoff)
        system = InPlaceNeighborList(x = x, cutoff = cutoff)
        cl = neighborlist!(system)
        @test compare_nb_lists(cl, nb, self=true)
        # Test system updating for self-lists
        cutoff = 0.05
        new_x = rand(N, 450)
        nb = nl_NN(BallTree, inrange, new_x, new_x, cutoff)
        update!(system, new_x; cutoff=cutoff)
        cl = neighborlist!(system)
        @test compare_nb_lists(cl, nb, self=true)

        # Test system updating for cross-lists
        x = rand(N, 500)
        y = rand(N, 1000)
        cutoff = 0.1
        nb = nl_NN(BallTree, inrange, x, y, cutoff)
        system = InPlaceNeighborList(x = x, y = y, cutoff = cutoff)
        cl = neighborlist!(system)
        @test compare_nb_lists(cl, nb, self=false)
        cutoff = 0.05
        new_x = rand(N, 500) 
        new_y = rand(N, 831) 
        nb = nl_NN(BallTree, inrange, new_x, new_y, cutoff)
        update!(system, new_x, new_y; cutoff=cutoff)
        cl = neighborlist!(system)
        @test compare_nb_lists(cl, nb, self=false)

    end

end

"""

```
neighborlist(x, cutoff; unitcell=nothing, parallel=true, show_progress=false)
```

Computes the list of pairs of particles in `x` which are closer to each other than `cutoff`.
If the keyword parameter `unitcell` is provided (as a vector of sides or a general unit cell
matrix, periodic boundary conditions are considered). 

## Example
```julia-repl
julia> using CellListMap

julia> x = [ rand(3) for i in 1:10_000 ];

julia> neighborlist(x,0.05)
24848-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1055, 0.022977369806392412)
 (1, 5086, 0.026650609138167428)
 ⋮
 (9989, 3379, 0.0467653507446483)
 (9989, 5935, 0.02432728985151653)

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

```
neighborlist(
    x, y, cutoff; 
    unitcell=nothing, 
    parallel=true, 
    show_progress=false, 
    autoswap=true,
    nbatches=(0,0))
```

Computes the list of pairs of particles of `x` which are closer than `r` to
the particles of `y`. The `autoswap` option will swap `x` and `y` to try to optimize
the cost of the construction of the cell list. 

## Example
```julia-repl
julia> x = [ rand(3) for i in 1:10_000 ];

julia> y = [ rand(3) for i in 1:1_000 ];

julia> CellListMap.neighborlist(x,y,0.05)
5006-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 269, 0.04770884036497686)
 (25, 892, 0.03850515231540869)
 ⋮
 (9952, 749, 0.048875643578313456)
 (9984, 620, 0.04101242499363183)

```

"""
function neighborlist(
    x, y, cutoff;
    unitcell=nothing,
    parallel=true,
    show_progress=false,
    autoswap=true,
    nbatches=(0, 0)
)
    system = InPlaceNeighborList(
        x=x,
        y=y,
        cutoff=cutoff,
        unitcell=unitcell,
        parallel=parallel,
        show_progress=show_progress,
        autoswap=autoswap,
        nbatches=nbatches
    )
    return neighborlist!(system)
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
        x = [rand(SVector{N,Float64}) for _ in 1:1000]
        y = [rand(SVector{N,Float64}) for _ in 1:500]

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(x, r)
        @test compare_nb_lists(cl, nb, self=true)

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r, autoswap=false)
        @test compare_nb_lists(cl, nb, self=false)
        cl = CellListMap.neighborlist(x, y, r, autoswap=true)
        @test compare_nb_lists(cl, nb, self=false)

        # with x smaller than y
        x = [rand(SVector{N,Float64}) for _ in 1:500]
        y = [rand(SVector{N,Float64}) for _ in 1:1000]
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r, autoswap=false)
        @test compare_nb_lists(cl, nb, self=false)
        cl = CellListMap.neighborlist(x, y, r, autoswap=true)
        @test compare_nb_lists(cl, nb, self=false)

        # Using matrices as input
        x = rand(N, 1000)
        y = rand(N, 500)

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(x, r)
        @test compare_nb_lists(cl, nb, self=true)

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r, autoswap=false)
        @test compare_nb_lists(cl, nb, self=false)
        cl = CellListMap.neighborlist(x, y, r, autoswap=true)
        @test compare_nb_lists(cl, nb, self=false)

        # with x smaller than y
        x = rand(N, 500)
        y = rand(N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r, autoswap=false)
        @test compare_nb_lists(cl, nb, self=false)
        cl = CellListMap.neighborlist(x, y, r, autoswap=true)
        @test compare_nb_lists(cl, nb, self=false)

        # Check random coordinates to test the limits more thoroughly
        check_random_NN = true
        for i in 1:500
            x = rand(SVector{N,Float64}, 100)
            y = rand(SVector{N,Float64}, 50)
            nb = nl_NN(BallTree, inrange, x, y, r)
            cl = CellListMap.neighborlist(x, y, r, autoswap=false)
            check_random_NN = compare_nb_lists(cl, nb, self=false)
        end
        @test check_random_NN

        # with different types
        x = rand(Float32, N, 500)
        y = rand(Float32, N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(x, y, r, autoswap=false)
        @test compare_nb_lists(cl, nb, self=false)
        cl = CellListMap.neighborlist(x, y, r, autoswap=true)
        @test compare_nb_lists(cl, nb, self=false)

    end

end

#
# some auxiliary functions for testing neighbor lists
#
module TestingNeighborLists

    export nl_NN
    export compare_nb_lists

    function nl_NN(BallTree, inrange, x, y, r)
        balltree = BallTree(x)
        return inrange(balltree, y, r, true)
    end

    function compare_nb_lists(list_CL, list_NN; self=false)
        if self
            for (i, list_original) in pairs(list_NN)
                list = filter(!isequal(i), list_original)
                cl = filter(tup -> (tup[1] == i || tup[2] == i), list_CL)
                if length(cl) != length(list)
                    @show i
                    @show length(list), list
                    @show length(cl), cl
                    return false
                end
                for j in list
                    length(findall(tup -> (tup[1] == j || tup[2] == j), cl)) == 1 || return false
                end
            end
        else
            for (i, list) in pairs(list_NN)
                cl = filter(tup -> tup[2] == i, list_CL)
                if length(cl) != length(list)
                    @show i
                    @show length(list), list
                    @show length(cl), cl
                    return false
                end
                for j in list
                    length(findall(tup -> tup[1] == j, cl)) == 1 || return false
                end
            end
        end
        return true
    end

end # module TestingNeighborLists
