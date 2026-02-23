@testitem "InPlaceNeighborList vs. NearestNeighbors" setup = [TestingNeighborLists] begin

    using CellListMap
    using CellListMap: get_cutoff, get_unitcell
    using NearestNeighbors

    for N in [2, 3]

        x = rand(N, 500)
        r = 0.1
        nb = nl_NN(BallTree, inrange, x, x, r)
        system = InPlaceNeighborList(xpositions = x, cutoff = r)
        cl = neighborlist!(system)
        @test is_unique(cl; self = true)
        @test compare_nb_lists(cl, nb, x, r)[1]
        # Test system updating for self-lists
        r = 0.05
        new_x = rand(N, 450)
        nb = nl_NN(BallTree, inrange, new_x, new_x, r)
        update!(system; xpositions=new_x, cutoff = r)
        cl = neighborlist!(system)
        @test is_unique(cl; self = true)
        @test compare_nb_lists(cl, nb, new_x, r)[1]

        # Test system updating for cross-lists
        x = rand(N, 500)
        y = rand(N, 1000)
        r = 0.1
        nb = nl_NN(BallTree, inrange, x, y, r)
        system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = r)
        cl = neighborlist!(system)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        r = 0.05
        new_x = rand(N, 500)
        new_y = rand(N, 831)
        nb = nl_NN(BallTree, inrange, new_x, new_y, r)
        update!(system; xpositions=new_x, ypositions=new_y, cutoff = r)
        cl = neighborlist!(system)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, new_x, new_y, r)[1]

    end

end

@testitem "InPlaceNeighborLists Updates" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    using LinearAlgebra: diag
    using CellListMap: _sides_from_limits,
                        get_cutoff,
                        get_unitcell

    # Non-periodic systems
    x = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1)
    update!(system; xpositions=x)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    @test diag(get_unitcell(system)) == _sides_from_limits(CellListMap.limits(PSP(x)), 0.05)

    x = rand(SVector{3, Float64}, 10^3)
    y = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = 0.1)
    @test diag(get_unitcell(system)) ≈ _sides_from_limits(CellListMap.limits(PSP(x), PSP(y)), 0.1)
    x = rand(SVector{3, Float64}, 10^3)
    y = rand(SVector{3, Float64}, 10^3)
    update!(system; xpositions=x, ypositions=y)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, ypositions=y, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    @test diag(get_unitcell(system)) ≈ _sides_from_limits(CellListMap.limits(PSP(x), PSP(y)), 0.05)

    # Orthorhombic systems
    x = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1, unitcell = [1, 1, 1])
    update!(system; xpositions=x)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    update!(system; xpositions=x, cutoff = 0.05, unitcell = [2, 2, 2])
    @test (get_cutoff(system), get_unitcell(system)) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1, 1, 1])
    update!(system; xpositions=x, ypositions=y)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, ypositions=y, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    update!(system; xpositions=x, ypositions=y, cutoff = 0.05, unitcell = [2, 2, 2])
    @test (get_cutoff(system), get_unitcell(system)) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    # Triclinic systems
    x = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1, unitcell = [1 0 0; 0 1 0; 0 0 1])
    update!(system; xpositions=x)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    update!(system; xpositions=x, cutoff = 0.05, unitcell = [2 0 0; 0 2 0; 0 0 2])
    @test (get_cutoff(system), get_unitcell(system)) == (0.05, [2 0 0; 0 2 0; 0 0 2])

    system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1 0 0; 0 1 0; 0 0 1])
    update!(system; xpositions=x, ypositions=y)
    @test get_cutoff(system) == 0.1
    update!(system; xpositions=x, ypositions=y, cutoff = 0.05)
    @test get_cutoff(system) == 0.05
    update!(system; xpositions=x, ypositions=y, cutoff = 0.05, unitcell = [2 0 0; 0 2 0; 0 0 2])
    @test (get_cutoff(system), get_unitcell(system)) == (0.05, [2 0 0; 0 2 0; 0 0 2])

end

@testitem "Allocations" setup = [AllocTest] begin
    using CellListMap
    using StaticArrays
    using BenchmarkTools

    #
    # Single set of particles
    #

    # Periodic systems
    x = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1, unitcell = [1, 1, 1], parallel = false)
    neighborlist!(system)
    x = rand(SVector{3, Float64}, 10^3)
    allocs = @ballocated update!($system; xpositions=$x) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, cutoff = 0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)

    # Non-Periodic systems
    x = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1, parallel = false)
    neighborlist!(system)
    x = rand(SVector{3, Float64}, 10^3)
    allocs = @ballocated update!($system; xpositions=$x) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, cutoff = 0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)

    #
    # Two sets of particles
    #

    # Periodic systems
    y = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1, 1, 1], parallel = false)
    neighborlist!(system)
    x = rand(SVector{3, Float64}, 10^3)
    y = rand(SVector{3, Float64}, 10^3)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, ypositions=$y) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, ypositions=$y, cutoff = 0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)

    # Non-Periodic systems
    y = rand(SVector{3, Float64}, 10^3)
    system = InPlaceNeighborList(xpositions = x, ypositions = y, cutoff = 0.1, parallel = false)
    neighborlist!(system)
    x = rand(SVector{3, Float64}, 10^3)
    y = rand(SVector{3, Float64}, 10^3)
    allocs = @ballocated neighborlist!($system) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, ypositions=$y) evals = 1 samples = 1
    @test allocs == Allocs(0)
    allocs = @ballocated update!($system; xpositions=$x, ypositions=$y, cutoff = 0.2) evals = 1 samples = 1
    @test allocs == Allocs(0)

end

@testitem "Neighborlist - pathological" setup = [TestingNeighborLists] begin
    using CellListMap
    using StaticArrays

    @test neighborlist(xpositions=[[0.0, 0.0, 1.0], [0.0, 0.0, 10.0], [0.0, 0.0, 7.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0, 1.0], [0.0, 0.0, 10.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 1.0], [0.0, 10.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 1.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0, 0.0]], cutoff=2.0) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0]], cutoff=1.0, unitcell=[2.0, 2.0] .+ nextfloat(1.0)) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0], [0.0, 1.0]], cutoff=1.0, unitcell=[2.0, 2.0] .+ nextfloat(1.0)) in ([(1, 2, 1.0)], [(2, 1, 1.0)])
    @test neighborlist(xpositions=[[0.0, 0.0], [0.0, 1.0]], cutoff=prevfloat(1.0), unitcell=[2.0, 2.0]) == Tuple{Int64, Int64, Float64}[]
    @test neighborlist(xpositions=[[0.0, 0.0], [0.0, 1.0] .+ nextfloat(1.0)], cutoff=prevfloat(1.0), unitcell=[2.0, 2.0]) in ([(1, 2, 0.9999999999999998)], [(2, 1, 0.9999999999999998)])

    # Some pathological cases related to bug 84
    l = SVector{3, Float32}[[0.0, 0.0, 0.0], [0.154, 1.136, -1.827], [-1.16, 1.868, 4.519], [-0.089, 2.07, 4.463], [0.462, -0.512, 5.473]]
    nl = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nl; self = true)
    lr = Ref(x_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)
    lr = Ref(y_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)
    lr = Ref(z_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)
    lr = Ref(z_rotation(π / 2) * y_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)
    lr = Ref(z_rotation(π / 2) * x_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)
    lr = Ref(y_rotation(π / 2) * x_rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)

    # in 2D
    rotation(x) = @SMatrix[cos(x) sin(x); -sin(x) cos(x)]

    l = SVector{2, Float32}[[0.0, 0.0], [0.0, -2.0], [-0.1, 5.0], [0.0, 5.5]]
    nl = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nl; self = true)
    lr = Ref(rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)

    l = SVector{2, Float32}[[0.0, 0.0], [-0.1, 5.0]]
    nl = neighborlist(xpositions=l, cutoff=7.0, unitcell=[14.01, 14.51])
    @test length(nl) == 1
    l = Ref(rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)

    l = SVector{2, Float64}[[0.0, 0.0], [-1, 0.0]]
    unitcell = [14.01, 14.02]
    nl = neighborlist(xpositions=l, cutoff=5.0, unitcell=unitcell)
    @test length(nl) == 1
    l = Ref(rotation(π / 2)) .* l
    nr = neighborlist(xpositions=l, cutoff=7.0)
    @test is_unique(nr; self = true)

    unitcell = [1.0, 1.0]
    for x in [nextfloat(0.1), prevfloat(0.9)]
        local l, nl, lr
        l = [[0.0, 0.0], [x, 0.0]]
        nl = neighborlist(xpositions=l, cutoff=0.1, unitcell=unitcell)
        @test length(nl) == 0
        lr = Ref(rotation(π / 2)) .* l
        nl = neighborlist(xpositions=l, cutoff=0.1, unitcell=unitcell)
        @test length(nl) == 0
    end
    for x in [-0.1, 0.1, 0.9]
        local l, nl, lr
        l = [[0.0, 0.0], [x, 0.0]]
        nl = neighborlist(xpositions=l, cutoff=0.1, unitcell=unitcell)
        @test length(nl) == 1
        lr = Ref(rotation(π / 2)) .* l
        nl = neighborlist(xpositions=l, cutoff=0.1, unitcell=unitcell)
        @test length(nl) == 1
    end

    # allow cutoff as an integer, promoting it to Float64
    x = [[1, 2], [3, 4]]
    nb = neighborlist(xpositions=x, cutoff=3)
    @test length(nb) == 1
    @test nb isa Vector{Tuple{Int64, Int64, Float64}}

end

@testitem "Neighborlist with units" begin
    using CellListMap
    using Unitful
    using StaticArrays

    positions = [SVector(0.1, 0.0, 0.0), SVector(0.11, 0.01, 0.01)]u"nm"
    cutoff = 0.1u"nm"
    nb = neighborlist(xpositions=positions, cutoff=cutoff)
    @test unit(nb[1][3]) == u"nm"

    # and with boundary coordinates (to test the fix for upper boundary shifts)
    l = [SVector(0.0, 0.0)u"nm", SVector(-1, 0.0)u"nm"]
    unitcell = [14.01, 14.02]u"nm"
    nl = neighborlist(xpositions=l, cutoff=7.0u"nm")
    @test length(nl) == 1
    @test nl[1][3] ≈ 1.0u"nm"

end

@testitem "Compare with NearestNeighbors" setup = [TestingNeighborLists] begin

    using CellListMap
    using StaticArrays
    using NearestNeighbors

    r = 0.1

    for N in [2, 3]

        #
        # Using vectors as input
        #

        # With y smaller than x
        x = [rand(SVector{N, Float64}) for _ in 1:500]
        y = [rand(SVector{N, Float64}) for _ in 1:250]

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(xpositions=x, cutoff=r)
        @test is_unique(cl; self = true)
        @test compare_nb_lists(cl, nb, x, r)[1]

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(xpositions=y, ypositions=x, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # with x smaller than y
        x = [rand(SVector{N, Float64}) for _ in 1:500]
        y = [rand(SVector{N, Float64}) for _ in 1:1000]
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(xpositions=y, ypositions=x, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # Using matrices as input
        x = rand(N, 1000)
        y = rand(N, 500)

        nb = nl_NN(BallTree, inrange, x, x, r)
        cl = CellListMap.neighborlist(xpositions=x, cutoff=r)
        @test is_unique(cl; self = true)
        @test compare_nb_lists(cl, nb, x, r)[1]

        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(xpositions=y, ypositions=x, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # with x smaller than y
        x = rand(N, 500)
        y = rand(N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(xpositions=y, ypositions=x, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

        # Check random coordinates to test the limits more thoroughly
        check_random_NN = true
        for i in 1:500
            x = rand(SVector{N, Float64}, 100)
            y = rand(SVector{N, Float64}, 50)
            nb = nl_NN(BallTree, inrange, x, y, r)
            cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
            @test is_unique(cl; self = false)
            check_random_NN = compare_nb_lists(cl, nb, x, y, r)[1]
        end
        @test check_random_NN

        # with different types
        x = rand(Float32, N, 500)
        y = rand(Float32, N, 1000)
        nb = nl_NN(BallTree, inrange, x, y, r)
        cl = CellListMap.neighborlist(xpositions=x, ypositions=y, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, x, y, r)[1]
        nb = nl_NN(BallTree, inrange, y, x, r)
        cl = CellListMap.neighborlist(xpositions=y, ypositions=x, cutoff=r)
        @test is_unique(cl; self = false)
        @test compare_nb_lists(cl, nb, y, x, r)[1]

    end

end

@testitem "list buffer reduction" begin
    using CellListMap
    using StaticArrays
    x = [SVector{3, Float64}(0, 0, 0), SVector{3, Float64}(0, 0, 0.05)]
    system = InPlaceNeighborList(xpositions = x, cutoff = 0.1, unitcell = [1, 1, 1], parallel = false)
    list0 = neighborlist!(system) # correct
    @test length(list0) == 1
    xnew = [SVector{3, Float64}(0, 0, 0), SVector{3, Float64}(0, 0, 0.2)]
    update!(system; xpositions=xnew)
    list1 = neighborlist!(system)
    @test length(list1) == 0
end
