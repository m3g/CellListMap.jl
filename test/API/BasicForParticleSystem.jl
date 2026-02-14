@testitem "BasicForParticleSystem" setup = [Testing] begin

    using StaticArrays
    using CellListMap

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    mass = rand(N)

    # Function to be evaluated for each pair: gravitational potential
    function potential(i, j, d2, u, mass)
        d = sqrt(d2)
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end

    # Some simple disjoint set properties
    system = ParticleSystem(
        xpositions = x,
        ypositions = y,
        cutoff = cutoff,
        unitcell = sides,
        output = 0.0,
        output_name = :gravitational_potential
    )
    naive = map_naive!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), 0.0, x, y, CellListMap.Box(sides, cutoff))
    system.parallel = false
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive
    system.parallel = true
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive
    # Test fetching by output name
    @test system.gravitational_potential ≈ naive

    # Same but for non-periodic systems
    system = ParticleSystem(
        xpositions = x,
        ypositions = y,
        cutoff = cutoff,
        output = 0.0,
    )
    naive = map_naive!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), 0.0, x, y, CellListMap.Box(CellListMap.limits(x, y), cutoff))
    system.parallel = false
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive
    system.parallel = true
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive

    # Use matrices as input coordinates
    xmat = zeros(3, N)
    ymat = zeros(3, N)
    for i in 1:N
        xmat[:, i] .= x[i]
        ymat[:, i] .= y[i]
    end
    system = ParticleSystem(
        xpositions = xmat,
        ypositions = ymat,
        cutoff = cutoff,
        unitcell = sides,
        output = 0.0,
    )
    naive = map_naive!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), 0.0, x, y, CellListMap.Box(sides, cutoff))
    system.parallel = false
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive
    system.parallel = true
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive

    # Check different lcell
    system = ParticleSystem(
        xpositions = x,
        ypositions = y,
        cutoff = cutoff,
        unitcell = sides,
        output = 0.0,
        lcell = 3,
    )
    system.parallel = false
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive
    system.parallel = true
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive

    # Test updating of the data on disjoint sets works fine
    for arrays in [
                # static vectors, length(x) > length(y)
                [rand(SVector{2, Float64}, 1000), rand(SVector{2, Float64}, 100)],
                # static vectors, length(x) < length(y)
                [rand(SVector{2, Float64}, 100), rand(SVector{2, Float64}, 1000)],
                # standard vectors, length(x) > length(y)
                [[rand(2) for _ in 1:1000], [rand(2) for _ in 1:100]], # with standard vectors
                # standard vectors, length(x) < length(y)
                [[rand(2) for _ in 1:100], [rand(2) for _ in 1:1000]], # with standard vectors
                # matrices, length(x) > length(y)
                [rand(2, 1000), rand(2, 100)],
                # matrices, length(x) < length(y)
                [rand(2, 100), rand(2, 1000)], # with standard vectors
            ], unitcell in [[1, 1], nothing]
        local x = PSP(arrays[1])
        local y = PSP(arrays[2])
        local system
        system = ParticleSystem(
            xpositions = x,
            ypositions = y,
            cutoff = 0.1,
            output = 0.0,
            unitcell = unitcell,
        )
        uc = isnothing(unitcell) ? CellListMap.limits(x, y) : [1, 1]
        box = CellListMap.Box(uc, 0.1)
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)
        if x isa AbstractVector
            x = PSP(rand(SVector{2, Float64}, length(x) + 100))
            resize!(system.xpositions, length(x))
        else
            # Cannot resize the matrices, so this interface is more limited
            x = PSP(rand(2, size(x, 2)))
        end
        cl = CellListMap.UpdateCellList!(x, y, box, cl)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system.xpositions .= x
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        if y isa AbstractVector
            y = PSP(rand(SVector{2, Float64}, length(y) + 100))
            resize!(system.ypositions, length(y))
        else
            y = PSP(rand(2, size(y, 2)))
        end
        cl = CellListMap.UpdateCellList!(x, y, box, cl)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system.ypositions .= y
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

    end

    # Check if the dimension of the input array is properly errored
    # when the dimension of the unitcell does not match that of SVectors
    # of the input array
    @test_throws DimensionMismatch ParticleSystem(
        positions = rand(SVector{3, Float64}, 100),
        unitcell = [1, 1],
        cutoff = 0.1,
        output = 0.0,
    )
    sys = ParticleSystem(
        positions = rand(SVector{3, Float64}, 100),
        unitcell = [1, 1, 1],
        cutoff = 0.1,
        output = 0.0,
    )
    @test size(sys.unitcell) == (3, 3)

end

@testitem "ParticleSystem - resize output" setup = [Testing] begin
    using StaticArrays
    using CellListMap
    function force(x, y, i, j, d2, f)
        fxy = (x - y) / d2
        f[i] -= fxy
        f[j] += fxy
        return f
    end
    x = rand(SVector{3, Float64}, 100)
    sides = [1, 1, 1]
    cutoff = 0.1
    sys = ParticleSystem(
        positions = x,
        unitcell = sides,
        cutoff = cutoff,
        output = zero(x),
        output_name = :force
    )
    naive = map_naive!(force, zero(x), x, CellListMap.Box(sides, cutoff))
    pairwise!((pair, f) -> force(pair.x, pair.y, pair.i, pair.j, pair.d2, f), sys)
    @test sys.force ≈ naive
    resize!(x, 120)
    resize!(sys.xpositions, 120)
    x[101:end] .= rand(SVector{3, Float64}, 20)
    sys.xpositions .= x
    naive = map_naive!(force, zero(x), x, CellListMap.Box(sides, cutoff))
    resize_output!(sys, 120)
    pairwise!((pair, f) -> force(pair.x, pair.y, pair.i, pair.j, pair.d2, f), sys)
    @test sys.force ≈ naive
    resize!(x, 90)
    resize!(sys.xpositions, 90)
    naive = map_naive!(force, zero(x), x, CellListMap.Box(sides, cutoff))
    resize_output!(sys, 90)
    pairwise!((pair, f) -> force(pair.x, pair.y, pair.i, pair.j, pair.d2, f), sys)
    @test sys.force ≈ naive
end

@testitem "ParticleSystem - automatic updating" setup=[Testing] begin
    using StaticArrays
    using CellListMap

    for unitcell in [[1, 1, 1], nothing]
        #
        # one set systems
        #
        x = PSP(rand(SVector{3, Float64}, 100))
        uc = isnothing(unitcell) ? CellListMap.limits(x) : unitcell
        box = CellListMap.Box(uc, 0.1)
        cl = CellListMap.CellList(x, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system = ParticleSystem(xpositions = x, cutoff = 0.1, output = 0.0, unitcell = unitcell)
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)
        # After pairwise!, the updated flag should be false
        @test system.xpositions.updated[] == false
        # Compute a different property from the same coordinates (lists not updated)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d, 0.0, box, cl)
        @test r ≈ pairwise!((pair, r) -> r += pair.d, system)

        #
        # two-set systems
        #

        #
        # x is smaller
        #
        x = PSP(rand(SVector{3, Float64}, 100))
        y = PSP(rand(SVector{3, Float64}, 1000))

        uc = isnothing(unitcell) ? CellListMap.limits(x, y) : unitcell
        box = CellListMap.Box(uc, 0.1)
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, output = 0.0, unitcell = unitcell)
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # change x coordinates (setindex! sets updated flag automatically)
        x = PSP(rand(SVector{3, Float64}, 100))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system.xpositions .= x
        @test system.xpositions.updated[] == true
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # increase x size (resize! sets updated flag automatically)
        x = PSP(rand(SVector{3, Float64}, 200))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        resize!(system.xpositions, length(x))
        system.xpositions .= x
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        #
        # x is greater
        #
        x = PSP(rand(SVector{3, Float64}, 1000))
        y = PSP(rand(SVector{3, Float64}, 100))

        uc = isnothing(unitcell) ? CellListMap.limits(x, y) : unitcell
        box = CellListMap.Box(uc, 0.1)
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, output = 0.0, unitcell = unitcell)
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # change x coordinates
        x = PSP(rand(SVector{3, Float64}, 1000))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system.xpositions .= x
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # increase x size
        x = PSP(rand(SVector{3, Float64}, 1100))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        resize!(system.xpositions, length(x))
        system.xpositions .= x
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # change y coordinates
        y = PSP(rand(SVector{3, Float64}, 100))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        system.ypositions .= y
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # increase y size
        y = PSP(rand(SVector{3, Float64}, 110))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        resize!(system.ypositions, length(y))
        system.ypositions .= y
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # increase y size beyond x size
        y = PSP(rand(SVector{3, Float64}, 1300))
        cl = CellListMap.CellList(x, y, box)
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += pair.d2, 0.0, box, cl)
        resize!(system.ypositions, length(y))
        system.ypositions .= y
        @test r ≈ pairwise!((pair, r) -> r += pair.d2, system)

        # compute a different property from the same coordinates (lists not updated automatically)
        @test system.xpositions.updated[] == false
        @test system.ypositions.updated[] == false
        r = CellListMap.CellListMap._pairwise!((pair, r) -> r += 2 * pair.d2, 0.0, box, cl)
        @test r ≈ pairwise!((pair, r) -> r += 2 * pair.d2, system)

    end
end

@testitem "cross_x_vs_sys" setup = [Testing] begin
    using CellListMap
    using StaticArrays

    f(_, _, _, _, d2, out) = out += d2
    f(pair::NeighborPair, out) = out += pair.d2
    cutoff = 0.1
    for uc in (nothing, [1, 1, 1], [1 0.2 0; 0.2 1.0 0; 0 0 1 ]), parallel in (false, true)
        x = rand(SVector{3, Float64}, 100)
        # y smaller than x
        y = rand(SVector{3, Float64}, 10)
        box = CellListMap.Box(isnothing(uc) ? CellListMap.limits(x, y) : uc, cutoff)
        naive = map_naive!(f, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, unitcell = uc, parallel = parallel)
        @test naive ≈ pairwise!(f, y, sys)
        # Call again without changing positions (lists should not be updated)
        @test naive ≈ pairwise!(f, y, sys)

        # Update x positions
        sys.xpositions .= rand(SVector{3, Float64}, 100)
        box = CellListMap.Box(isnothing(uc) ? CellListMap.limits(x, y) : uc, cutoff)
        naive = map_naive!(f, 0.0, x, y, box)
        @test naive ≈ pairwise!(f, y, sys)

        # Matrices as inputs
        xmat = stack(x)
        ymat = stack(y)
        sys = ParticleSystem(xpositions = xmat, cutoff = cutoff, output = 0.0, unitcell = uc, parallel = parallel)
        @test naive ≈ pairwise!(f, ymat, sys)
        # Call again without changing positions
        @test naive ≈ pairwise!(f, ymat, sys)

        # y greater than x, and possibly spanning a region outside the box of x
        y = 2 .* rand(SVector{3, Float64}, 100)
        box = CellListMap.Box(isnothing(uc) ? CellListMap.limits(x, y) : uc, cutoff)
        naive = map_naive!(f, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, unitcell = uc, parallel = parallel)
        @test naive ≈ pairwise!(f, y, sys)
        # Call again without changing positions
        @test naive ≈ pairwise!(f, y, sys)

        # here y is completely outside the range of x, thus without PBCs, this is zero
        y = rand(SVector{3, Float64}, 10) .+ Ref(SVector(10.0, 10.0, 10.0))
        box = CellListMap.Box(isnothing(uc) ? CellListMap.limits(x, y) : uc, cutoff)
        naive = map_naive!(f, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, unitcell = uc, parallel = parallel)
        @test naive ≈ pairwise!(f, y, sys)
        # Call again without changing positions
        @test naive ≈ pairwise!(f, y, sys)
    end

    # Test if output_threaded was not provided
    x = rand(SVector{3, Float64}, 100)
    y = rand(SVector{3, Float64}, 10)
    box = CellListMap.Box([1, 1, 1], cutoff)
    cl = CellListMap.CellList(x, box)
    naive = map_naive!(f, 0.0, x, y, box)
    @test naive ≈ CellListMap.CellListMap._pairwise!(f, 0.0, box, y, cl; output_threaded = nothing, parallel = true)

end

# Regression tests for NonPeriodicCell cross-computation with negative/shifted coordinates
@testitem "cross_x_vs_sys_nonperiodic_negative_coords" setup = [Testing] begin
    using CellListMap
    using StaticArrays

    f_sum_d2(_, _, _, _, d2, out) = out += d2
    f_sum_d2(pair::NeighborPair, out) = out += pair.d2
    cutoff = 0.1

    # Test 1: Coordinates spanning negative values (original bug reproduction)
    for parallel in (false, true)
        x = [SVector{3,Float64}(-0.1, -0.2, -0.3) .+ rand(SVector{3,Float64}) for _ in 1:200]
        y = [SVector{3,Float64}(-0.1, -0.2, -0.3) .+ rand(SVector{3,Float64}) for _ in 1:50]
        box = CellListMap.Box(CellListMap.limits(x, y), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        result = pairwise!(f_sum_d2, y, sys)
        @test naive ≈ result
    end

    # Test 2: Coordinates entirely in negative range
    for parallel in (false, true)
        x = [SVector(-5.0, -5.0, -5.0) .+ rand(SVector{3,Float64}) for _ in 1:200]
        y = [SVector(-5.0, -5.0, -5.0) .+ rand(SVector{3,Float64}) for _ in 1:50]
        box = CellListMap.Box(CellListMap.limits(x, y), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test naive ≈ pairwise!(f_sum_d2, y, sys)
    end

    # Test 3: Large positive coordinate offset (far from origin)
    for parallel in (false, true)
        x = [SVector(1000.0, 1000.0, 1000.0) .+ rand(SVector{3,Float64}) for _ in 1:200]
        y = [SVector(1000.0, 1000.0, 1000.0) .+ rand(SVector{3,Float64}) for _ in 1:50]
        box = CellListMap.Box(CellListMap.limits(x, y), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test naive ≈ pairwise!(f_sum_d2, y, sys)
    end

    # Test 4: Cross particles far away from system particles should find no pairs
    for parallel in (false, true)
        x = rand(SVector{3,Float64}, 100)
        y = [SVector(100.0, 100.0, 100.0) .+ rand(SVector{3,Float64}) for _ in 1:20]
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test pairwise!(f_sum_d2, y, sys) == 0.0
    end

    # Test 5: Verify that NeighborPair reports original (untranslated) coordinates
    let # just to avoid warning about scope 
        offset = SVector(-10.0, -20.0, -30.0)
        x = [offset .+ SVector(0.5, 0.5, 0.5)]
        y = [offset .+ SVector(0.5, 0.5, 0.55)]  # close to x[1]
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = Tuple{SVector{3,Float64}, SVector{3,Float64}}[], parallel = false)
        function collect_coords(pair::NeighborPair, out)
            push!(out, (pair.x, pair.y))
            return out
        end
        result = pairwise!(collect_coords, y, sys)
        @test length(result) == 1
        # pair.x should be the cross particle (y[1]) in original coordinates
        # pair.y should be the system particle (x[1]) in original coordinates
        @test result[1][1] ≈ y[1]
        @test result[1][2] ≈ x[1]
    end

    # Test 6: 2D NonPeriodicCell with negative coordinates
    for parallel in (false, true)
        x2d = [SVector(-3.0, -3.0) .+ rand(SVector{2,Float64}) for _ in 1:200]
        y2d = [SVector(-3.0, -3.0) .+ rand(SVector{2,Float64}) for _ in 1:50]
        box2d = CellListMap.Box(CellListMap.limits(x2d, y2d), cutoff)
        naive2d = map_naive!(f_sum_d2, 0.0, x2d, y2d, box2d)
        sys2d = ParticleSystem(xpositions = x2d, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test naive2d ≈ pairwise!(f_sum_d2, y2d, sys2d)
    end

    # Test 7: Self-computation with NonPeriodicCell and negative coordinates (non-regression)
    for parallel in (false, true)
        x = [SVector(-2.0, -2.0, -2.0) .+ rand(SVector{3,Float64}) for _ in 1:200]
        box = CellListMap.Box(CellListMap.limits(x), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test naive ≈ pairwise!(f_sum_d2, sys)
    end

    # Test 8: Mixed signs - x in negative, y in positive, with some pairs within cutoff
    for parallel in (false, true)
        x = [SVector(-0.05, 0.0, 0.0) .+ 0.01 .* rand(SVector{3,Float64}) for _ in 1:50]
        y = [SVector(0.0, 0.0, 0.0) .+ 0.01 .* rand(SVector{3,Float64}) for _ in 1:50]
        box = CellListMap.Box(CellListMap.limits(x, y), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = parallel)
        @test naive ≈ pairwise!(f_sum_d2, y, sys)
    end

    # Test 9: Matrix input with negative coordinates
    let # just to avoid warning about scope 
        x = [SVector(-1.0, -1.0, -1.0) .+ rand(SVector{3,Float64}) for _ in 1:100]
        y = [SVector(-1.0, -1.0, -1.0) .+ rand(SVector{3,Float64}) for _ in 1:30]
        ymat = stack(y)
        box = CellListMap.Box(CellListMap.limits(x, y), cutoff)
        naive = map_naive!(f_sum_d2, 0.0, x, y, box)
        sys = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, parallel = false)
        @test naive ≈ pairwise!(f_sum_d2, ymat, sys)
    end
end

@testitem "lcell and nbatches" setup=[ Testing ] begin
    using StaticArrays
    using CellListMap

    if Threads.nthreads() == 1
        println(
            """

                 WARNING: Ideally, run a multi-threaded test to check the parallel versions.

            """
        )
    end

    # Function to be evaluated for each pair: sum of displacements on x
    f(x, y, avg_dx) = avg_dx + abs(x[1] - y[1])

    x, y, sides, cutoff = pathological_coordinates(2000)
    box = CellListMap.Box(sides, cutoff)
    naive = map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, x, box)

    # Check if changing lcell breaks something
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, lcell = 1)
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, lcell = 3)
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, lcell = 5)
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive

    # Test if changing the number of batches breaks anything
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (3, 5))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (1, 1))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (1, 7))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (7, 1))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (13, 17))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, unitcell = sides, cutoff = cutoff, output = 0.0, nbatches = (4, 16))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive

    # non-periodic system
    box = CellListMap.Box(CellListMap.limits(x), cutoff)
    naive = map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, x, box)
    system = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, lcell = 1)
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, lcell = 3)
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, nbatches = (3, 5))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system = ParticleSystem(xpositions = x, cutoff = cutoff, output = 0.0, nbatches = (1, 1))
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive

end

@testitem "ParticleSystem updating lists" setup = [Testing] begin

    using Test
    using StaticArrays
    using CellListMap

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    box = CellListMap.Box(sides, cutoff)

    # Initialize auxiliary linked lists
    system = ParticleSystem(
        xpositions = x,
        cutoff = cutoff,
        unitcell = sides,
        output = 0.0,
    )

    # Function to be evaluated for each pair: sum of displacements on x
    f(x, y, avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, x, box)
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ naive

    # Orthorhombic cell
    new_x = copy(x) .+ [rand(SVector{3, Float64}) for _ in 1:N]
    new_sides = sides + rand(SVector{3, Float64})
    new_cutoff = cutoff + rand()
    new_box = CellListMap.Box(new_sides, new_cutoff)
    new_naive = map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, new_x, new_box)
    # Update system
    system.xpositions .= new_x
    update_unitcell!(system, new_sides)
    update_cutoff!(system, new_cutoff)
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_naive
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_naive

    # If the number of particles and box change
    new_x, new_box = xatomic(10^5)
    new_cl = CellListMap.CellList(new_x, new_box)
    new_val = CellListMap.CellListMap._pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), 0.0, new_box, new_cl)
    resize!(system.xpositions, length(new_x))
    system.xpositions .= new_x
    update_unitcell!(system, [new_box.input_unit_cell.matrix[i, i] for i in 1:3])
    update_cutoff!(system, new_box.cutoff)
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val

    #
    # Triclinic cell
    #
    unitcell = [250 0 10; 10 250 0; 0 0 250]
    new_box = CellListMap.Box(unitcell, cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.CellListMap._pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), 0.0, new_box, new_cl, parallel = true)
    system = ParticleSystem(
        xpositions = system.xpositions,
        unitcell = unitcell,
        cutoff = cutoff,
        output = 0.0,
    )
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val

    # If the number of particles and box change
    cutoff = cutoff + rand()
    new_x, new_box = xatomic(10^4)
    unitcell = [250 0 10; 10 250 0; 0 0 250]
    new_box = CellListMap.Box(unitcell, cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.CellListMap._pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), 0.0, new_box, new_cl)
    resize!(system.xpositions, length(new_x))
    system.xpositions .= new_x
    update_unitcell!(system, unitcell)
    update_cutoff!(system, cutoff)
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val

    #
    # Non-periodic system
    #
    new_box = CellListMap.Box(CellListMap.limits(new_x), cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.CellListMap._pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), 0.0, new_box, new_cl, parallel = true)
    system = ParticleSystem(
        xpositions = system.xpositions,
        cutoff = cutoff,
        output = 0.0,
    )
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val

    # If the number of particles and box change
    cutoff = cutoff + rand()
    new_x, new_box = xatomic(10^4)
    new_box = CellListMap.Box(CellListMap.limits(new_x), cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.CellListMap._pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), 0.0, new_box, new_cl)
    resize!(system.xpositions, length(new_x))
    system.xpositions .= new_x
    update_cutoff!(system, cutoff)
    system.parallel = false
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test pairwise!((pair, avg_dx) -> f(pair.x, pair.y, avg_dx), system) ≈ new_val
end

@testitem "ParticleSystem applications" setup = [Testing] begin

    using Test
    using StaticArrays
    using CellListMap

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    box = CellListMap.Box(sides, cutoff)
    cl = CellListMap.CellList(x, box)

    # Function to be evaluated for each pair: build distance histogram
    function build_histogram!(d2, hist)
        d = sqrt(d2)
        ibin = floor(Int, d) + 1
        hist[ibin] += 1
        return hist
    end
    naive = map_naive!((x, y, i, j, d2, hist) -> build_histogram!(d2, hist), zeros(Int, 10), x, box)
    system = ParticleSystem(xpositions = x, cutoff = cutoff, unitcell = sides, output = zeros(Int, 10))
    @test naive == pairwise!((pair, hist) -> build_histogram!(pair.d2, hist), system)

    # Function to be evaluated for each pair: gravitational potential
    function potential(i, j, d2, u, mass)
        d2 == 0.0 && return u
        d = sqrt(d2)
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end
    mass = rand(N)
    naive = CellListMap.CellListMap._pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), 0.0, box, cl)
    system = ParticleSystem(xpositions = x, cutoff = cutoff, unitcell = sides, output = 0.0)
    @test pairwise!((pair, u) -> potential(pair.i, pair.j, pair.d2, u, mass), system) ≈ naive

    # Check the functionality of computing a different function from the same coordinates
    # (after pairwise!, the updated flag is false, so cell lists are not recomputed)
    naive = CellListMap.CellListMap._pairwise!((pair, u) -> u += pair.d2, 0.0, box, cl)
    @test pairwise!((pair, u) -> u += pair.d2, system) ≈ naive

    # Function to be evaluated for each pair: gravitational force
    function calc_forces!(x, y, i, j, d2, mass, forces)
        d2 == 0.0 && return forces
        G = 9.8 * mass[i] * mass[j] / d2
        d = sqrt(d2)
        df = (G / d) * (x - y)
        forces[i] = forces[i] - df
        forces[j] = forces[j] + df
        return forces
    end
    forces = [zeros(SVector{3, Float64}) for i in 1:N]
    naive = CellListMap.CellListMap._pairwise!((pair, forces) -> calc_forces!(pair.x, pair.y, pair.i, pair.j, pair.d2, mass, forces), copy(forces), box, cl)
    system = ParticleSystem(xpositions = x, cutoff = cutoff, unitcell = sides, output = forces)
    @test pairwise!((pair, forces) -> calc_forces!(pair.x, pair.y, pair.i, pair.j, pair.d2, mass, forces), system) ≈ naive

end

@testitem "ParticleSystemPositions interface" begin
    using StaticArrays
    using CellListMap

    # Construction from vector of vectors
    x = [rand(3) for _ in 1:10]
    p = CellListMap.ParticleSystemPositions(x)
    @test p isa AbstractVector
    @test length(p) == 10
    @test size(p) == (10,)
    @test p.updated[] == true

    # Construction from matrix
    xmat = rand(3, 10)
    p2 = CellListMap.ParticleSystemPositions(xmat)
    @test length(p2) == 10
    @test p2.updated[] == true

    # getindex
    @test p[1] isa SVector{3, Float64}
    @test p2[1] ≈ SVector{3, Float64}(xmat[:, 1])

    # setindex! sets updated flag
    p.updated[] = false
    p[1] = SVector(1.0, 2.0, 3.0)
    @test p[1] == SVector(1.0, 2.0, 3.0)
    @test p.updated[] == true

    # iterate
    count = 0
    for _ in p
        count += 1
    end
    @test count == length(p)

    # keys
    @test keys(p) == 1:10

    # resize! sets updated flag and changes length
    p.updated[] = false
    resize!(p, 15)
    @test length(p) == 15
    @test p.updated[] == true

    # empty! sets updated flag and clears
    p.updated[] = false
    empty!(p)
    @test length(p) == 0
    @test p.updated[] == true

    # Broadcasting assignment (setindex!) sets updated flag
    p3 = CellListMap.ParticleSystemPositions([rand(3) for _ in 1:5])
    p3.updated[] = false
    p3 .= [SVector(0.0, 0.0, 0.0) for _ in 1:5]
    @test p3.updated[] == true
    @test all(p3[i] == SVector(0.0, 0.0, 0.0) for i in 1:5)

    # Integration with ParticleSystem: flag is reset after pairwise!
    x = rand(SVector{3, Float64}, 100)
    sys = ParticleSystem(xpositions = x, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    @test sys.xpositions.updated[] == true
    pairwise!((pair, r) -> r += pair.d2, sys)
    @test sys.xpositions.updated[] == false
    # Modifying a position triggers the flag
    sys.xpositions[1] = SVector(0.5, 0.5, 0.5)
    @test sys.xpositions.updated[] == true
    pairwise!((pair, r) -> r += pair.d2, sys)
    @test sys.xpositions.updated[] == false

    # Integration with ParticleSystem2: both flags reset
    x = rand(SVector{3, Float64}, 100)
    y = rand(SVector{3, Float64}, 50)
    sys2 = ParticleSystem(xpositions = x, ypositions = y, cutoff = 0.1, unitcell = [1, 1, 1], output = 0.0)
    @test sys2.xpositions.updated[] == true
    @test sys2.ypositions.updated[] == true
    pairwise!((pair, r) -> r += pair.d2, sys2)
    @test sys2.xpositions.updated[] == false
    @test sys2.ypositions.updated[] == false
    # Modifying only y triggers only y flag
    sys2.ypositions[1] = SVector(0.5, 0.5, 0.5)
    @test sys2.xpositions.updated[] == false
    @test sys2.ypositions.updated[] == true

    # update_unitcell! sets updated flag
    sys.xpositions.updated[] = false
    update_unitcell!(sys, [2, 2, 2])
    @test sys.xpositions.updated[] == true

    # update_cutoff! sets updated flag
    sys.xpositions.updated[] = false
    update_cutoff!(sys, 0.2)
    @test sys.xpositions.updated[] == true

    # update_cutoff! sets both flags for ParticleSystem2
    sys2.xpositions.updated[] = false
    sys2.ypositions.updated[] = false
    update_cutoff!(sys2, 0.2)
    @test sys2.xpositions.updated[] == true
    @test sys2.ypositions.updated[] == true

    # update_unitcell! sets both flags for ParticleSystem2
    sys2.xpositions.updated[] = false
    sys2.ypositions.updated[] = false
    update_unitcell!(sys2, [2, 2, 2])
    @test sys2.xpositions.updated[] == true
    @test sys2.ypositions.updated[] == true
end
