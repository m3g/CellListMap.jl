module BasicTestsForPeriodicSystems

using Test
using StaticArrays
using CellListMap
using CellListMap.PeriodicSystems

@testset "PeriodicSystems disjoint sets" begin

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    mass = rand(N)

    # Function to be evalulated for each pair: gravitational potential
    function potential(i, j, d2, u, mass)
        d = sqrt(d2)
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end

    # Some simple disjoint set properties
    system = PeriodicSystem(
        xpositions=x,
        ypositions=y,
        cutoff=cutoff,
        unitcell=sides,
        output=0.0,
    )
    naive = CellListMap.map_naive!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), 0.0, x, y, Box(sides, cutoff))
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), system) ≈ naive
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), system) ≈ naive

    # Check different lcell
    system = PeriodicSystem(
        xpositions=x,
        ypositions=y,
        cutoff=cutoff,
        unitcell=sides,
        output=0.0,
        lcell=3,
    )
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), system) ≈ naive
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), system) ≈ naive

    # Test updating of the data on disjoint sets works fine
    for arrays in [
        [rand(SVector{2,Float64}, 1000), rand(SVector{2,Float64}, 100)], # with static vectors
        [[rand(2) for _ in 1:1000], [rand(2) for _ in 1:100]], # with standard vectors
    ]
        x = arrays[1]
        y = arrays[2]
        system = PeriodicSystem(
            xpositions=x,
            ypositions=y,
            cutoff=0.1,
            output=0.0,
            unitcell=[1, 1],
        )
        box = Box([1, 1], 0.1)
        cl = CellList(x, y, box)
        r = CellListMap.map_pairwise!((x, y, i, j, d2, r) -> r += d2, 0.0, box, cl)
        @test r ≈ PeriodicSystems.map_pairwise!((x, y, i, j, d2, r) -> r += d2, system)
        x = rand(SVector{2,Float64}, 1100)
        cl = UpdateCellList!(x, y, box, cl)
        r = CellListMap.map_pairwise!((x, y, i, j, d2, r) -> r += d2, 0.0, box, cl)
        system.xpositions = x
        @test r ≈ PeriodicSystems.map_pairwise!((x, y, i, j, d2, r) -> r += d2, system)

        y = rand(SVector{2,Float64}, 200)
        cl = UpdateCellList!(x, y, box, cl)
        r = CellListMap.map_pairwise!((x, y, i, j, d2, r) -> r += d2, 0.0, box, cl)
        system.ypositions = y
        @test r ≈ PeriodicSystems.map_pairwise!((x, y, i, j, d2, r) -> r += d2, system)

    end

end

@testset "PeriodicSystems parallelization" begin

    if Threads.nthreads() == 1
        println("""

             WARNING: Ideally, run a multi-threaded test to check the parallel versions.

        """)
    end

    # Function to be evalulated for each pair: sum of displacements on x
    f(x, y, avg_dx) = avg_dx + abs(x[1] - y[1])

    x, y, sides, cutoff = CellListMap.pathological_coordinates(2000)
    box = Box(sides, cutoff)
    naive = CellListMap.map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, x, box)

    # Check if changing lcell breaks something
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, lcell=1)
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, lcell=3)
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, lcell=5)
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive

    # Test if changing the number of batches breaks anything
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(3, 5))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(1, 1))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(1, 7))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(7, 1))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(13, 17))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system = PeriodicSystem(xpositions=x, unitcell=sides, cutoff=cutoff, output=0.0, nbatches=(4, 16))
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive

end

@testset "PeriodicSystems updating lists" begin

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    box = Box(sides, cutoff)

    # Initialize auxiliary linked lists
    system = PeriodicSystem(
        xpositions=x,
        cutoff=cutoff,
        unitcell=sides,
        output=0.0,
    )

    # Function to be evalulated for each pair: sum of displacements on x
    f(x, y, avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = CellListMap.map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, x, box)
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ naive

    # Orthorhombic cell
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = Box(new_sides, new_cutoff)
    new_naive = CellListMap.map_naive!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, new_x, new_box)
    # Update system
    system.xpositions .= new_x
    update_unitcell!(system, new_sides)
    update_cutoff!(system, new_cutoff)
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_naive
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_naive

    # If the number of particles and box change
    new_x, new_box = CellListMap.xatomic(10^5)
    new_cl = CellList(new_x, new_box)
    new_val = CellListMap.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, new_box, new_cl)
    system.xpositions = new_x
    update_unitcell!(system, [new_box.unit_cell.matrix[i, i] for i in 1:3])
    update_cutoff!(system, new_box.cutoff)
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val

    #
    # Triclinic cell
    #
    unitcell = [250 0 10; 10 250 0; 0 0 250]
    new_box = Box(unitcell, cutoff)
    new_cl = UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, new_box, new_cl, parallel=true)
    system = PeriodicSystem(
        xpositions=system.xpositions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=0.0,
    )
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val

    # If the number of particles and box change
    cutoff = cutoff + rand()
    new_x, new_box = CellListMap.xatomic(10^4)
    unitcell = [250 0 10; 10 250 0; 0 0 250]
    new_box = Box(unitcell, cutoff)
    new_cl = UpdateCellList!(new_x, new_box, new_cl)
    new_val = CellListMap.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), 0.0, new_box, new_cl)
    system.xpositions = new_x
    update_unitcell!(system, unitcell)
    update_cutoff!(system, cutoff)
    system.parallel = false
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val
    system.parallel = true
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, avg_dx) -> f(x, y, avg_dx), system) ≈ new_val

end

@testset "ParticleSystems applications" begin

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    box = Box(sides, cutoff)
    cl = CellList(x, box)

    # Function to be evalulated for each pair: build distance histogram
    function build_histogram!(d2, hist)
        d = sqrt(d2)
        ibin = floor(Int, d) + 1
        hist[ibin] += 1
        return hist
    end
    naive = CellListMap.map_naive!((x, y, i, j, d2, hist) -> build_histogram!(d2, hist), zeros(Int, 10), x, box)
    system = PeriodicSystem(xpositions=x, cutoff=cutoff, unitcell=sides, output=zeros(Int, 10))
    @test naive == PeriodicSystems.map_pairwise!((x, y, i, j, d2, hist) -> build_histogram!(d2, hist), system)

    # Function to be evalulated for each pair: gravitational potential
    function potential(i, j, d2, u, mass)
        d = sqrt(d2)
        u = u - 9.8 * mass[i] * mass[j] / d
        return u
    end
    mass = rand(N)
    naive = CellListMap.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), 0.0, box, cl)
    system = PeriodicSystem(xpositions=x, cutoff=cutoff, unitcell=sides, output=0.0)
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, u) -> potential(i, j, d2, u, mass), system) ≈ naive

    # Function to be evalulated for each pair: gravitational force
    function calc_forces!(x, y, i, j, d2, mass, forces)
        G = 9.8 * mass[i] * mass[j] / d2
        d = sqrt(d2)
        df = (G / d) * (x - y)
        forces[i] = forces[i] - df
        forces[j] = forces[j] + df
        return forces
    end
    forces = [zeros(SVector{3,Float64}) for i in 1:N]
    naive = CellListMap.map_pairwise!((x, y, i, j, d2, forces) -> calc_forces!(x, y, i, j, d2, mass, forces), copy(forces), box, cl)
    system = PeriodicSystem(xpositions=x, cutoff=cutoff, unitcell=sides, output=forces)
    @test PeriodicSystems.map_pairwise!((x, y, i, j, d2, forces) -> calc_forces!(x, y, i, j, d2, mass, forces), system) ≈ naive

end # testset

end # module
