using CellListMap
using StaticArrays
using Test

@testset "CellListMap.jl" begin

    if Threads.nthreads() == 1
        println("""

             WARNING: Ideally, run a multi-threaded test to check the parallel versions.

        """)
    end

    # Number of particles, sides and cutoff
    sides = @SVector [250,250,250]
    cutoff = 10.
    box = Box(sides,cutoff)

    # Particle positions
    N = 2000
    x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    # Add some pathological coordinates
    x[1] = -sides/2
    x[10] = sides
    x[100] = @SVector [sides[1]/2, -sides[2]/2, 2*sides[3]]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]

    # Initialize auxiliary linked lists
    cl = CellList(x,box)

    # Function to be evalulated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=false) ≈ naive

    # Check if changing lcell breaks something
    box = Box(sides,cutoff,lcell=1); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=2); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=3); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=5); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive

    # Test random triclinic boxes
    for N in 2:3, 
        M in rand(10:20), 
        UnitCellType in [ TriclinicCell, OrthorhombicCell],
        parallel in [ false, true ],
        lcell in 1:3
        @test CellListMap.check_random_cells(
            N,M,
            UnitCellType=UnitCellType,
            parallel=parallel,
            lcell=lcell,
            show_progress=false
        )[1] 
    end

    # Function to be evalulated for each pair: build distance histogram
    function build_histogram!(d2,hist)
        d = sqrt(d2)
        ibin = floor(Int,d) + 1
        hist[ibin] += 1
        return hist
    end

    naive = CellListMap.map_naive!(
        (x,y,i,j,d2,hist) -> build_histogram!(d2,hist),
        zeros(Int,10),x,box
    )
    @test map_pairwise!(
         (x,y,i,j,d2,hist) -> build_histogram!(d2,hist),
         zeros(Int,10),box,cl,parallel=true
    ) ≈ naive
    @test map_pairwise!(
         (x,y,i,j,d2,hist) -> build_histogram!(d2,hist),
         zeros(Int,10),box,cl,parallel=false
    ) ≈ naive

    # Function to be evalulated for each pair: build distance histogram
    function potential(i,j,d2,u,mass)
        d = sqrt(d2)
        u = u - 9.8*mass[i]*mass[j]/d
        return u
    end

    # Run pairwise computation
    mass = rand(N)
    naive = CellListMap.map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,box)
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,box,cl,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,box,cl,parallel=false) ≈ naive

    # Function to be evalulated for each pair: build distance histogram
    function calc_forces!(x,y,i,j,d2,mass,forces)
        G = 9.8*mass[i]*mass[j]/d2
        d = sqrt(d2)
        df = (G/d)*(x - y)
        forces[i] = forces[i] - df
        forces[j] = forces[j] + df
        return forces
    end

    # forces
    forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

    # Run pairwise computation
    naive = CellListMap.map_naive!(
        (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
        copy(forces),x,box
    )
    @test map_pairwise!(
        (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
        copy(forces),box,cl,parallel=true
    ) ≈ naive
    @test map_pairwise!(
        (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
        copy(forces),box,cl,parallel=false
    ) ≈ naive

    # Compute some properties of disjoint sets 
    box = Box(sides,cutoff,lcell=1)
    cl = CellList(x,y,box)

    naive = CellListMap.map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,y,box)
    @test map_pairwise!(
        (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),
        0.0,box,cl,parallel=false
    ) ≈ naive
    @test map_pairwise!(
        (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),
        0.0,box,cl,parallel=true
    ) ≈ naive

    # Check different lcell
    box = Box(sides,cutoff,lcell=3)
    cl = CellList(x,y,box)
    @test map_pairwise!(
        (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),
        0.0,box,cl,parallel=false
    ) ≈ naive
    @test map_pairwise!(
        (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),
        0.0,box,cl,parallel=true
    ) ≈ naive

    # Test the examples, to check further if the parallelization didn't break something
    N = 100_000
    x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    @test CellListMap.test1(parallel=true,x=x)[2] ≈ CellListMap.test1(parallel=false,x=x)[2]
    @test CellListMap.test2(parallel=true,x=x) ≈ CellListMap.test2(parallel=false,x=x)
    @test CellListMap.test3(parallel=true,x=x) ≈ CellListMap.test3(parallel=false,x=x)
    @test CellListMap.test4(parallel=true,x=x) ≈ CellListMap.test4(parallel=false,x=x)
    @test count(CellListMap.test5(parallel=true,x=x,y=y) .≈ CellListMap.test5(parallel=false,x=x,y=y)) == 3

    pairs1 = sort!(CellListMap.test7(parallel=true,x=x),by=x->x[3])
    pairs2 = sort!(CellListMap.test7(parallel=false,x=x),by=x->x[3])
    @test count([ count(pairs1[i] .≈ pairs2[i]) == 3 for i in 1:length(pairs1) ]) == length(pairs1)

    x = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500 ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500_000 ]
    @test count(CellListMap.test6(parallel=true,x=x,y=y) .≈ CellListMap.test6(parallel=false,x=x,y=y)) == 3

    # invert x and y to test swap
    ixy = CellListMap.test6(parallel=false,x=x,y=y) 
    iyx = CellListMap.test6(parallel=false,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 
    ixy = CellListMap.test6(parallel=true,x=x,y=y) 
    iyx = CellListMap.test6(parallel=true,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 

    # Test some fractional box lengths with the packmol test
    @test CellListMap.packmol_test(parallel=false,
        sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=TriclinicCell)[1]
    @test CellListMap.packmol_test(parallel=true,
        sides=[18.4,30.1,44], tol=2, UnitCellType=TriclinicCell)[1]
    @test CellListMap.packmol_test(parallel=false,
        sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=OrthorhombicCell)[1]
    @test CellListMap.packmol_test(parallel=true,
        sides=[18.4,30.1,44], tol=2, UnitCellType=OrthorhombicCell)[1]

#    # Test resizing of the cell lists
#    x = [ rand(SVector{3,Float64}) for i in 1:1000 ]
#    box = Box([0.83,0.41,0.97],0.1)
#    cl = CellList(x,box) 
#    @test length(cl.cwp) == 924
#
#    box = Box([0.33,0.41,0.97],0.1)
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 924 
#
#    box = Box([0.83,0.81,0.97],0.1)
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 1598
#
#    x .= 0.9*x
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 1598
#
#    x .= 1.2*x
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 1598
#
#    box = Box([0.83,0.81,0.97],0.2)
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 1598
#
#    box = Box([0.83,0.81,0.97],0.05)
#    cl = UpdateCellList!(x,box,cl)   
#    @test length(cl.cwp) == 8737

end

include("./namd/compare_with_namd.jl")
