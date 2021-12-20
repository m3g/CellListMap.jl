using CellListMap
using StaticArrays
using Test

# Loads Unitful and ForwardDiff
include("../src/examples/generic_types.jl")

@testset "CellListMap.jl" begin

    if Threads.nthreads() == 1
        println("""

             WARNING: Ideally, run a multi-threaded test to check the parallel versions.

        """)
    end

    # Number of particles, sides and cutoff
    sides = @SVector [250.,250.,250.]
    cutoff = 10.
    box = Box(sides,cutoff)

    # Particle positions
    N = 2000
    x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    # Add some pathological coordinates
    x[1] = -sides/2
    x[2] = -sides/2 + @SVector [ nextfloat(0.), nextfloat(0.), sides[3]*rand() ]
    x[3] = -sides/2 + @SVector [ prevfloat(0.), prevfloat(0.), sides[3]*rand() ]
    x[4] = sides/2 + @SVector [ nextfloat(0.), nextfloat(0.), sides[3]*rand() ]
    x[5] = sides/2 + @SVector [ prevfloat(0.), prevfloat(0.), sides[3]*rand() ]
    x[10] = sides
    x[11] = sides + @SVector [ nextfloat(0.), nextfloat(0.), sides[3]*rand() ]
    x[12] = sides + @SVector [ prevfloat(0.), prevfloat(0.), sides[3]*rand() ]
    x[13] = @SVector [ nextfloat(0.), nextfloat(0.), sides[3]*rand() ]
    x[14] = @SVector [ prevfloat(0.), prevfloat(0.), sides[3]*rand() ]
    x[100] = @SVector [sides[1]/2, -sides[2]/2, 2*sides[3]]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]

    # Initialize auxiliary linked lists
    cl = CellList(x,box)

    # Function to be evalulated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=false) ≈ naive

    # Test updating the cell lists
    new_cl = deepcopy(cl)
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N ]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = Box(sides,cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    new_box = Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive

    # Test the input as a matrix
    xmat = zeros(3,N) 
    for i in 1:N
        for j in 1:3
            xmat[j,i] = x[i][j] 
        end
    end
    cl_mat = CellList(xmat,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl_mat,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl_mat,parallel=false) ≈ naive

    # Check if changing lcell breaks something
    box = Box(sides,cutoff,lcell=1); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=2); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=3); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = Box(sides,cutoff,lcell=5); cl = CellList(x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive

    # Test if changing the number of batches breaks anything
    cl = CellList(x,box,nbatches=(3,5))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellList(x,box,nbatches=(1,1))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellList(x,box,nbatches=(1,7))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellList(x,box,nbatches=(7,1))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellList(x,box,nbatches=(13,17))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellList(x,box,nbatches=(4,16))
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive

    # Test random cells of all possible types
    for N in 2:3, 
        M in rand(10:20), 
        UnitCellType in [ TriclinicCell, OrthorhombicCell ],
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

    naive = CellListMap.map_naive!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10),x,box)
    @test map_pairwise!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10),box,cl,parallel=true) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10),box,cl,parallel=false) ≈ naive

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
    naive = CellListMap.map_naive!( (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces), copy(forces),x,box)
    @test map_pairwise!( (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces), copy(forces),box,cl,parallel=true) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces), copy(forces),box,cl,parallel=false) ≈ naive

    #
    # Compute some properties of disjoint sets 
    #
    box = Box(sides,cutoff,lcell=1)
    cl = CellList(x,y,box)

    naive = CellListMap.map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,y,box)
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=true) ≈ naive

    # Test disjoint sets, with matrices
    xmat = zeros(3,length(x))
    for i in 1:length(x)
        for j in 1:3
            xmat[j,i] = x[i][j]
        end
    end
    ymat = zeros(3,length(y))
    for i in 1:length(y)
        for j in 1:3
            ymat[j,i] = y[i][j]
        end
    end
    cl_mat = CellList(xmat,ymat,box)
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0, box, cl_mat, parallel=false) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0, box, cl_mat, parallel=true) ≈ naive

    # Check different lcell
    box = Box(sides,cutoff,lcell=3)
    cl = CellList(x,y,box)
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=true) ≈ naive

    # Test if changing the number of batches breaks anything
    cl = CellList(x,y,box,nbatches=(1,1))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellList(x,y,box,nbatches=(3,5))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellList(x,y,box,nbatches=(7,1))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellList(x,y,box,nbatches=(1,7))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellList(x,y,box,nbatches=(4,16))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellList(x,y,box,nbatches=(13,17))
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive

    # Test the examples, to check further if the parallelization didn't break something
    N = 100_000
    x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    @test CellListMap.Examples.average_displacement(parallel=true) ≈ CellListMap.Examples.average_displacement(parallel=false)[1]
    @test CellListMap.Examples.distance_histogram(parallel=true,x=x) ≈ CellListMap.Examples.distance_histogram(parallel=false,x=x)
    @test CellListMap.Examples.gravitational_potential(parallel=true,x=x) ≈ CellListMap.Examples.gravitational_potential(parallel=false,x=x)
    @test CellListMap.Examples.gravitational_force(parallel=true,x=x) ≈ CellListMap.Examples.gravitational_force(parallel=false,x=x)
    @test count(CellListMap.Examples.nearest_neighbor(parallel=true,x=x,y=y) .≈ CellListMap.Examples.nearest_neighbor(parallel=false,x=x,y=y)) == 3
    
    function pair_match(p1,p2) 
        p1[3] ≈ p2[3] || return false 
        p1[1] == p2[1] && p1[2] == p2[2] && return true
        p1[1] == p2[2] && p1[2] == p2[1] && return true
    end
    pairs1 = sort!(CellListMap.Examples.neighborlist(parallel=true,x=x),by=x->x[3])
    pairs2 = sort!(CellListMap.Examples.neighborlist(parallel=false,x=x),by=x->x[3])
    @test length(pairs1) == length(pairs2)
    @test count(pair_match.(pairs1,pairs2)) == length(pairs1)

    x = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500 ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500_000 ]
    @test count(CellListMap.Examples.nearest_neighbor_nopbc(parallel=true,x=x,y=y) .≈ CellListMap.Examples.nearest_neighbor_nopbc(parallel=false,x=x,y=y)) == 3

    # invert x and y to test swap
    ixy = CellListMap.Examples.nearest_neighbor_nopbc(parallel=false,x=x,y=y) 
    iyx = CellListMap.Examples.nearest_neighbor_nopbc(parallel=false,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 
    ixy = CellListMap.Examples.nearest_neighbor_nopbc(parallel=true,x=x,y=y) 
    iyx = CellListMap.Examples.nearest_neighbor_nopbc(parallel=true,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 

    # Test some fractional box lengths with the packmol test
    @test CellListMap.Examples.packmol(parallel=false, sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=TriclinicCell)[1]
    @test CellListMap.Examples.packmol(parallel=true, sides=[18.4,30.1,44], tol=2, UnitCellType=TriclinicCell)[1]
    @test CellListMap.Examples.packmol(parallel=false, sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=OrthorhombicCell)[1]
    @test CellListMap.Examples.packmol(parallel=true, sides=[18.4,30.1,44], tol=2, UnitCellType=OrthorhombicCell)[1]
    
    # Testing the propagation of types in automatic differentiation
    @test generic_types(false) == (true,u"nm^2",Measurement{Float64})

    # Test when we have pathologically few number of particles
    x = [ Float64[1,1,1] ]
    y = [ Float64[1.05,1,1], Float64[0,0,0]  ]
    @test CellListMap.neighborlist(x,y,0.1)[1] == (1, 1, 0.050000000000000044)
    z = [ Float64[1,1,1], Float64[1.05,1,1], Float64[0,0,0]  ]
    @test CellListMap.neighborlist(z,0.1)[1] == (1, 2, 0.050000000000000044)

end

include("./namd/compare_with_namd.jl")
include("./neighborlist.jl")
