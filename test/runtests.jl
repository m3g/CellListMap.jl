using TestItemRunner
using TestItems

@testmodule AllocTest begin
    # This module defines the Allocs struct and the comparison operators
    # to conditionally compare the number of allocations based on the
    # BUILD_IS_PRODUCTION_BUILD environment variable.
    export Allocs
    @kwdef struct Allocs
        prodbuild::Bool = haskey(ENV, "BUILD_IS_PRODUCTION_BUILD") && ENV["BUILD_IS_PRODUCTION_BUILD"] == "true"
        allocs::Int
    end
    Allocs(allocs::Int) = Allocs(; allocs)
    import Base: ==, >, <
    ==(a::Int, b::Allocs) = b.prodbuild ? a == b.allocs : true
    <(a::Int, b::Allocs) = b.prodbuild ? a < b.allocs : true
    ==(a::Allocs, b::Int) = a.prodbuild ? a.allocs == b : true
    <(a::Allocs, b::Int) = a.prodbuild ? a.allocs < b : true
end

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(CellListMap)
end

@testitem "disjoint sets" begin

    using CellListMap
    using StaticArrays

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    mass = rand(N)

    # Function to be evalulated for each pair: gravitational potential
    function potential(i,j,d2,u,mass)
        d = sqrt(d2)
        u = u - 9.8*mass[i]*mass[j]/d
        return u
    end

    # Some simple disjoint set properties
    box = Box(sides,cutoff,lcell=1)
    cl = CellList(x,y,box)

    naive = CellListMap.map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,y,box)
    @test map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,u) -> potential(i,j,d2,u,mass), 0.0,box,cl,parallel=true) ≈ naive

    # Test disjoint sets, with matrices
    xmat = zeros(3,length(x))
    for i in eachindex(x)
        for j in 1:3
            xmat[j,i] = x[i][j]
        end
    end
    ymat = zeros(3,length(y))
    for i in eachindex(x)
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

    # Test updating of the data on disjoint sets works fine
    for arrays in [ 
       [ rand(SVector{2,Float64},1000), rand(SVector{2,Float64},100) ], # with static vectors
       [ [ rand(2) for _ in 1:1000 ], [ rand(2) for _ in 1:100 ] ], # with standard vectors
    ]
        local x = arrays[1]
        local y = arrays[2]
        local box = Box([1,1],0.1)
        local cl = CellList(x,y,box)
        r_naive = CellListMap.map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = map_pairwise!((x,y,i,j,d2,r) -> r += d2, 0., box, cl)
        @test r_naive ≈ r
        x = rand(SVector{2,Float64},1100)
        cl = UpdateCellList!(x, y, box, cl)
        r_naive = CellListMap.map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = map_pairwise!((x,y,i,j,d2,r) -> r += d2, 0., box, cl)
        @test r_naive ≈ r
        y = rand(SVector{2,Float64},200)
        cl = UpdateCellList!(x, y, box, cl)
        r_naive = CellListMap.map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = map_pairwise!((x,y,i,j,d2,r) -> r += d2, 0., box, cl)
        @test r_naive ≈ r
    end

end

@testitem "matrix inputs" begin

    using CellListMap
    using StaticArrays

    # Function to be evalulated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    # Number of particles, sides and cutoff
    N = 2000
    sides = @SVector [250.,250.,250.]
    cutoff = 10.
    box = Box(sides,cutoff)

    # Test the input as a matrix
    x = rand(SVector{3,Float64}, N)
    xmat = zeros(3,N) 
    for i in 1:N
        for j in 1:3
            xmat[j,i] = x[i][j] 
        end
    end
    cl_mat = CellList(xmat,box)
    naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl_mat,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl_mat,parallel=false) ≈ naive

end

@testitem "parallelization" begin

    using CellListMap
    using StaticArrays

    if Threads.nthreads() == 1
        println("""

             WARNING: Ideally, run a multi-threaded test to check the parallel versions.

        """)
    end

    # Function to be evalulated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    x, y, sides, cutoff = CellListMap.pathological_coordinates(2000)
    box = Box(sides, cutoff)
    naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)

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

end

@testitem "updating lists" begin

    using CellListMap
    using StaticArrays

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    box = Box(sides, cutoff)

    # Initialize auxiliary linked lists
    cl = CellList(x,box)

    # Function to be evalulated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.0,x,box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=false) ≈ naive

    # Orthorhombic cell
    new_cl = deepcopy(cl)
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N ]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = Box(new_sides,new_cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_x, new_box = CellListMap.xatomic(10^4)
    new_cl = CellList(new_x,new_box)
    new_x, new_box = CellListMap.xatomic(10^4 + 10^3)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box) # slow
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    #
    # Using auxiliary preallocated arrays
    #
    new_cl = deepcopy(cl)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N ]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = Box(new_sides,new_cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change (more)
    new_x, new_box = CellListMap.xatomic(10^4)
    new_cl = CellList(new_x,new_box)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x, new_box = CellListMap.xatomic(10^4 + 10^3)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl, new_aux)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change (less)
    new_x, new_box = CellListMap.xatomic(10^4)
    new_cl = CellList(new_x,new_box)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x, new_box = CellListMap.xatomic(10^4 - 10^3)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl, new_aux)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    #
    # Triclinic cell
    #
    new_box = Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_box = Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_x, new_box = CellListMap.xatomic(10^4)
    new_box = Box([ 200   0  10 
                     15 200   0 
                      0   0 200 ],cutoff)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # Using auxiliary preallocated aux arrays
    new_box = Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_box = Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_x, new_box = CellListMap.xatomic(10^4)
    new_box = Box([ 200   0  10 
                     15 200   0 
                      0   0 200 ],cutoff)
    new_naive = CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

end



@testitem "applications" begin

    using CellListMap
    using StaticArrays

    # Loads Unitful and ForwardDiff
    dir = @__DIR__
    include(dir*"/../src/examples/generic_types.jl")

    N = 2000
    x, y, sides, cutoff = CellListMap.pathological_coordinates(N)
    box = Box(sides, cutoff)
    cl = CellList(x,box)

    # Function to be evalulated for each pair: build distance histogram
    function build_histogram!(d2,hist)
        d = sqrt(d2)
        ibin = floor(Int,d) + 1
        hist[ibin] += 1
        return hist
    end

    naive = CellListMap.map_naive!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10), x, box)
    @test map_pairwise!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10), box, cl, parallel=true) ≈ naive
    @test map_pairwise!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10), box, cl, parallel=false) ≈ naive

    # Function to be evalulated for each pair: gravitational potential
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

    # Function to be evalulated for each pair: gravitational force
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
    @test generic_types_triclinic(false) == (true,u"nm^2",Measurement{Float64})

    # Test when we have pathologically few number of particles
    x = [ Float64[1,1,1] ]
    y = [ Float64[1.05,1,1], Float64[0,0,0]  ]
    @test all(CellListMap.neighborlist(x,y,0.1)[1] .≈ (1, 1, 0.05))
    z = [ Float64[1,1,1], Float64[1.05,1,1], Float64[0,0,0]  ]
    @test all(CellListMap.neighborlist(z,0.1)[1] .≈ (1, 2, 0.05))
    x = SVector{3,Float64}[]
    @test CellListMap.neighborlist(x,0.1,unitcell=[1,1,1]) == Tuple{Int64, Int64, Float64}[]
    x = Vector{Float64}[]
    @test CellListMap.neighborlist(x,0.1,unitcell=[1,1,1]) == Tuple{Int64, Int64, Float64}[]

end

@testitem "pathological cells" begin
    using StaticArrays
    # This function is an interesting function because it sort of counts 
    # the indexes of the particles within the cutoff. However, we need to 
    # check for numerical precision innacuracies before adding them, otherwise
    # pairs that match in one test can be skiped in another. This kind of issue
    # exists for any property that does not goes to zero when the distance
    # approaches that of the cutoff.
    function f(i,j,d2,out)
        d = sqrt(d2)
        if !(d ≈ 0.2) 
            out += (i + j)
        end
        return out
    end
    l = sqrt(2)/2
    for m in [
        @SMatrix[1 0; 0 1],
        @SMatrix[l 0; l 1],
        @SMatrix[1.1 0; 0 1],
        @SMatrix[1.2 0; 0 1],
        @SMatrix[1 0; 0 1.1],
        @SMatrix[1 0; 0 1.2],
        @SMatrix[1 0.2; 0 1.2],
        @SMatrix[1 0.2; 0.2 1.2],
        @SMatrix[1.2 0.2; 0.2 1.2],
    ]
        box = Box(m,0.2)
        for x in [
            100 .* rand(SVector{2,Float64}, 100),
            [SVector{2,Float64}(0.1 * i + 0.1*j, 0.1 * j + 0.1*j) for i in 0:5 for j in 0:5],
            [SVector{2,Float64}(0.1 * i, 0.1 * j) for i in 0:5 for j in 0:5]
        ]
            cl = CellList(x, box)
            d0 = map_pairwise((x, y, i, j, d2, out) -> f(i,j,d2,out), 0, box, cl)
            d1 = CellListMap.map_naive!((x, y, i, j, d2, out) -> f(i,j,d2,out), 0, x, box)
            @test d0 ≈ d1
        end
    end

    # sph test (https://github.com/m3g/CellListMap.jl/issues/95)
    using StaticArrays: SVector
    using CellListMap: neighborlist
    p2d = SVector{2, Float64}[[0.0, 2.52], [0.02, 2.56], [3.98, 2.96], [4.0, 0.26], [4.0, 2.5]]
    r = 0.06788225099390856
    @test length(neighborlist(p2d, r)) == 1

end

include("$(@__DIR__)/namd/compare_with_namd.jl")
include("$(@__DIR__)/gromacs/compare_with_gromacs.jl")
include("$(@__DIR__)/BasicForParticleSystem.jl")
include("$(@__DIR__)/namd/ParticleSystem_vs_NAMD.jl")

@run_package_tests
