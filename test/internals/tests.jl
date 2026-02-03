
@testitem "disjoint sets" setup=[Testing] begin

    using CellListMap
    using StaticArrays

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    mass = rand(N)

    # Function to be evaluated for each pair: gravitational potential
    function potential(i,j,d2,u,mass)
        d = sqrt(d2)
        u = u - 9.8*mass[i]*mass[j]/d
        return u
    end

    # Some simple disjoint set properties
    box = CellListMap.Box(sides,cutoff,lcell=1)
    cl = CellListMap.CellList(x,y,box)

    naive = map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,y,box)
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=true) ≈ naive

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
    cl_mat = CellListMap.CellList(xmat,ymat,box)
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0, box, cl_mat, parallel=false) ≈ naive
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0, box, cl_mat, parallel=true) ≈ naive

    # Check different lcell
    box = CellListMap.Box(sides,cutoff,lcell=3)
    cl = CellListMap.CellList(x,y,box)
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=true) ≈ naive

    # Test if changing the number of batches breaks anything
    cl = CellListMap.CellList(x,y,box,nbatches=(1,1))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellListMap.CellList(x,y,box,nbatches=(3,5))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellListMap.CellList(x,y,box,nbatches=(7,1))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellListMap.CellList(x,y,box,nbatches=(1,7))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellListMap.CellList(x,y,box,nbatches=(4,16))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive
    cl = CellListMap.CellList(x,y,box,nbatches=(13,17))
    @test pairwise!((pair,u) -> potential(pair.i,pair.j,pair.d2,u,mass), 0.0,box,cl,parallel=false) ≈ naive

    # Test updating of the data on disjoint sets works fine
    for arrays in [ 
       [ rand(SVector{2,Float64},1000), rand(SVector{2,Float64},100) ], # with static vectors
       [ [ rand(2) for _ in 1:1000 ], [ rand(2) for _ in 1:100 ] ], # with standard vectors
    ]
        local x = arrays[1]
        local y = arrays[2]
        local box = CellListMap.Box([1,1],0.1)
        local cl = CellListMap.CellList(x,y,box)
        r_naive = map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = pairwise!((pair,r) -> r += pair.d2, 0., box, cl)
        @test r_naive ≈ r
        x = rand(SVector{2,Float64},1100)
        cl = CellListMap.UpdateCellList!(x, y, box, cl)
        r_naive = map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = pairwise!((pair,r) -> r += pair.d2, 0., box, cl)
        @test r_naive ≈ r
        y = rand(SVector{2,Float64},200)
        cl = CellListMap.UpdateCellList!(x, y, box, cl)
        r_naive = map_naive!((x,y,i,j,d2,r) -> r += d2, 0., x, y, box)
        r = pairwise!((pair,r) -> r += pair.d2, 0., box, cl)
        @test r_naive ≈ r
    end

end

@testitem "matrix inputs" setup=[Testing] begin

    using CellListMap
    using StaticArrays

    # Function to be evaluated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    # Number of particles, sides and cutoff
    N = 2000
    sides = @SVector [250.,250.,250.]
    cutoff = 10.
    box = CellListMap.Box(sides,cutoff)

    # Test the input as a matrix
    x = rand(SVector{3,Float64}, N)
    xmat = zeros(3,N) 
    for i in 1:N
        for j in 1:3
            xmat[j,i] = x[i][j] 
        end
    end
    cl_mat = CellListMap.CellList(xmat,box)
    naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl_mat,parallel=true) ≈ naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl_mat,parallel=false) ≈ naive

end

@testitem "parallelization" setup=[Testing] begin

    using CellListMap
    using StaticArrays

    if Threads.nthreads() == 1
        println("""

             WARNING: Ideally, run a multi-threaded test to check the parallel versions.

        """)
    end

    # Function to be evaluated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    x, y, sides, cutoff = pathological_coordinates(2000)
    box = CellListMap.Box(sides, cutoff)
    naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box)

    # Check if changing lcell breaks something
    box = CellListMap.Box(sides,cutoff,lcell=1); cl = CellListMap.CellList(x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = CellListMap.Box(sides,cutoff,lcell=2); cl = CellListMap.CellList(x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = CellListMap.Box(sides,cutoff,lcell=3); cl = CellListMap.CellList(x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    box = CellListMap.Box(sides,cutoff,lcell=5); cl = CellListMap.CellList(x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive

    # Test if changing the number of batches breaks anything
    cl = CellListMap.CellList(x,box,nbatches=(3,5))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellListMap.CellList(x,box,nbatches=(1,1))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellListMap.CellList(x,box,nbatches=(1,7))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellListMap.CellList(x,box,nbatches=(7,1))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellListMap.CellList(x,box,nbatches=(13,17))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    cl = CellListMap.CellList(x,box,nbatches=(4,16))
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive

end

@testitem "updating lists" setup=[Testing] begin

    using CellListMap
    using StaticArrays

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    box = CellListMap.Box(sides, cutoff)

    # Initialize auxiliary linked lists
    cl = CellListMap.CellList(x,box)

    # Function to be evaluated for each pair: sum of displacements on x
    f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

    naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.0,x,box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=true) ≈ naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,box,cl,parallel=false) ≈ naive

    # Orthorhombic cell
    new_cl = deepcopy(cl)
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N ]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = CellListMap.Box(new_sides,new_cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_x, new_box = xatomic(10^4)
    new_cl = CellListMap.CellList(new_x,new_box)
    new_x, new_box = xatomic(10^4 + 10^3)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box) # slow
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    #
    # Using auxiliary preallocated arrays
    #
    new_cl = deepcopy(cl)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x = copy(x) .+ [rand(SVector{3,Float64}) for _ in 1:N ]
    new_sides = sides + rand(SVector{3,Float64})
    new_cutoff = cutoff + rand()
    new_box = CellListMap.Box(new_sides,new_cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change (more)
    new_x, new_box = xatomic(10^4)
    new_cl = CellListMap.CellList(new_x,new_box)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x, new_box = xatomic(10^4 + 10^3)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl, new_aux)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change (less)
    new_x, new_box = xatomic(10^4)
    new_cl = CellListMap.CellList(new_x,new_box)
    new_aux = CellListMap.AuxThreaded(new_cl)
    new_x, new_box = xatomic(10^4 - 10^3)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x, new_box, new_cl, new_aux)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    #
    # Triclinic cell
    #
    new_box = CellListMap.Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_box = CellListMap.Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    new_x, new_box = xatomic(10^4)
    new_box = CellListMap.Box([ 200   0  10 
                     15 200   0 
                      0   0 200 ],cutoff)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # Using auxiliary preallocated aux arrays
    new_box = CellListMap.Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive
    # If the number of particles and box change
    new_box = CellListMap.Box([ 250   0  10 
                     10 250   0 
                      0   0 250 ],cutoff)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    new_x, new_box = xatomic(10^4)
    new_box = CellListMap.Box([ 200   0  10 
                     15 200   0 
                      0   0 200 ],cutoff)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    # Same as above, but with parallel=false on the update
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux; parallel=false)
    new_x, new_box = xatomic(10^4)
    new_box = CellListMap.Box([ 200   0  10 
                     15 200   0 
                      0   0 200 ],cutoff)
    new_naive = map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,new_x,new_box)
    new_cl = CellListMap.UpdateCellList!(new_x,new_box,new_cl,new_aux; parallel=false)
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=false) ≈ new_naive
    @test pairwise!((pair,avg_dx) -> f(pair.x,pair.y,avg_dx),0.,new_box,new_cl,parallel=true) ≈ new_naive

    # Internal argument-error test: the should never reach this test 
    x = rand(SVector{3,Float64},100)
    cl1 = CellListMap.CellList(x,CellListMap.Box([1,1,1],0.1))
    cl2 = CellListMap.CellList(x,CellListMap.Box([1.2,1.2,1.2],0.1))
    @test_throws ArgumentError CellListMap.merge_cell_lists!(cl1,cl2)

end

@testitem "applications" setup=[Testing, Examples] begin

    using CellListMap
    using StaticArrays
    using Unitful
    using Measurements

    N = 2000
    x, y, sides, cutoff = pathological_coordinates(N)
    box = CellListMap.Box(sides, cutoff)
    cl = CellListMap.CellList(x,box)

    # Function to be evaluated for each pair: build distance histogram
    function build_histogram!(d2,hist)
        d = sqrt(d2)
        ibin = floor(Int,d) + 1
        hist[ibin] += 1
        return hist
    end

    naive = map_naive!( (x,y,i,j,d2,hist) -> build_histogram!(d2,hist), zeros(Int,10), x, box)
    @test pairwise!((pair,hist) -> build_histogram!(pair.d2,hist), zeros(Int,10), box, cl, parallel=true) ≈ naive
    @test pairwise!((pair,hist) -> build_histogram!(pair.d2,hist), zeros(Int,10), box, cl, parallel=false) ≈ naive

    # Function to be evaluated for each pair: gravitational potential
    function potential(pair,u,mass)
        (; i, j, d2) = pair
        d2 == 0.0 && return u
        u = u - 9.8*mass[i]*mass[j]/pair.d
        return u
    end

    function potential(i,j,d2,u,mass)
        d2 == 0.0 && return u
        d = sqrt(d2)
        u = u - 9.8*mass[i]*mass[j]/d
        return u
    end

    # Run pairwise computation
    mass = rand(N)
    naive = map_naive!((x,y,i,j,d2,u) -> potential(i,j,d2,u,mass),0.0,x,box)
    @test pairwise!((pair,u) -> potential(pair,u,mass),0.0,box,cl,parallel=true) ≈ naive
    @test pairwise!((pair,u) -> potential(pair,u,mass),0.0,box,cl,parallel=false) ≈ naive

    # Function to be evaluated for each pair: gravitational force
    function calc_forces!(x,y,i,j,d2,mass,forces)
        d2 == 0.0 && return forces
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
    naive = map_naive!( (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces), copy(forces),x,box)
    @test pairwise!((pair,forces) -> calc_forces!(pair.x,pair.y,pair.i,pair.j,pair.d2,mass,forces), copy(forces),box,cl,parallel=true) ≈ naive
    @test pairwise!((pair,forces) -> calc_forces!(pair.x,pair.y,pair.i,pair.j,pair.d2,mass,forces), copy(forces),box,cl,parallel=false) ≈ naive

    # Test the examples, to check further if the parallelization didn't break something
    N = 100_000
    x = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:N ]
    @test Examples.average_displacement(parallel=true) ≈ Examples.average_displacement(parallel=false)[1]
    @test Examples.distance_histogram(parallel=true,x=x) ≈ Examples.distance_histogram(parallel=false,x=x)
    @test Examples.gravitational_potential(parallel=true,x=x) ≈ Examples.gravitational_potential(parallel=false,x=x)
    @test Examples.gravitational_force(parallel=true,x=x) ≈ Examples.gravitational_force(parallel=false,x=x)
    @test count(Examples.nearest_neighbor(parallel=true,x=x,y=y) .≈ Examples.nearest_neighbor(parallel=false,x=x,y=y)) == 3
    
    function pair_match(p1,p2) 
        p1[3] ≈ p2[3] || return false 
        p1[1] == p2[1] && p1[2] == p2[2] && return true
        p1[1] == p2[2] && p1[2] == p2[1] && return true
    end
    pairs1 = sort!(Examples.neighborlist(parallel=true,x=x),by=x->x[3])
    pairs2 = sort!(Examples.neighborlist(parallel=false,x=x),by=x->x[3])
    @test length(pairs1) == length(pairs2)
    @test count(pair_match.(pairs1,pairs2)) == length(pairs1)

    x = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500 ]
    y = [ sides .* rand(SVector{3,Float64}) for i in 1:1_500_000 ]
    @test count(Examples.nearest_neighbor_nopbc(parallel=true,x=x,y=y) .≈ Examples.nearest_neighbor_nopbc(parallel=false,x=x,y=y)) == 3

    # invert x and y to test swap
    ixy = Examples.nearest_neighbor_nopbc(parallel=false,x=x,y=y) 
    iyx = Examples.nearest_neighbor_nopbc(parallel=false,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 
    ixy = Examples.nearest_neighbor_nopbc(parallel=true,x=x,y=y) 
    iyx = Examples.nearest_neighbor_nopbc(parallel=true,x=y,y=x) 
    @test ( ixy[1] == iyx[2] && ixy[2] == iyx[1] && ixy[3] ≈ iyx[3] ) 

    # Test some fractional box lengths with the packmol test
    @test Examples.packmol(parallel=false, sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=CellListMap.TriclinicCell)[1]
    @test Examples.packmol(parallel=true, sides=[18.4,30.1,44], tol=2, UnitCellType=CellListMap.TriclinicCell)[1]
    @test Examples.packmol(parallel=false, sides=[46.4,32.1,44.7], tol=3.14, UnitCellType=CellListMap.OrthorhombicCell)[1]
    @test Examples.packmol(parallel=true, sides=[18.4,30.1,44], tol=2, UnitCellType=CellListMap.OrthorhombicCell)[1]
    
    # Testing the propagation of types in automatic differentiation
    @test Examples.generic_types(false) == (true,u"nm^2",Measurement{Float64})
    @test Examples.generic_types_triclinic(false) == (true,u"nm^2",Measurement{Float64})

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

@testitem "pathological cells" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    # This function is an interesting function because it sort of counts 
    # the indexes of the particles within the cutoff. However, we need to 
    # check for numerical precision inaccuracies before adding them, otherwise
    # pairs that match in one test can be skipped in another. This kind of issue
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
        box = CellListMap.Box(m,0.2)
        for x in [
            100 .* rand(SVector{2,Float64}, 100),
            [SVector{2,Float64}(0.1 * i + 0.1*j, 0.1 * j + 0.1*j) for i in 0:5 for j in 0:5],
            [SVector{2,Float64}(0.1 * i, 0.1 * j) for i in 0:5 for j in 0:5]
        ]
            cl = CellListMap.CellList(x, box)
            d0 = pairwise!((pair, out) -> f(pair.i, pair.j, pair.d2, out), 0, box, cl)
            d1 = map_naive!((x, y, i, j, d2, out) -> f(i,j,d2,out), 0, x, box)
            @test d0 ≈ d1
        end
    end

    # sph test (https://github.com/m3g/CellListMap.jl/issues/95)
    using StaticArrays: SVector
    p2d = SVector{2, Float64}[[0.0, 2.52], [0.02, 2.56], [3.98, 2.96], [4.0, 0.26], [4.0, 2.5]]
    r = 0.06788225099390856
    @test length(neighborlist(p2d, r)) == 1

end

@testitem "particle index uniqueness and periodic images" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    
    # Test 1: Verify particle pairs near periodic boundaries are counted exactly once
    # Two particles separated such that they interact through PBC
    count_interactions(x, y, i, j, d2, counter) = counter + 1
    
    sides = [10.0, 10.0, 10.0]
    cutoff = 2.0
    x = [SVector{3,Float64}(0.5, 5.0, 5.0), SVector{3,Float64}(9.6, 5.0, 5.0)]
    box = CellListMap.Box(sides, cutoff)
    
    # These particles are ~0.9 apart through PBC, should interact exactly once
    n_cell = pairwise!((pair, counter) -> counter + 1, 0, box, CellListMap.CellList(x, box))
    n_naive = map_naive!(count_interactions, 0, x, box)
    
    @test n_cell == n_naive
    @test n_cell == 1  # Should find exactly one pair
    
    # Test 2: Triclinic cell - verify interactions counted exactly once
    # even with complex periodic wrapping
    unit_cell = [10.0 3.0 0.0; 0.0 10.0 2.0; 0.0 0.0 10.0]
    cutoff = 2.5
    box = CellListMap.Box(unit_cell, cutoff)
    
    # Place particles near boundaries where triclinic wrapping matters
    x = [SVector{3,Float64}(0.5, 0.5, 0.5), 
         SVector{3,Float64}(2.0, 1.0, 1.0),
         SVector{3,Float64}(9.5, 9.5, 9.5),
         SVector{3,Float64}(8.5, 9.0, 9.0)]
    
    cl = CellListMap.CellList(x, box)
    n_cell = pairwise!((pair, counter) -> counter + 1, 0, box, cl)
    n_naive = map_naive!(count_interactions, 0, x, box)
    
    @test n_cell == n_naive
    
    # Test 3: Verify no duplicate pairs are computed
    # Collect all pairs and check for duplicates (regardless of order)
    function collect_pairs_naive(x, y, i, j, d2, pairs)
        # Normalize pair to (min, max) for duplicate detection
        p = i < j ? (i, j) : (j, i)
        push!(pairs, p)
        return pairs
    end
    
    x = [SVector{3,Float64}(rand(3) .* 20.0...) for _ in 1:50]
    box = CellListMap.Box([20.0, 20.0, 20.0], 3.0)
    cl = CellListMap.CellList(x, box)
    
    pairs = pairwise!((pair, pairs) -> begin
        p = pair.i < pair.j ? (pair.i, pair.j) : (pair.j, pair.i)
        push!(pairs, p)
        return pairs
    end, Tuple{Int,Int}[], box, cl, parallel=false)
    
    # Verify no duplicate pairs (each unique pair appears exactly once)
    @test length(pairs) == length(unique(pairs))
    
    # Verify count matches naive implementation
    n_pairs = length(pairs)
    n_naive = pairwise!((pair,c) -> c + 1, 0, box, cl, parallel=false)
    @test n_pairs == n_naive
    
    # Test 4: Same test with triclinic cell
    unit_cell = [20.0 5.0 0.0; 0.0 20.0 3.0; 0.0 0.0 20.0]
    box = CellListMap.Box(unit_cell, 3.0)
    x = [SVector{3,Float64}((unit_cell * rand(3))...) for _ in 1:50]
    cl = CellListMap.CellList(x, box)
    
    pairs = pairwise!((pair, pairs) -> begin
        p = pair.i < pair.j ? (pair.i, pair.j) : (pair.j, pair.i)
        push!(pairs, p)
        return pairs
    end, Tuple{Int,Int}[], box, cl, parallel=false)
    
    @test length(pairs) == length(unique(pairs))
    
    # Test 5: Dense system near boundaries - stress test for duplicate prevention
    # Many particles near boundaries where PBC creates many potential duplicates
    sides = [15.0, 15.0, 15.0]
    cutoff = 2.5
    box = CellListMap.Box(sides, cutoff)
    
    # Create particles clustered near boundaries
    x = SVector{3,Float64}[]
    for i in 1:5, j in 1:5, k in 1:5
        push!(x, SVector{3,Float64}(0.3 * i, 0.3 * j, 0.3 * k))
        push!(x, SVector{3,Float64}(sides[1] - 0.3 * i, sides[2] - 0.3 * j, sides[3] - 0.3 * k))
    end
    
    cl = CellListMap.CellList(x, box)
    n_cell = pairwise!((pair, counter) -> counter + 1, 0, box, cl, parallel=false)
    n_naive = map_naive!(count_interactions, 0, x, box)
    
    @test n_cell == n_naive
end

@testitem "margin calculation with lcell > 1" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    
    # Test that particles at cell corners are not missed with lcell > 1
    # This tests the margin calculation in project_particles!
    sum_distances_naive(x, y, i, j, d2, s) = s + sqrt(d2)
    
    # Test various lcell values
    for lcell in [1, 2, 3, 5]
        # Create a system where particles are positioned at strategic locations
        # including near cell corners
        sides = [20.0, 20.0, 20.0]
        cutoff = 2.5
        box = CellListMap.Box(sides, cutoff, lcell=lcell)
        
        # Generate particles including some at corners of cells
        N = 200
        x = [SVector{3,Float64}(sides .* rand(3)...) for _ in 1:N]
        
        # Add specific particles near cell boundaries to test corner cases
        cell_size = sides ./ (2 * lcell + 1)
        for i in 1:lcell+1
            for j in 1:lcell+1
                for k in 1:lcell+1
                    # Place particle at corner of cell (i,j,k)
                    corner = SVector{3,Float64}(
                        cell_size[1] * i + 0.1,
                        cell_size[2] * j + 0.1,
                        cell_size[3] * k + 0.1
                    )
                    if all(corner .< sides)
                        push!(x, corner)
                    end
                end
            end
        end
        
        cl = CellListMap.CellList(x, box)
        
        # Compare with naive implementation
        result_cell = pairwise!((pair, s) -> s + pair.d, 0.0, box, cl, parallel=false)
        result_naive = map_naive!(sum_distances_naive, 0.0, x, box)
        
        @test result_cell ≈ result_naive rtol=1e-10
    end
    
    # Test 2: Triclinic cells with lcell > 1
    for lcell in [1, 2, 3]
        unit_cell = [15.0 5.0 0.0; 0.0 15.0 3.0; 0.0 0.0 15.0]
        cutoff = 2.0
        box = CellListMap.Box(unit_cell, cutoff, lcell=lcell)
        
        N = 150
        # Generate random fractional coordinates and convert to Cartesian
        x = [SVector{3,Float64}((unit_cell * rand(3))...) for _ in 1:N]
        
        cl = CellListMap.CellList(x, box)
        
        result_cell = pairwise!((pair, s) -> s + pair.d, 0.0, box, cl, parallel=false)
        result_naive = map_naive!(sum_distances_naive, 0.0, x, box)
        
        @test result_cell ≈ result_naive rtol=1e-10
    end
    
    # Test 3: Verify no interactions are missed at exactly the cutoff distance
    # with lcell > 1 and particles positioned to stress-test the margin calculation
    cutoff = 3.0
    lcell = 3
    sides = [30.0, 30.0, 30.0]
    box = CellListMap.Box(sides, cutoff, lcell=lcell)
    
    # Create pairs of particles at exactly cutoff distance apart
    # positioned across cell boundaries
    x = SVector{3,Float64}[]
    push!(x, SVector{3,Float64}(5.0, 5.0, 5.0))
    push!(x, SVector{3,Float64}(5.0 + cutoff - 0.001, 5.0, 5.0))  # Just within cutoff
    push!(x, SVector{3,Float64}(10.0, 10.0, 10.0))
    push!(x, SVector{3,Float64}(10.0, 10.0 + cutoff - 0.001, 10.0))  # Just within cutoff
    push!(x, SVector{3,Float64}(15.0, 15.0, 15.0))
    push!(x, SVector{3,Float64}(15.0, 15.0, 15.0 + cutoff - 0.001))  # Just within cutoff
    
    cl = CellListMap.CellList(x, box)
    
    n_cell = pairwise!((pair,c) -> c + 1, 0, box, cl, parallel=false)
    n_naive = map_naive!((x,y,i,j,d2,c) -> c + 1, 0, x, box)
    
    @test n_cell == n_naive
    @test n_cell >= 3  # At least the 3 pairs we explicitly created
end

@testitem "random cells - Orthorhombic" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    # Test random cells of all possible types
    @test test_random_cells(CellListMap.OrthorhombicCell) == (nothing, nothing, nothing)
end

@testitem "random cells - Triclinic" setup=[Testing] begin
    using CellListMap
    using StaticArrays
    # Test random cells of all possible types
    @test test_random_cells(CellListMap.TriclinicCell) == (nothing, nothing, nothing)
end

@testitem "test_pathological2D" setup=[Testing] begin
    @test test_pathological() == (nothing, nothing)
end
