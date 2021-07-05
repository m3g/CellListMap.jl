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
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  N = 2000
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  # Add some pathological coordinates
  x[1] = @SVector([-sides[1], -sides[2], -sides[3]/2])
  x[10] = @SVector([sides[1], sides[2], 3*sides[3]])
  x[100] = @SVector([sides[1]/2, -sides[2]/2, 2*sides[3]])
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # Initialize auxiliary linked lists
  cl = CellList(x,box)

  # Function to be evalulated for each pair: sum of displacements on x
  f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

  naive = abs(CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box))
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=false)) ≈ naive

  # Check if changing lcell breaks something
  box = Box(sides,cutoff,lcell=1); cl = CellList(x,box)
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive
  box = Box(sides,cutoff,lcell=2); cl = CellList(x,box)
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive
  box = Box(sides,cutoff,lcell=3); cl = CellList(x,box)
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive
  box = Box(sides,cutoff,lcell=15); cl = CellList(x,box)
  @test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive

  # Function to be evalulated for each pair: build distance histogram
  function build_histogram!(x,y,d2,hist)
    d = sqrt(d2)
    ibin = floor(Int,d) + 1
    hist[ibin] += 1
    return hist
  end

  naive = CellListMap.map_naive!(
    (x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),
    zeros(Int,10),x,box
  )
  @test map_pairwise!(
    (x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),
    zeros(Int,10),box,cl,parallel=true
  ) ≈ naive
  @test map_pairwise!(
    (x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),
    zeros(Int,10),box,cl,parallel=false
  ) ≈ naive

  # Function to be evalulated for each pair: build distance histogram
  function potential(x,y,i,j,d2,u,mass)
    d = sqrt(d2)
    u = u - 9.8*mass[i]*mass[j]/d
    return u
  end

  # Run pairwise computation
  mass = rand(N)
  naive = CellListMap.map_naive!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box)
  @test map_pairwise!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,box,cl,parallel=true) ≈ naive
  @test map_pairwise!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,box,cl,parallel=false) ≈ naive

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

  naive = CellListMap.map_naive_two!((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,y,box)
  @test map_pairwise!(
    (x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),
    0.0,box,cl,parallel=false
  ) ≈ naive
  @test map_pairwise!(
    (x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),
    0.0,box,cl,parallel=true
  ) ≈ naive

  # Check different lcell
  box = Box(sides,cutoff,lcell=3)
  cl = CellList(x,y,box)
  @test map_pairwise!(
    (x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),
    0.0,box,cl,parallel=false
  ) ≈ naive
  @test map_pairwise!(
    (x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),
    0.0,box,cl,parallel=true
  ) ≈ naive

  # Test the examples, to check further if the parallelization didn't break something
  N = 100_000
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  @test CellListMap.test1(parallel=true,x=x) ≈ CellListMap.test1(parallel=false,x=x)
  @test CellListMap.test2(parallel=true,x=x) ≈ CellListMap.test2(parallel=false,x=x)
  @test CellListMap.test3(parallel=true,x=x) ≈ CellListMap.test3(parallel=false,x=x)
  @test CellListMap.test4(parallel=true,x=x) ≈ CellListMap.test4(parallel=false,x=x)
  @test count(CellListMap.test5(parallel=true,x=x,y=y) .≈ CellListMap.test5(parallel=false,x=x,y=y)) == 3

  pairs1 = sort!(CellListMap.test7(parallel=true,x=x),by=x->x[3])
  pairs2 = sort!(CellListMap.test7(parallel=false,x=x),by=x->x[3])
  @test count([ count(pairs1[i] .≈ pairs2[i]) == 3 for i in 1:length(pairs1) ]) == length(pairs1)

  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:1_500 ]
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:1_500_000 ]
  @test count(CellListMap.test6(parallel=true,x=x,y=y) .≈ CellListMap.test6(parallel=false,x=x,y=y)) == 3

  # invert to see if swap is working as expected
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:1_500_000 ]
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:1_500 ]
  @test count(CellListMap.test6(parallel=true,x=x,y=y) .≈ CellListMap.test6(parallel=false,x=x,y=y)) == 3

  # Test resizing of the cell lists
  x = [ rand(SVector{3,Float64}) for i in 1:1000 ]
  box = Box([0.83,0.41,0.97],0.1)
  cl = CellList(x,box) 
  @test length(cl.cwp) == 317

  box = Box([0.33,0.41,0.97],0.1)
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 317

  box = Box([0.83,0.81,0.97],0.1)
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 634 

  x .= 0.9*x
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 634 

  x .= 1.2*x
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 634 

  box = Box([0.83,0.81,0.97],0.2)
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 634 

  box = Box([0.83,0.81,0.97],0.05)
  cl = UpdateCellList!(x,box,cl)   
  @test length(cl.cwp) == 5351 

end
