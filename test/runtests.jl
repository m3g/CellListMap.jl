using CellListMap
using StaticArrays
using Test

@testset "CellListMap.jl" begin

  # Number of particles, sides and cutoff
  N = 2000
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists
  lc = LinkedLists(N)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # Initializing linked cells with these positions
  initcells!(x,box,lc)

  # Function to be evalulated for each pair: sum of displacements on x
  f(x,y,avg_dx) = avg_dx + x[1] - y[1]

  naive = abs(CellListMap.map_naive((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box))
  @test abs(map_pairwise((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box,lc,parallel=true)) ≈ naive
  @test abs(map_pairwise((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box,lc,parallel=false)) ≈ naive

  # Function to be evalulated for each pair: build distance histogram
  function build_histogram!(x,y,d2,hist)
    d = sqrt(d2)
    ibin = floor(Int,d) + 1
    hist[ibin] += 1
    return hist
  end

  naive = CellListMap.map_naive((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),zeros(Int,10),x,box)
  @test map_pairwise((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),zeros(Int,10),x,box,lc,parallel=true) ≈ naive
  @test map_pairwise((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),zeros(Int,10),x,box,lc,parallel=false) ≈ naive

  # Function to be evalulated for each pair: build distance histogram
  function potential(x,y,i,j,d2,u,mass)
    d = sqrt(d2)
    u = u - 9.8*mass[i]*mass[j]/d
    return u
  end
  mass = rand(N)

  # Run pairwise computation
  naive = CellListMap.map_naive((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box)
  @test map_pairwise((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc,parallel=true) ≈ naive
  @test map_pairwise((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc,parallel=false) ≈ naive

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
  naive = CellListMap.map_naive((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),copy(forces),x,box)
  @test map_pairwise((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),copy(forces),x,box,lc,parallel=true) ≈ naive
  @test map_pairwise((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),copy(forces),x,box,lc,parallel=false) ≈ naive

end
