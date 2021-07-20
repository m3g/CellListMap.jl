import Random

#
# In this test we compute the average displacement of the x coordinates of the atoms
#              
function test1(;N=100_000,parallel=true,x=nothing)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ SVector{3,Float64}(sides...) .* rand(SVector{3,Float64}) for i in 1:N ]
  end

  # Initialize auxiliary linked lists
  cl = CellList(x,box,parallel=parallel)

  # Function to be evalulated for each pair: sum of displacements on x
  f(x,y,avg_dx) = avg_dx + x[1] - y[1]

  avg_dx = (N/(N*(N-1)/2)) * map_pairwise!(
    (x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),
    0.,box,cl,
    parallel=parallel
  )
  return avg_dx

end

#
# In this test we compute the histogram of distances, expected to follow the
# function f(f) = ρ(4/3)π(r[i+1]^3 - r[i]^3) with ρ being the density of the system.
#
function test2(;N=100_000,parallel=true,x=nothing)

  # Number of particles, sides and cutoff
  sides = SVector{3,Float64}(250,250,250)
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  end

  # Initialize auxiliary linked lists
  cl = CellList(x,box,parallel=parallel)

  # Function to be evalulated for each pair: build distance histogram
  function build_histogram!(x,y,d2,hist) 
    d = sqrt(d2)
    ibin = floor(Int,d) + 1
    hist[ibin] += 1
    return hist
  end

  # Preallocate and initialize histogram
  hist = zeros(Int,10)

  # Run pairwise computation
  hist = (N/(N*(N-1)/2)) * map_pairwise!(
    (x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),
    hist,box,cl,
    parallel=parallel
  )
  return hist

end

#
# In this test we compute the "gravitational potential", pretending that each particle
# has a different mass. In this case, the closure is used to pass the masses to the
# function that computes the potential.
#
function test3(;N=100_000,parallel=true,x=nothing)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  end

  # Initialize auxiliary linked lists
  cl = CellList(x,box,parallel=parallel)

  # masses
  mass = [ 5 * x[i][1] for i in 1:N ]

  # Function to be evalulated for each pair: build distance histogram
  function potential(x,y,i,j,d2,u,mass) 
    d = sqrt(d2)
    u = u - 9.8*mass[i]*mass[j]/d
    return u
  end

  # Run pairwise computation
  u = map_pairwise!(
    (x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),
    0.0,box,cl,
    parallel=parallel
  )
  return u

end

#
# In this test we compute the "gravitational force", pretending that each particle
# has a different mass. In this case, the closure is used to pass the masses and
# the force vector to the function that computes the potential.
#
function test4(;N=100_000,parallel=true,x=nothing)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  end

  # masses
  mass = [ 5 * x[i][1] for i in 1:N ]

  # Initialize auxiliary linked lists
  cl = CellList(x,box,parallel=parallel)

  # Function to be evalulated for each pair
  function calc_forces!(x,y,i,j,d2,mass,forces) 
    G = 9.8*mass[i]*mass[j]/d2
    d = sqrt(d2)
    df = (G/d) * (x - y)
    forces[i] = forces[i] - df
    forces[j] = forces[j] + df
    return forces
  end

  # Preallocate and initialize forces
  forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

  # Run pairwise computation
  forces = map_pairwise!(
    (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
    forces,box,cl,
    parallel=parallel
  )
  return forces

end

#
# In this test we compute the minimum distance between two independent sets of particles
#
function test5(;N1=1_500,N2=1_500_000,parallel=true,x=nothing,y=nothing)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N1 ]
  end
  if y === nothing 
    y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N2 ]
  end

  # Initialize auxiliary linked lists (largest set!)
  cl = CellList(x,y,box,parallel=parallel)

  # Function that keeps the minimum distance
  f(i,j,d2,mind) = d2 < mind[3] ? (i,j,d2) : mind

  # We have to define our own reduce function here (for the parallel version)
  function reduce_mind(output,output_threaded)
    mind = output_threaded[1]
    for i in 2:nthreads()
      if output_threaded[i][3] < mind[3]
        mind = output_threaded[i]
      end
    end
    return mind
  end 

  # Initialize
  mind = ( 0, 0, +Inf )

  # Run pairwise computation
  mind = map_pairwise!(
    (x,y,i,j,d2,mind) -> f(i,j,d2,mind),
    mind,box,cl;reduce=reduce_mind, parallel=parallel
  )
  return (mind[1],mind[2],sqrt(mind[3]))

end

#
# In this test we compute the minimum distance between two independent sets of particles,
# without periodic conditions
#
function test6(;N1=1_500,N2=1_500_000,parallel=true,x=nothing,y=nothing)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ rand(SVector{3,Float64}) for i in 1:N1 ]
  else
    N1 = length(x)
  end
  if y === nothing 
    y = [ rand(SVector{3,Float64}) for i in 1:N2 ]
  else
    N2 = length(y)
  end

  # Boundaries
  xmin = [ +Inf, +Inf, +Inf ]
  xmax = [ -Inf, -Inf, -Inf ]
  for v in x
    @. xmin = min(xmin,v)
    @. xmax = max(xmax,v)
  end
  for v in y
    @. xmin = min(xmin,v)
    @. xmax = max(xmax,v)
  end
   
  # Obtain one upper bound for dmin by computing one distance for each element
  # of the smallest vector
  cutoff = +Inf
  for v in x
    iy = rand(1:N2)
    cutoff = min(CellListMap.distance(v,y[iy]),cutoff)
  end 
   
  # Define box sides
  sides = (xmax - xmin) .+ cutoff
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists (largest set!)
  cl = CellList(x,y,box,parallel=parallel)

  # Function that keeps the minimum distance
  f(i,j,d2,mind) = d2 < mind[3] ? (i,j,d2) : mind

  # We have to define our own reduce function here (for the parallel version)
  function reduce_mind(output,output_threaded)
    mind = output_threaded[1]
    for i in 2:Threads.nthreads()
      if output_threaded[i][3] < mind[3]
        mind = output_threaded[i]
      end
    end
    return mind
  end 

  # Initialize 
  mind = ( 0, 0, +Inf )

  # Run pairwise computation
  mind = map_pairwise!(
    (x,y,i,j,d2,mind) -> f(i,j,d2,mind),
    mind,box,cl;reduce=reduce_mind,parallel=parallel
  )
  return (mind[1],mind[2],sqrt(mind[3]))

end

#
# In this test we compute the complete neighbor list of particles, meaning all the pairs
# that are within the cutoff distance
#
function test7(;N=100_000,parallel=true,x=nothing)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Particle positions
  Random.seed!(321)
  if x === nothing 
    x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]
  end

  # Initialize auxiliary linked lists
  cl = CellList(x,box,parallel=parallel)

  # Function to be evalulated for each pair: push pair if d<cutoff
  function push_pair!(i,j,d2,pairs,cutoff) 
    d = sqrt(d2)
    if d < cutoff
      push!(pairs,(i,j,d))
    end
    return pairs
  end

  # Reduction function
  function reduce_pairs(pairs,pairs_threaded)
    pairs = pairs_threaded[1]
    for i in 2:nthreads()
      append!(pairs,pairs_threaded[i])
    end
    return pairs
  end

  # Initialize output
  pairs = Tuple{Int,Int,Float64}[]

  # Run pairwise computation
  pairs = map_pairwise!(
    (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs,cutoff),
    pairs,box,cl,
    reduce=reduce_pairs,
    parallel=parallel
  )
  return pairs

end

#
# florpi
#
function florpi(;N=100_000,cd=true,parallel=true)

  @inline dot(x::SVector{3,Float64},y::SVector{3,Float64}) = x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
  
  function compute_pairwise_mean_cell_lists!(x,y,i,j,d2,hist,velocities,rbins,sides)
    d = x - y
    r = sqrt(d2)
    ibin = searchsortedfirst(rbins, r) - 1
    hist[1][ibin] += 1
    hist[2][ibin] += dot(velocities[i]-velocities[j],d)/r
    return hist
  end

  function reduce_hist(hist,hist_threaded)
    hist = hist_threaded[1]
    for i in 2:Threads.nthreads()
      hist[1] .+= hist_threaded[i][1]
      hist[2] .+= hist_threaded[i][2]
    end
    return hist
  end

  n_halos = N

  if cd
    density = 10^5/250^3  # density of the original problem
    boxsize = (n_halos / density)^(1/3)
  else
    boxsize = 250.
  end
  
  Random.seed!(321)
  Lbox = [boxsize,boxsize,boxsize]
  positions = boxsize .* rand(Float64, 3, n_halos)
  velocities = rand(Float64, 3, n_halos)
  rbins = [0.,2.,4.,6.,8.,10.]
  r_max = maximum(rbins)

  n = size(positions)[2]
  positions = reshape(reinterpret(SVector{3,Float64},positions),n)
  velocities = reshape(reinterpret(SVector{3,Float64},velocities),n)

  box = Box(Lbox, r_max)
  cl = CellList(positions,box,parallel=parallel)
  hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))

  # Needs this to stabilize the type of velocities and hist, probably
  function barrier(f,velocities,rbins,Lbox,hist,positions,box,cl,reduce_hist,parallel)
    hist = map_pairwise!(
      (x,y,i,j,d2,hist) -> compute_pairwise_mean_cell_lists!(
         x,y,i,j,d2,hist,velocities,rbins,Lbox
      ),
      hist, box, cl,
      reduce=reduce_hist,
      parallel=parallel
    )
    return hist
  end

  hist = barrier(compute_pairwise_mean_cell_lists!,
    velocities,rbins,Lbox,hist,positions,box,cl,reduce_hist,parallel)

  n_pairs = hist[1]
  mean_v_r = hist[2]
  mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0]./n_pairs[n_pairs .> 0]
  return mean_v_r

end
