#
# To compare with halotools
#
function halotools(;N=100_000,cd=true,parallel=true)

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
  
  Lbox = [boxsize,boxsize,boxsize]
  rbins = [0.,2.,4.,6.,8.,10.]
  r_max = maximum(rbins)

  positions = boxsize*rand(SVector{3,Float64},n_halos)
  velocities = rand(SVector{3,Float64},n_halos)

  box = Box(Lbox, r_max)
  cl = CellList(positions,box,parallel=parallel)
  hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))

  hist = map_pairwise!(
    (x,y,i,j,d2,hist) -> compute_pairwise_mean_cell_lists!(
       x,y,i,j,d2,hist,velocities,rbins,Lbox
    ),
    hist, box, cl,
    reduce=reduce_hist,
    parallel=parallel
  )

  n_pairs = hist[1]
  mean_v_r = hist[2]
  mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0]./n_pairs[n_pairs .> 0]
  return mean_v_r

end


