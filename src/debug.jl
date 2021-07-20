import Random
using StaticArrays, Plots, BenchmarkTools, Revise; using CellListMap

Random.seed!(321)

# 2D
box = Box([ 100   0 
              0 100 ],10)
p = [ 10*box.unit_cell_max .* rand(SVector{2,Float64}) for i in 1:1000 ];
cl = CellList(p,box)
x, y = CellListMap.view_celllist_particles(cl);
scatter(x,y,label=nothing,xlims=(-10,210),ylims=(-10,210),markersize=0.1,aspect_ratio=1)
  
# florpi 
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
    cl = CellList(positions,box)
    hist = (zeros(Int,length(rbins)-1), zeros(Float64,length(rbins)-1))
  
    # Needs this to stabilize the type of velocities and hist, probably
    function barrier(f,velocities,rbins,Lbox,hist,positions,box,cl,reduce_hist,parallel)
      hist = CellListMap.map_pairwise_serial!(
        (x,y,i,j,d2,hist) -> compute_pairwise_mean_cell_lists!(
           x,y,i,j,d2,hist,velocities,rbins,Lbox
        ),
        hist, box, cl,
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

x, y, z = CellListMap.view_celllist_particles(cl);
scatter(x,y,z,label=nothing,lims=(-20,270),markersize=0.1,aspect_ratio=1)
  