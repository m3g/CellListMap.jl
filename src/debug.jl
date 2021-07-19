import Random
using StaticArrays, Plots, BenchmarkTools, Revise; using CellListMap

Random.seed!(321)

# 2D
box = Box([ 100  90 
              0 100 ],10);
p = [ 10*box.unit_cell_max .* rand(SVector{2,Float64}) for i in 1:1000 ];
cl = CellList(p,box)
x, y = CellListMap.view_celllist_particles(cl);
scatter(x,y,label=nothing,xlims=(-10,210),ylims=(-10,210),markersize=0.1)
  
# 3D
box = Box([ 100   0   0 
              0 100   0
              0  90 100],10)
p = [ 10*box.unit_cell_max .* rand(SVector{3,Float64}) for i in 1:1000 ];
cl = CellList(p,box)
x, y, z = CellListMap.view_celllist_particles(cl);
scatter(x,y,z,label=nothing,xlims=(-10,210),ylims=(-10,210),markersize=1.0)
  