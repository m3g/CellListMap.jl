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

# Number of particles, sides and cutoff
sides = [250,250,250]
cutoff = 10.
box = Box(sides,cutoff)

# Particle positions
N = 2000
x = [ box.unit_cell_max .* rand(SVector{3,Float64}) for i in 1:N ]
# Add some pathological coordinates
x[1] = @SVector([-sides[1], -sides[2], -sides[3]/2])
x[10] = @SVector([sides[1], sides[2], 3*sides[3]])
x[100] = @SVector([sides[1]/2, -sides[2]/2, 2*sides[3]])

# Initialize auxiliary linked lists
cl = CellList(x,box)

# Function to be evalulated for each pair: sum of displacements on x
f(x,y,avg_dx) = avg_dx + abs(x[1] - y[1])

naive = abs(CellListMap.map_naive!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box))
CellListMap.map_pairwise_serial!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl)

@test abs(map_pairwise!((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,box,cl,parallel=true)) ≈ naive
