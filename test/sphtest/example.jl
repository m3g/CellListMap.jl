import Pkg; Pkg.activate(".")
using CSV
using CellListMap
using StaticArrays
using DataFrames
using LinearAlgebra: norm

path    = @__DIR__
DF_FLUID = CSV.read(joinpath(path, "./FluidPoints_Dp0.02.csv"), DataFrame)
DF_BOUND = CSV.read(joinpath(path, "./BoundaryPoints_Dp0.02.csv"), DataFrame)
P1F = DF_FLUID[!,"Points:0"]
P2F = DF_FLUID[!,"Points:1"]
P3F = DF_FLUID[!,"Points:2"]
P1B = DF_BOUND[!,"Points:0"]
P2B = DF_BOUND[!,"Points:1"]
P3B = DF_BOUND[!,"Points:2"]
points = Vector{SVector{3,Float64}}()
for i = 1:length(P1F)
    push!(points,SVector(P1F[i],P3F[i],P2F[i]))
end
for i = 1:length(P1B)
    push!(points,SVector(P1B[i],P3B[i],P2B[i]))
end

H = 0.06788225099390856
system  = InPlaceNeighborList(x=points, cutoff = H, parallel=true)
update!(system, points)
list = neighborlist!(system)

function vprop(p, r)
    box = Box(limits(p), r)
    cl = CellList(p, box)
    vprop = map_pairwise(
        (x,y,i,j,d2,u) -> u += sqrt(d2) - r,
        0.0, box, cl
    )
    return vprop
end

function brute_force(p, r)
    ps = Vector{Tuple{Int, Int}}()
    n = length(p)
    vprop = 0.0
    for i in 1:(n-1)
        for j in (i+1):n
            d = norm(p[j] - p[i])
            if d <= r
                push!(ps, (i, j))
                vprop += d - r
            end
        end
    end
    return ps, vprop
end

list_brute_force, vprop_brute_force = brute_force(points, H)

println("length(list) = ", length(list))
println("length(list_brute_force) = ", length(list_brute_force))

vprop_celllists = vprop(points, H)

println("vprop by CellListMap = ", vprop_celllists)
println("vprop by brute force= ", vprop_brute_force)











