#
# Function that uses the naive algorithm, for testing
#
function map_naive!(f,output,x,box)
  @unpack unit_cell, cutoff_sq = box
  for i in 1:length(x)-1
    xᵢ = x[i]
    for j in i+1:length(x)
      xⱼ = wrap_relative_to(x[j],xᵢ,unit_cell)
      d2 = norm_sqr(xᵢ - xⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,xⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

#
# Function that uses the naive algorithm, for testing
#
function map_naive_two!(f,output,x,y,box)
  @unpack unit_cell, cutoff_sq = box
  for i in 1:length(x)
    xᵢ = x[i]
    for j in 1:length(y)
      yⱼ = wrap_relative_to(y[j],xᵢ,unit_cell)
      d2 = norm_sqr(xᵢ - yⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,yⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

"""

```
view_celllist_particles(cl::CellList)
```

Auxiliary function to view the particles of a computing box, including images created
for computing purposes.

### Example
```julia
julia> box = Box([ 100 50; 50 100 ],10);

julia> p = [ box.unit_cell_max .* rand(SVector{2,Float64}) for i in 1:1000 ];

julia> cl = CellList(p,box);

julia> x, y = CellListMap.view_celllist_particles(cl);

julia> using Plots

julia> scatter(x,y,label=nothing,xlims=(-10,180),ylims=(-10,180))

```

"""
function view_celllist_particles(cl::CellList{N,T},box::Box) where {N,T}
  @unpack nc = box
  @unpack cwp, ncp, fp, np = cl
  x = Vector{SVector{N,T}}(undef,ncp)
  ip = 0
  for p in cl.fp
    p.index == 0 && continue
    ip += 1
    x[ip] = p.coordinates
    while np[p.index].index > 0
      ip += 1
      x[ip] = p.coordinates
      p = np[p.index]
    end
  end
  return [SVector{N,T}(ntuple(j -> x[i][j],N)) for i in 1:ncp]
end

