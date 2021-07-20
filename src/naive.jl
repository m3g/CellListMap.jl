#
# Function that uses the naive algorithm, for testing
#
function map_naive!(f,output,x,box)
  @unpack unit_cell, cutoff_sq = box
  for i in 1:length(x)-1
    xᵢ = x[i]
    for j in i+1:length(x)
      xⱼ = wrap_relative_to(x[j],xᵢ,unit_cell)
      d2 = distance_sq(xᵢ,xⱼ)
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
      d2 = distance_sq(xᵢ,yⱼ)
      if d2 <= cutoff_sq
        output = f(xᵢ,yⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

#
# naive filling of computing box, for testing
#
function replicate_particle0!(p::SVector{2,T},box,cl) where T
  ri = -100:100
  rj = -100:100 
  for i in ri
    for j in rj
      x = translation_image(p,box.unit_cell,(i,j))
      if ! out_of_bounding_box(x,box)
        add_particle_to_celllist!(x,box,cl) 
      end
    end
  end 
  return nothing
end

function replicate_particle0!(p::SVector{3,T},box,cl) where T
  ri = -100:100
  rj = -100:100 
  rk = -100:100 
  for i in ri
    for j in rj
      for k in rk
        x = translation_image(p,box.unit_cell,(i,j,k))
        if ! out_of_bounding_box(x,box)
          add_particle_to_celllist!(x,box,cl) 
        end
      end
    end
  end 
  return nothing
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
function view_celllist_particles(cl::CellList{V,N,T}) where {V,N,T}
  @unpack ncwp, cwp, ncp, fp, np = cl
  x = Vector{SVector{N,T}}(undef,ncp[1])
  ip = 0
  for i in 1:ncwp[1]
    ip += 1
    p = fp[cwp[i].icell]
    x[ip] = p.coordinates
    while np[p.index].index > 0
      ip += 1
      x[ip] = p.coordinates
      p = np[p.index]
    end
  end
  return ([x[i][j] for i in 1:ncp[1]] for j in 1:N)
end

