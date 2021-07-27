#
# Function that checks if the particule is outside the computation bounding box
#
function out_of_bounding_box(x::SVector{N,T},box::Box{N,T,M}) where {N,T,M}
  @unpack cutoff, unit_cell = box
  unit_cell_max = sum(@view(unit_cell[:,i]) for i in 1:N) 
  for i in 1:N
    (x[i] < -cutoff) && return true
    (x[i] >= unit_cell_max[i]+cutoff) && return true
  end
  return false
end
out_of_bounding_box(p::AtomWithIndex,box::Box{N,T,M}) where {N,T,M} =
  out_of_bounding_box(p.coordinates,box)

"""

```
replicate_particle!(ip,p::T,box,cl) where {T <: SVector{2,S} where S}
```

Replicates the particle as many times as necessary to fill the computing box.

"""
function replicate_particle!(ip,p::T,box,cl) where {T <: SVector{2,S} where S}
  @unpack ranges = box
  for i in ranges[1]
    for j in ranges[2]
      i == 0 && j == 0 && continue
      x = translation_image(p,box.unit_cell,(i,j))
      if ! out_of_bounding_box(x,box)
        cl = add_particle_to_celllist!(ip,x,box,cl;real_particle=false) 
      end
    end
  end 
  return cl
end


"""

```
replicate_particle!(ip,p::T,box,cl) where {T <: SVector{3,S} where S}
```

Replicates the particle as many times as necessary to fill the computing box.

"""
function replicate_particle!(ip,p::T,box,cl) where {T <: SVector{3,S} where S}
  @unpack ranges = box
  for i in ranges[1]
    for j in ranges[2]
      for k in ranges[3]
        i == 0 && j == 0 && k == 0 && continue
        x = translation_image(p,box.unit_cell,(i,j,k))
        if ! out_of_bounding_box(x,box)
          cl = add_particle_to_celllist!(ip,x,box,cl;real_particle=false) 
        end
      end
    end
  end 
  return cl
end

function ranges_of_replicas(cell_size, lcell, nc, unit_cell::SMatrix{3,3,T,9}) where T
  V = SVector{3,T}
  cmin = ntuple(i->-lcell*cell_size,3)
  cmax = ntuple(i->(nc[i]-lcell)*cell_size,3)
  cell_vertices = SVector{8,V}( 
    V(   cmin[1],   cmin[2],   cmin[3] ), 
    V(   cmin[1],   cmin[2],   cmax[3] ),
    V(   cmin[1],   cmax[2],   cmin[3] ),
    V(   cmin[1],   cmax[2],   cmax[3] ),
    V(   cmax[1],   cmin[2],   cmin[3] ), 
    V(   cmax[1],   cmin[2],   cmax[3] ), 
    V(   cmax[1],   cmax[2],   cmin[3] ), 
    V(   cmax[1],   cmax[2],   cmax[3] ) 
  )
  r_min, r_max = _ranges_of_replicas(
    nc,
    SVector{3,Int}(-lcell,-lcell,-lcell),
    unit_cell,
    cell_vertices
  )
  ranges = SVector{3,UnitRange{Int}}(
    r_min[1]:r_max[1],
    r_min[2]:r_max[2],
    r_min[3]:r_max[3]
  )
  return ranges
end

function ranges_of_replicas(lcell,cell_size,nc,unit_cell::SMatrix{2,2,T,4}) where T
  V = SVector{2,T}
  cmin = ntuple(i->-lcell*cell_size,2)
  cmax = ntuple(i->(nc[i]-lcell)*cell_size,2)
  cell_vertices = SVector{4,V}( 
    V(   cmin[1],   cmin[2] ), 
    V(   cmin[1],   cmax[2] ),
    V(   cmax[1],   cmin[2] ), 
    V(   cmax[1],   cmax[2] ), 
  )
  r_min, r_max = _ranges_of_replicas(
    nc,
    SVector{2,Int}(-lcell,-lcell),
    unit_cell,
    cell_vertices
  )
  ranges = SVector{2,UnitRange{Int}}(
    r_min[1]:r_max[1],
    r_min[2]:r_max[2]
  )
  return ranges
end

function _ranges_of_replicas(r_min,r_max,unit_cell,cell_vertices)
  for vert in cell_vertices
    r = unit_cell \ vert 
    ri = @. ceil(Int,abs(r))
    for (i,el) in pairs(r)
      if el < 0
        @set! ri[i] = -ri[i]
      end
    end
    r_min = min.(ri,r_min)
    r_max = max.(ri,r_max)
  end
  return r_min, r_max
end

"""

```
wrap_cell_fraction(x,cell)
```

`x` is a vector of dimension `N` and `cell` a matrix of dimension `NxN`

"""
@inline function wrap_cell_fraction(x,cell)
  p = rem.(cell\x,1)
  for i in eachindex(p)
    if p[i] < 0
      @set! p[i] += 1
    end
  end
  return p
end

"""

```
wrap_to_first(x::SVector{N,T},cell) where {N,T}
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates. 

"""
@inline function wrap_to_first(x,cell)
  p = wrap_cell_fraction(x,cell)
  return cell*p
end
"""

```
wrap_to_first(x::SVector{N,T},box::Box)
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates, given the `Box` structure.

"""
@inline wrap_to_first(x,box::Box) = wrap_to_first(x,box.unit_cell)

"""

```
wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

"""
function wrap_relative_to(x::T, xref, cell) where T<:AbstractVector
  p = SVector{length(x),eltype(x)}(rem.(cell\(x-xref),1))
  for i in eachindex(p)
    if p[i] > 0.5 
      @set! p[i] -= 1
    elseif p[i] < -0.5
      @set! p[i] += 1
    end
  end
  return cell*p + xref
end

"""

```
translation_image(x::SVector{N,T},unit_cell,indices) where {N,T}
```

Translate vector `x` according to the `unit_cell` lattice vectors and the `indices`
provided.

"""
translation_image(x::SVector{N,T},unit_cell,indices) where {T,N} =
  x + unit_cell*SVector{N,T}(ntuple(i -> indices[i],N))

"""

```
neighbour_cells(box::Box{N,T,M}) where {N,M}
```

Function that returns the iterator of the cartesian indices of the cells that must be 
evaluated (forward, i. e. to avoid repeated interactions) 
if the cells have sides of length `box.cell_size`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells(box::Box{3,T,9}) where T
  @unpack lcell = box
  nb = Iterators.flatten((
    CartesianIndices((1:lcell,-lcell:lcell,-lcell:lcell)),
    CartesianIndices((0:0,1:lcell,-lcell:lcell)),
    CartesianIndices((0:0,0:0,1:lcell))
  ))
  return nb
end
function neighbour_cells(box::Box{2,T,4}) where T
  @unpack lcell = box
  nb = Iterators.flatten((
    CartesianIndices((1:lcell,-lcell:lcell)),
    CartesianIndices((0:0,1:lcell))
  ))
  return nb
end

"""

```
neighbour_cells_all(box::Box{N,T,M}) where {N,M}
```

Function that returns the iterator of the cartesian indices of all neighbouring
cells of a cell if the cells have sides of `box.cell_size`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells_all(box::Box{3,T,9}) where T  
  @unpack lcell = box
  return CartesianIndices((-lcell:lcell,-lcell:lcell,-lcell:lcell))
end
function neighbour_cells_all(box::Box{2,T,4}) where T  
  @unpack lcell = box
  return CartesianIndices((-lcell:lcell,-lcell:lcell))
end

"""

```
particle_cell(x::SVector{N,T}, box::Box) where {N,T}
```

Returns the coordinates of the computing cell to which a particle belongs, given its coordinates
and the cell_size. 

"""
@inline particle_cell(x::SVector{N,T}, box::Box) where {N,T} =
  CartesianIndex(ntuple(i -> ceil(Int,x[i]/box.cell_size) + box.lcell, N))

"""

```
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}
```

Given the linear index of the cell in the cell list, returns the cartesian indices 
of the cell (for arbitrary dimension N).

"""
@inline cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} = 
  CartesianIndices(ntuple(i -> nc[i],N))[i1D]

"""

```
cell_linear_index(nc::SVector{N,Int}, indices) where N
```

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

"""
@inline cell_linear_index(nc::SVector{N,Int}, indices) where N =
  LinearIndices(ntuple(i -> nc[i],N))[ntuple(i->indices[i],N)...]




