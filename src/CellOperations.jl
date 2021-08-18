"""

```
wrap_cell_fraction(x,unit_cell_matrix)
```

Obtaint the coordinates of `x` as a fraction of unit cell vectors, first
positive cell. `x` is a vector of dimension `N` and `cell` a matrix of 
dimension `NxN`

### Example

```jldoctest
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_cell_fraction(x,unit_cell_matrix)
2-element Vector{Float64}:
 0.5
 0.3
```

"""
@inline function wrap_cell_fraction(x,unit_cell_matrix)
    p = mod.(unit_cell_matrix\x,1)
    p = @. ifelse(p < 0, p + 1, p)
    return p
end

"""

```
wrap_to_first(x,unit_cell_matrix)
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates. 

### Example

```
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_to_first(x,unit_cell_matrix)
2-element Vector{Float64}:
 5.0
 3.0000000000000004
 ```

"""
@inline function wrap_to_first(x,unit_cell_matrix)
    p = wrap_cell_fraction(x,unit_cell_matrix)
    return unit_cell_matrix*p
end

"""

```
wrap_to_first(x,box::Box)
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates, given the `Box` structure.

"""
@inline wrap_to_first(x,box::Box) = wrap_to_first(x,box.unit_cell.matrix)

"""

```
wrap_to_first(x,box::Box{OrthorhombicCell,N,T}) where {N,T}
```

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates, given an Orthorhombic cell. 
This is slightly cheaper than for general cells.  

"""
@inline function wrap_to_first(x,box::Box{OrthorhombicCell,N,T}) where {N,T}
    sides = SVector{N,T}(ntuple(i->box.unit_cell.matrix[i,i],N))
    x = mod.(x,sides)
    x = @. ifelse(x < 0, x + sides, x)
    return x
end

"""

```
wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

"""
@inline function wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    x_f = wrap_cell_fraction(x,unit_cell_matrix)
    xref_f = wrap_cell_fraction(xref,unit_cell_matrix)
    xw = wrap_relative_to(x_f,xref_f,SVector{N,T}(ntuple(i->1,N)))
    return unit_cell_matrix * (xw - xref_f) + xref
end

"""

```
wrap_relative_to(x,xref,box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`,
given a general `Box` structure.

"""
@inline wrap_relative_to(x,xref,box::Box{UnitCellType,N,T}) where {UnitCellType,N,T} =
    wrap_relative_to(x,xref,box.unit_cell.matrix)

"""

```
wrap_relative_to(x,xref,box::Box{OrthorhombicCell,N,T}) where {N,T}
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`,
given an Orthorhombic cell. This is slightly cheaper than for general cells.

"""
@inline function wrap_relative_to(x,xref,box::Box{OrthorhombicCell,N,T}) where {N,T}
    sides = SVector{N,T}(ntuple(i->box.unit_cell.matrix[i,i],N))
    return wrap_relative_to(x,xref,sides)
end

"""

```
wrap_relative_to(x,xref,sides::AbstractVector)
```

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`,
for an Orthorhombic cell of which only the side lengths are provided.

"""
@inline function wrap_relative_to(x,xref,sides::AbstractVector)
    xw = mod.(x-xref,sides)
    xw = _wrap_single_coordinate.(xw,sides)
    return xw + xref
end

#
# Wrap a single coordinate
#
@inline function _wrap_single_coordinate(x,s)
    if x >= s/2
        x = x - s
    elseif x < -s/2
        x = x + s
    end
    return x
end


"""

```
translation_image(x::SVector{N,T},unit_cell_matrix,indices) where {N,T}
```

Translate vector `x` according to the `unit_cell_matrix` lattice vectors and the `indices`
provided.

"""
@inline translation_image(x::SVector{N,T},unit_cell_matrix,indices) where {N,T} =
    x + unit_cell_matrix*SVector{N,T}(ntuple(i -> indices[i],N))

"""

```
neighbour_cells_forward(box::Box{UnitCellType,N}) where UnitCellType 
```

Function that returns the iterator of the cartesian indices of the cells that must be 
evaluated (forward, i. e. to avoid repeated interactions) 
if the cells have sides of length `box.cell_size`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells_forward(box::Box{UnitCellType,3}) where UnitCellType 
    @unpack lcell = box
    nb = Iterators.flatten((
        CartesianIndices((1:lcell,-lcell:lcell,-lcell:lcell)),
        CartesianIndices((0:0,1:lcell,-lcell:lcell)),
        CartesianIndices((0:0,0:0,1:lcell))
    ))
    return nb
end

function neighbour_cells_forward(box::Box{UnitCellType,2}) where UnitCellType 
    @unpack lcell = box
    nb = Iterators.flatten((
        CartesianIndices((1:lcell,-lcell:lcell)),
        CartesianIndices((0:0,1:lcell))
    ))
    return nb
end

"""

```
neighbour_cells(box::Box{UnitCellType,N}) where N
```

Function that returns the iterator of the cartesian indices of all neighbouring
cells of a cell if the cells have sides of `box.cell_size`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

"""
function neighbour_cells(box::Box{UnitCellType,3}) where UnitCellType
    @unpack lcell = box
    return CartesianIndices((-lcell:lcell,-lcell:lcell,-lcell:lcell))
end

function neighbour_cells(box::Box{UnitCellType,2}) where UnitCellType  
    @unpack lcell = box
    return CartesianIndices((-lcell:lcell,-lcell:lcell))
end

"""

```
particle_cell(x::SVector{N,T}, box::Box) where {N,T}
```

Returns the coordinates of the *computing cell* to which a particle belongs, given its coordinates
and the `cell_size` vector. The computing box is always Orthorhombic, and the first
computing box with positive coordinates has indexes `Box.lcell + 1`.

"""
@inline particle_cell(x::SVector{N}, box::Box) where N =
    CartesianIndex(ntuple(i -> floor(Int,x[i]/box.cell_size[i]) + box.lcell + 1, N))

"""

```
cell_center(c::CartesianIndex{N},box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}
```

Computes the geometric center of a computing cell, to be used in the projection
of points. Returns a `SVector{N,T}`

"""
@inline cell_center(c::CartesianIndex{N},box::Box{UnitCellType,N,T}) where {UnitCellType,N,T} =
    SVector{N,T}(ntuple(i -> box.cell_size[i]*(c[i] - 0.5 - box.lcell), N))


"""

```
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}
```

Given the linear index of the cell in the cell list, returns the cartesian indices 
of the cell (for arbitrary dimension N).

"""
@inline cell_cartesian_indices(nc::SVector{N,Int}, i1D) where N = 
    CartesianIndices(ntuple(i -> nc[i],N))[i1D]

"""

```
cell_linear_index(nc::SVector{N,Int}, indices) where N
```

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

"""
@inline cell_linear_index(nc::SVector{N,Int}, indices) where N =
    LinearIndices(ntuple(i -> nc[i],N))[ntuple(i->indices[i],N)...]

"""

```
out_of_bounding_box(x::SVector{N},box::Box) where N
```

Function that evaluates if a particle is outside the computing bounding box,
defined by the maximum and minimum unit cell coordinates. 

"""
function out_of_bounding_box(x::SVector{N},box::Box) where N
    @unpack cutoff, unit_cell_max = box
    for i in 1:N
        (x[i] < -cutoff) && return true
        (x[i] >= unit_cell_max[i]+cutoff) && return true
    end
    return false
end
out_of_bounding_box(p::ParticleWithIndex,box::Box) =
    out_of_bounding_box(p.coordinates,box)

"""

```
replicate_particle!(ip,p::SVector{N},box,cl) where N
```

Replicates the particle as many times as necessary to fill the computing box.

"""
function replicate_particle!(ip,p::SVector{N},box,cl) where N
    itr = Iterators.product(ntuple(i->box.ranges[i],N)...)
    for indexes in itr
        (count(isequal(0),indexes) == N) && continue
        x = translation_image(p,box.unit_cell.matrix,indexes)
        if ! out_of_bounding_box(x,box)
            cl = add_particle_to_celllist!(ip,x,box,cl;real_particle=false) 
        end
    end 
    return cl
end

"""

```
check_unit_cell(box::Box)
```

Checks if the unit cell satisfies the conditions for using the minimum-image
convention, with the same criteria used gy Gromacs. 

"""
check_unit_cell(box::Box) = check_unit_cell(box.unit_cell.matrix,box.cutoff)

function check_unit_cell(unit_cell_matrix::SMatrix{3},cutoff;printerr=true)
    a = @view(unit_cell_matrix[:,1])
    b = @view(unit_cell_matrix[:,2])
    c = @view(unit_cell_matrix[:,3])
    check = true

    if size(unit_cell_matrix) != (3,3) 
        printerr && println("UNIT CELL CHECK FAILED: unit cell matrix must have dimenions (3,3).")
        check = false
    end

    if count(el -> el < 0, unit_cell_matrix) != 0
         printerr && println("UNIT CELL CHECK FAILED: unit cell matrix components be strictly positive.")
         check = false
    end

    bc = cross(b,c)
    bc = bc / norm(bc)
    aproj = dot(a,bc) 

    ab = cross(a,b)
    ab = ab / norm(ab)
    cproj = dot(c,ab) 

    ca = cross(c,a)
    ca = ca / norm(ca)
    bproj = dot(b,ca) 

    if (aproj < 2*cutoff) || (bproj < 2*cutoff) || (cproj < 2*cutoff)
        printerr && println("UNIT CELL CHECK FAILED: distance between cell planes too small relative to cutoff.")
        check = false
    end

    return check
end

function check_unit_cell(unit_cell_matrix::SMatrix{2},cutoff;printerr=true)
    a = @view(unit_cell_matrix[:,1])
    b = @view(unit_cell_matrix[:,2])
    check = true

    if size(unit_cell_matrix) != (2,2) 
        printerr && println("UNIT CELL CHECK FAILED: unit cell matrix must have dimenions (2,2).")
        check = false
    end

    if count(el -> el < 0, unit_cell_matrix) != 0
         printerr && println("UNIT CELL CHECK FAILED: unit cell matrix components must be strictly positive.")
         check = false
    end

    i = a / norm(a)
    bproj = sqrt(norm_sqr(b) - dot(b,i)^2)

    j = b / norm(b)
    aproj = sqrt(norm_sqr(a) - dot(a,j)^2)

    if (aproj < 2*cutoff) || (bproj < 2*cutoff)
         printerr && println("UNIT CELL CHECK FAILED: distance between cell planes too small relative to cutoff.")
         check = false
    end

    return check
end

"""

Structure that contains the maximum lengths on each direction,
to dispatch on the construction of boxes without periodic boundary
conditions.

"""
struct Limits{T<:AbstractVector} 
    limits::T
end

"""

```
limits(x::AbstractVector{<:AbstractVector})
```

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`, 
to be used to set a box without effective periodic boundary conditions.

"""
function limits(x::AbstractVector{<:AbstractVector})
    xmin = similar(x[1])
    xmax = similar(x[1])
    xmin .= +Inf
    xmax .= -Inf
    for v in x
       @. xmin = min(xmin,v)       
       @. xmax = max(xmax,v)       
    end
    return Limits(xmax .- xmin)
end
"""

```
limits(x::T,y::T) where T<:AbstractVector{<:AbstractVector})
```

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`
and `y`, to used to set a box without effective periodic boundary conditions.

"""
function limits(x::T,y::T) where T<:AbstractVector{<:AbstractVector}
    xlims = limits(x)
    ylims = limits(y)
    return Limits(max.(xlims.limits,ylims.limits))
end

#
# The following functions are not being used anymore in CellListMap, but may
# be useful if this is split in into a sepparate package for playing with
# periodic systems.
#
"""

```
ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix::SMatrix{3,3,T}) where T
```

Function that sets which is the range of periodic images necessary to fill
the computing box, in 3D.

"""
function ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix::SMatrix{3,3,T}) where T
    V = SVector{3,T}
    cmin = -lcell*cell_size
    cmax = (nc .- lcell) .* cell_size
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
        SVector{3,Int}(10^6,10^6,10^6), #min
        SVector{3,Int}(-1,-1,-1),       #max
        unit_cell_matrix,
        cell_vertices
    )
    ranges = SVector{3,UnitRange{Int}}(
        r_min[1]:r_max[1],
        r_min[2]:r_max[2],
        r_min[3]:r_max[3]
    )
    return ranges
end

"""

```
ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix::SMatrix{3,3,T}) where T
```

Function that sets which is the range of periodic images necessary to fill
the computing box, in 2D.

"""
function ranges_of_replicas(cell_size, lcell, nc,unit_cell_matrix::SMatrix{2,2,T}) where T
    V = SVector{2,T}
    cmin = -lcell*cell_size
    cmax = (nc .- lcell) .* cell_size
    cell_vertices = SVector{4,V}( 
        V(   cmin[1],   cmin[2] ), 
        V(   cmin[1],   cmax[2] ),
        V(   cmax[1],   cmin[2] ), 
        V(   cmax[1],   cmax[2] ), 
    )
    r_min, r_max = _ranges_of_replicas(
        SVector{2,Int}(10^6,10^6), #min
        SVector{2,Int}(-1,-1),     #max
        unit_cell_matrix,
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

