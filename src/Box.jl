#
# This difference is important because for Orthorhombic cells it is
# possible to run over only half of the cells, and wrapping coordinates
# in Orthorhombic cells is slightly cheaper. 
#
struct TriclinicCell end
struct OrthorhombicCell end

struct UnitCell{UnitCellType,N,T,M}
    matrix::SMatrix{N,N,T,M}
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that contains the maximum lengths on each direction,
to dispatch on the construction of boxes without periodic boundary
conditions.

"""
struct Limits{T<:AbstractVector} 
    limits::T
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

$(TYPEDFIELDS)

Structure that contains some data required to compute the linked cells. To
be initialized with the box size and cutoff. 

## Examples

```julia-repl
julia> using CellListMap

julia> sides = [250,250,250];

julia> cutoff = 10;

julia> box = Box(sides,cutoff)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [250.0 0.0 0.0; 0.0 250.0 0.0; 0.0 0.0 250.0]
  cutoff: 10.0
  number of computing cells on each dimension: [27, 27, 27]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 19683

```

```julia-repl
julia> box = Box([ 10  0  0 
                    0 10  5
                    0  0 10 ], 1)
Box{TriclinicCell, 3, Float64, 9}
unit cell matrix: [10.0 0.0 0.0; 0.0 10.0 5.0; 0.0 0.0 10.0]
cutoff: 1.0
number of computing cells on each dimension: [12, 17, 12]
computing cell sizes: [1.0, 1.0, 1.0] (lcell: 1)
Total number of cells: 2448

```

"""
Base.@kwdef struct Box{UnitCellType,N,T,TSQ,M}
    unit_cell::UnitCell{UnitCellType,N,T,M}
    lcell::Int
    nc::SVector{N,Int}
    cutoff::T
    cutoff_sq::TSQ
    ranges::SVector{N,UnitRange{Int}}
    cell_size::SVector{N,T}
    unit_cell_max::SVector{N,T}
end

"""

```
_promote_types(cell,cutoff)
```

$(INTERNAL)

# Extended help

Promotes the types of the unit cell matrix (or sides) and cutoff to floats if one or both were input as integers. 

"""
function _promote_types(cell,cutoff)
    if eltype(cell) <: Integer || typeof(cutoff) <: Integer
        if cutoff isa Integer && eltype(cell) <: Integer
            cell = convert.(Float64,cell)
            cutoff = convert(eltype(cell),cutoff)
        elseif eltype(cell) <: Integer
            cell = convert.(typeof(cutoff),cell)
        else
            cutoff = convert(eltype(cell),cutoff)
        end
    end
    return cell, cutoff
end

"""

```
Box(
  unit_cell_matrix::AbstractMatrix, 
  cutoff, 
  lcell::Int=1,
  UnitCellType=TriclinicCell
)
```

Construct box structure given the cell matrix of lattice vectors. This 
constructor will always return a `TriclinicCell` box type, unless the
`UnitCellType` parameter is set manually to `OrthorhombicCell`

## Example
```julia-repl
julia> unit_cell = [ 100   50    0 
                       0  120    0
                       0    0  130 ];

julia> box = Box(unit_cell,10)
Box{TriclinicCell, 3, Float64, 9}
  unit cell matrix: [100.0 50.0 0.0; 0.0 120.0 0.0; 0.0 0.0 130.0]
  cutoff: 10.0
  number of computing cells on each dimension: [17, 14, 15]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 3570

```

"""
function Box(
    unit_cell_matrix::AbstractMatrix, 
    cutoff, 
    lcell::Int,
    ::Type{UnitCellType}
) where {UnitCellType}
    unit_cell_matrix, cutoff = _promote_types(unit_cell_matrix, cutoff)
    T = eltype(unit_cell_matrix)

    s = size(unit_cell_matrix)
    unit_cell_matrix = SMatrix{s[1],s[2],T,s[1]*s[2]}(unit_cell_matrix)

    lcell >= 1 || throw(ArgumentError("lcell must be greater or equal to 1"))
    N = size(unit_cell_matrix)[1]
    N == size(unit_cell_matrix)[2] || throw(ArgumentError("Unit cell matrix must be square."))
    check_unit_cell(unit_cell_matrix,cutoff) || throw(ArgumentError(" Unit cell matrix does not satisfy required conditions."))

    unit_cell = UnitCell{UnitCellType,N,T,N*N}(unit_cell_matrix)
    unit_cell_max = sum(unit_cell_matrix[:,i] for i in 1:N) 

    # To use the advantages of Orthorhombic cells, box size must
    # be  multiple of the cell size. In some pathological cases this
    # may be bad for performance, because we are increasing the 
    # effective cell cutoff
    if UnitCellType == OrthorhombicCell
        nc = floor.(Int,lcell*unit_cell_max/cutoff)
        cell_size = unit_cell_max ./ nc 
        nc = nc .+ 2*lcell 
    # For triclinic cells we have to use the ceil here, because we need
    # to major the number of cells, when the box size is not a multiple of
    # the cell size
    else
        cell_size = SVector{N,T}(ntuple(i->cutoff/lcell,N)...)
        nc = ceil.(Int,(unit_cell_max .+ 2*cutoff) ./ cell_size)
    end

    #ranges = ranges_of_replicas(cell_size, lcell, nc, unit_cell_matrix)
    ranges = SVector{N,UnitRange{Int}}(ntuple(i->-1:1,N))
    return Box{UnitCellType,N,T,typeof(cutoff^2),N*N}(
        unit_cell,
        lcell, 
        nc,
        cutoff,
        cutoff^2,
        ranges,
        cell_size,
        unit_cell_max
    )
end
Box(
    unit_cell_matrix::AbstractMatrix,
    cutoff;
    lcell::Int=1,
    UnitCellType=TriclinicCell
) = Box(unit_cell_matrix,cutoff,lcell,UnitCellType)

function Base.show(io::IO,::MIME"text/plain",box::Box{UnitCellType,N}) where {UnitCellType,N}
    println(io,typeof(box))
    print(io,"  unit cell matrix = [ ") 
    print(io,join(box.unit_cell.matrix[1:N,1],", "))
    for i in 2:N
        print(io,"; ",join(box.unit_cell.matrix[1:N,i],", "))
    end
    println(io," ]")
    println(io,"  cutoff = ", box.cutoff)
    println(io,"  number of computing cells on each dimension = ",box.nc)
    println(io,"  computing cell sizes = [", join(box.cell_size,", "), "] (lcell: ",box.lcell,")")
    print(io,"  Total number of cells = ", prod(box.nc))
end

"""

```
cell_matrix_from_sides(sides::AbstractVector)
```

$(INTERNAL)

# Extended help

Function that returns the Orthorhombic unit cell matrix given a sides vector.

## Example

```
julia> CellListMap.cell_matrix_from_sides([1,1,1])
3×3 SMatrix{3, 3, Int64, 9} with indices SOneTo(3)×SOneTo(3):
 1  0  0
 0  1  0
 0  0  1
```

"""
function cell_matrix_from_sides(sides::AbstractVector)
    T = eltype(sides)
    N = length(sides)
    cart_idxs = CartesianIndices((1:N,1:N))
    unit_cell_matrix = SMatrix{N,N,T,N*N}( 
        ntuple(N*N) do i
            c = cart_idxs[i]
            if c[1] == c[2] 
                return sides[c[1]] 
            else
                return zero(T)
            end
        end
    )
    return unit_cell_matrix
end

"""

```
Box(
  sides::AbstractVector, 
  cutoff, 
  lcell::Int=1,
  UnitCellType=OrthorhombicCell
)
```

For orthorhombic unit cells, `Box` can be initialized with a vector of the 
length of each side. 

## Example
```julia-repl
julia> box = Box([120,150,100],10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [120.0 0.0 0.0; 0.0 150.0 0.0; 0.0 0.0 100.0]
  cutoff: 10.0
  number of computing cells on each dimension: [14, 17, 12]
  computing cell sizes: [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells: 2856

```

"""
function Box(
    sides::AbstractVector, 
    cutoff, 
    lcell::Int,
    ::Type{UnitCellType}
) where {UnitCellType}
    sides, cutoff = _promote_types(sides, cutoff)
    unit_cell_matrix = cell_matrix_from_sides(sides)
    return Box(unit_cell_matrix,cutoff,lcell,UnitCellType) 
end
Box(
    sides::AbstractVector,
    cutoff;
    lcell::Int=1,
    UnitCellType=OrthorhombicCell
) = Box(sides,cutoff,lcell,UnitCellType)

"""

```
Box(
    limits::Limits,
    cutoff;
    lcell::Int=1
)
```

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus twice the cutoff. This is used to 
emulate pairwise interactions in non-periodic boxes. The output box is always an `Orthorhombic`
cell.

## Examples

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x),10)
Box{OrthorhombicCell, 3, Float64, Float64, 9}
  unit cell matrix = [ 119.99907193208746, 0.0, 0.0; 0.0, 119.99968623301143, 0.0; 0.0, 0.0, 119.99539603156498 ]
  cutoff = 10.0
  number of computing cells on each dimension = [13, 13, 13]
  computing cell sizes = [10.909006539280679, 10.90906238481922, 10.908672366505908] (lcell: 1)
  Total number of cells = 2197

julia> y = [ [150,150,50] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x,y),10)
Box{OrthorhombicCell, 3, Float64, Float64, 9}
  unit cell matrix = [ 169.99914503548962, 0.0, 0.0; 0.0, 169.9990736881799, 0.0; 0.0, 0.0, 119.99726063023918 ]
  cutoff = 10.0
  number of computing cells on each dimension = [18, 18, 13]
  computing cell sizes = [10.624946564718101, 10.624942105511243, 10.90884187547629] (lcell: 1)
  Total number of cells = 4212

```

"""
function Box(limits::Limits, cutoff::T; lcell::Int=1) where T
    sides = max.(limits.limits .+ cutoff, 2 * cutoff)
    return Box(sides, cutoff, lcell, OrthorhombicCell) 
end
