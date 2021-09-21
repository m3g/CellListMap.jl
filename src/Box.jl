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
Base.@kwdef struct Box{UnitCellType,N,T,M}
    unit_cell::UnitCell{UnitCellType,N,T,M}
    lcell::Int
    nc::SVector{N,Int}
    cutoff::T
    cutoff_sq::T
    ranges::SVector{N,UnitRange{Int}}
    cell_size::SVector{N,T}
    unit_cell_max::SVector{N,T}
end

"""

```
Box(
  unit_cell_matrix::AbstractMatrix, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  UnitCellType=TriclinicCell
)
```

Construct box structure given the cell matrix of lattice vectors. This 
constructor will always return a `TriclinicCell` box type, unless the
`UnitCellType` parameter is set manually to `OrthorhombicCell`

### Example
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
    ::Type{T},
    lcell::Int,
    ::Type{UnitCellType}
) where {T,UnitCellType}

    s = size(unit_cell_matrix)
    unit_cell_matrix = SMatrix{s[1],s[2],T,s[1]*s[2]}(unit_cell_matrix)

    @assert lcell >= 1 "lcell must be greater or equal to 1"

    N = size(unit_cell_matrix)[1]
    @assert N == size(unit_cell_matrix)[2] "Unit cell matrix must be square."
    @assert check_unit_cell(unit_cell_matrix,cutoff) " Unit cell matrix does not satisfy required conditions."

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
    return Box{UnitCellType,N,T,N*N}(
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
    T::DataType=Float64,
    lcell::Int=1,
    UnitCellType=TriclinicCell
) = Box(unit_cell_matrix,cutoff,T,lcell,UnitCellType)

function Base.show(io::IO,::MIME"text/plain",box::Box)
    println(io,typeof(box))
    println(io,"  unit cell matrix: ", box.unit_cell.matrix) 
    println(io,"  cutoff: ", box.cutoff)
    println(io,"  number of computing cells on each dimension: ",box.nc)
    println(io,"  computing cell sizes: ", box.cell_size, " (lcell: ",box.lcell,")")
    print(io,"  Total number of cells: ", prod(box.nc))
end

"""

```
Box(
  sides::AbstractVector, 
  cutoff, 
  T::DataType, 
  lcell::Int=1,
  UnitCellType=OrthorhombicCell
)
```

For orthorhombic unit cells, `Box` can be initialized with a vector of the 
length of each side. 

### Example
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
    ::Type{T},
    lcell::Int,
    ::Type{UnitCellType}
) where {T,UnitCellType}
    N = length(sides)
    cart_idxs = CartesianIndices((1:N,1:N))
    # Build unit cell matrix from lengths
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
    return Box(
        unit_cell_matrix,
        cutoff,
        T,
        lcell,
        UnitCellType
    ) 
end
Box(
    sides::AbstractVector,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1,
    UnitCellType=OrthorhombicCell
) = Box(sides,cutoff,T,lcell,UnitCellType)

"""

```
Box(
    limits::Limits,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1
)
```

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus the cutoff. This is used to 
emulate pairwise interactions in non-periodic boxes. The output box is always an `Orthorhombic`
cell.

### Examples

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x),10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [109.99633932875878 0.0 0.0; 0.0 109.99780283179763 0.0; 0.0 0.0 109.99587254766517]
  cutoff: 10.0
  number of computing cells on each dimension: [12, 12, 12]
  computing cell sizes: [10.999633932875877, 10.999780283179764, 10.999587254766517] (lcell: 1)
  Total number of cells: 1728

julia> y = [ [150,150,50] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x,y),10)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [159.99787690924168 0.0 0.0; 0.0 159.98878289444897 0.0; 0.0 0.0 109.99587254766517]
  cutoff: 10.0
  number of computing cells on each dimension: [18, 17, 12]
  computing cell sizes: [10.666525127282778, 10.665918859629931, 10.999587254766517] (lcell: 1)
  Total number of cells: 3672

```

"""
Box(
    limits::Limits,
    cutoff;
    T::DataType=Float64,
    lcell::Int=1
) = Box(limits.limits .+ cutoff,cutoff, T,lcell, OrthorhombicCell) 
