#
# This difference is important because for Orthorhombic cells it is possible to run over 
# only half of the cells, and wrapping coordinates in Orthorhombic cells is slightly cheaper. 
#
struct TriclinicCell end
struct OrthorhombicCell end
struct NonPeriodicCell end

const OrthorhombicCellType = Union{OrthorhombicCell,NonPeriodicCell}
const PeriodicCellType = Union{OrthorhombicCell,TriclinicCell}

# Wrapper for the unitcell matrix, to be able to dispatch on the unitcell type that contains it
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
struct Limits{N,T}
    limits::SVector{N,T}
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
function _promote_types(cell, cutoff)
    input_type = promote_type(eltype(cell), typeof(cutoff))
    float_type = input_type == Int ? Float64 : input_type
    return float_type.(cell), float_type(cutoff)
end

@testitem "promote types" begin
    import CellListMap: _promote_types
    cell, cutoff = _promote_types([1, 1, 1], 0.1)
    @test (eltype(cell), typeof(cutoff)) == (Float64, Float64)
    cell, cutoff = _promote_types([1, 1, 1], 0.1f0)
    @test (eltype(cell), typeof(cutoff)) == (Float32, Float32)
    cell, cutoff = _promote_types([1.0, 1, 1], 0.1f0)
    @test (eltype(cell), typeof(cutoff)) == (Float64, Float64)
    cell, cutoff = _promote_types([1.0f0, 1, 1], 0.1f0)
    @test (eltype(cell), typeof(cutoff)) == (Float32, Float32)
    cell, cutoff = _promote_types([10, 10, 10.0], 1)
    @test (eltype(cell), typeof(cutoff)) == (Float64, Float64)
    cell, cutoff = _promote_types([10.0f0, 10, 10], 1)
    @test (eltype(cell), typeof(cutoff)) == (Float32, Float32)
    cell, cutoff = _promote_types([10, 10, 10], 1)
    @test (eltype(cell), typeof(cutoff)) == (Float64, Float64)
end

"""

```
Box(unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int=1, UnitCellType=TriclinicCell)
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
function Box(unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int, ::Type{UnitCellType}) where {UnitCellType}
    unit_cell_matrix, cutoff = _promote_types(unit_cell_matrix, cutoff)
    T = eltype(unit_cell_matrix)

    s = size(unit_cell_matrix)
    unit_cell_matrix = SMatrix{s[1],s[2],T,s[1] * s[2]}(unit_cell_matrix)

    lcell >= 1 || throw(ArgumentError("lcell must be greater or equal to 1"))
    N = size(unit_cell_matrix)[1]
    N == size(unit_cell_matrix)[2] || throw(ArgumentError("Unit cell matrix must be square."))
    check_unit_cell(unit_cell_matrix, cutoff) || throw(ArgumentError(" Unit cell matrix does not satisfy required conditions."))

    unit_cell = UnitCell{UnitCellType,N,T,N * N}(unit_cell_matrix)
    unit_cell_max, nc, cell_size, ranges = _set_unitcell_ranges(unit_cell, lcell, cutoff)

    return Box{UnitCellType,N,T,typeof(cutoff^2),N * N}(unit_cell, lcell, nc, cutoff, cutoff^2, ranges, cell_size, unit_cell_max)
end
Box(unit_cell_matrix::AbstractMatrix, cutoff; lcell::Int=1, UnitCellType=TriclinicCell) = 
    Box(unit_cell_matrix, cutoff, lcell, UnitCellType)

#
# The following functions define some internal parameters required for the calculations. In particular,
# the maximum lengths of the unit cell in each dimension, the number of cells, the cell size, and the
# ranges of vicinal cells where it is necessary to create imaginary particles to cope with periodic
# boundary conditions. 
#
function _set_unitcell_ranges(unitcell::UnitCell{TriclinicCell,N,T}, lcell, cutoff) where {N,T}
    unit_cell_max = sum(unitcell.matrix[:, i] for i in 1:N)
    cell_size = SVector{N,T}(ntuple(i -> cutoff / lcell, N))
    nc = ceil.(Int, (unit_cell_max .+ 2 * cutoff) ./ cell_size)
    ranges = SVector{N,UnitRange{Int}}(ntuple(i -> -1:1, N))
    return unit_cell_max, nc, cell_size, ranges
end
# For Ortorhombic cells, the cell size must be a multiple of the cell size. In some pathological cases 
# this may be bad for performance, because we are increasing the effective cell cutoff.
function _set_unitcell_ranges(unitcell::UnitCell{<:OrthorhombicCellType,N,T}, lcell, cutoff) where {N,T}
    unit_cell_max = sum(unitcell.matrix[:, i] for i in 1:N)
    nc = floor.(Int, lcell * unit_cell_max / cutoff)
    cell_size = unit_cell_max ./ nc
    nc = nc .+ 2 * lcell
    ranges = SVector{N,UnitRange{Int}}(ntuple(i -> -1:1, N))
    return unit_cell_max, nc, cell_size, ranges
end

function Base.show(io::IO, ::MIME"text/plain", box::Box{UnitCellType,N}) where {UnitCellType,N}
    _println(io, "Box{$UnitCellType, $N}")
    _print(io, "  unit cell matrix = [ ")
    print(io, join(_uround.(box.unit_cell.matrix[1:N, 1]), ", "))
    for i in 2:N
        print(io, "; ", join(_uround.(box.unit_cell.matrix[1:N, i]), ", "))
    end
    println(io, " ]")
    _println(io, "  cutoff = ", box.cutoff)
    _println(io, "  number of computing cells on each dimension = ", box.nc)
    _println(io, "  computing cell sizes = [",
        join(_uround.(box.cell_size), ", "), "] (lcell: ", box.lcell, ")"
    )
    _print(io, "  Total number of cells = ", prod(box.nc))
end

"""

```
cell_matrix_from_sides(sides::AbstractVector)
```

$(INTERNAL)

# Extended help

Function that returns the Orthorhombic unit cell matrix given a sides vector.
This function is type-unstable if the input is not static.

## Example

```
julia> CellListMap.cell_matrix_from_sides([1,1,1])
3×3 SMatrix{3, 3, Int64, 9} with indices SOneTo(3)×SOneTo(3):
 1  0  0
 0  1  0
 0  0  1
```

"""
cell_matrix_from_sides(sides::AbstractVector) = oneunit(eltype(sides)) * diagm(sides ./ oneunit(eltype(sides)))
cell_matrix_from_sides(sides::Tuple) = cell_matrix_from_sides(SVector(sides))

@testitem "cell_matrix_from_sides: units" begin
    using StaticArrays
    using Unitful
    import CellListMap: cell_matrix_from_sides

    sides = SVector(1.0, 1.0, 1.0)u"nm"
    cell_matrix = cell_matrix_from_sides(sides)
    @test eltype(sides) == eltype(cell_matrix) 

    sides = (1.0u"nm", 1.0u"nm", 1.0u"nm")
    cell_matrix = cell_matrix_from_sides(sides)
    @test eltype(sides) == eltype(cell_matrix) 
end

"""

```
Box(sides::AbstractVector, cutoff, lcell::Int=1, UnitCellType=OrthorhombicCell)
```

For orthorhombic unit cells, `Box` can be initialized with a vector of the length of each side. 

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
function Box(sides::AbstractVector, cutoff, lcell::Int, ::Type{UnitCellType}) where {UnitCellType}
    sides, cutoff = _promote_types(sides, cutoff)
    unit_cell_matrix = cell_matrix_from_sides(sides)
    return Box(unit_cell_matrix, cutoff, lcell, UnitCellType)
end
Box(sides::AbstractVector, cutoff; lcell::Int=1, UnitCellType=OrthorhombicCell) = Box(sides, cutoff, lcell, UnitCellType)

"""

```
Box(unitcell::Limits, cutoff; lcell::Int=1)
```

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus twice the cutoff. This is used to 
emulate pairwise interactions in non-periodic boxes. The output box is an `NonPeriodicCell`
box type, which internally is treated as Orthorhombic with boundaries that guarantee that
particles do not see images of each other. 

## Examples

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x),10)
Box{NonPeriodicCell, 3}
  unit cell matrix = [ 110.0, 0.0, 0.0; 0.0, 110.0, 0.0; 0.0, 0.0, 110.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [12, 12, 12]
  computing cell sizes = [11.0, 11.0, 11.0] (lcell: 1)
  Total number of cells = 1728

julia> y = [ [150,150,50] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x,y),10)
Box{NonPeriodicCell, 3}
  unit cell matrix = [ 160.0, 0.0, 0.0; 0.0, 160.0, 0.0; 0.0, 0.0, 110.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [17, 17, 12]
  computing cell sizes = [10.67, 10.67, 11.0] (lcell: 1)
  Total number of cells = 3468
```

"""
function Box(unitcell::Limits, cutoff::T; lcell::Int=1) where {T}
    sides = max.(unitcell.limits .+ cutoff, 2 * cutoff)
    return Box(sides, cutoff, lcell, NonPeriodicCell)
end


# Types of input variables that are acceptable as a unitcell, to construct
# the unitcell matrix. 
const InputUnitCellTypes = Union{Nothing,AbstractVector,AbstractMatrix,Limits,Tuple}

"""
update_box(
    box::Box{UnitCellType,N,T,TSQ,M};
    unitcell::Union{Nothing,AbstractVector{T},AbstractMatrix{T},Limits,Tuple}=nothing,
    cutoff::Union{Nothing,T}=nothing,
    lcell::Union{Nothing,Int}=nothing
)

$(INTERNAL)

Function that returns an updated system box in a type-stable manner, given possible 
variations in the `unitcell`, `cutoff`, or `lcell` parameters. 

"""
function update_box(
    box::Box{UnitCellType,N,T,TSQ,M};
    unitcell::InputUnitCellTypes=nothing,
    cutoff::Union{Nothing,T}=nothing,
    lcell::Union{Nothing,Int}=nothing
) where {UnitCellType,N,T,TSQ,M}
    _lcell = isnothing(lcell) ? box.lcell : lcell
    _cutoff = isnothing(cutoff) ? box.cutoff : cutoff
    _unitcell = if isnothing(unitcell)
        box.unit_cell
    elseif unitcell isa AbstractVector
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa Tuple
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa AbstractMatrix
        UnitCell{UnitCellType,N,T,M}(SMatrix{N,N,T,M}(unitcell))
    elseif unitcell isa Limits
        sides = SVector(max.(unitcell.limits .+ _cutoff, 2 * _cutoff))
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(sides))
    end
    unit_cell_max, nc, cell_size, ranges = _set_unitcell_ranges(_unitcell, _lcell, _cutoff)
    return Box{UnitCellType,N,T,TSQ,M}(
        _unitcell, _lcell, nc, _cutoff, _cutoff^2, ranges, cell_size, unit_cell_max
    )
end

@testitem "Stable Box update" begin
    using CellListMap
    using StaticArrays
    using BenchmarkTools
    using LinearAlgebra: diag

    # update with tuples
    box = Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell=(2, 2, 2), cutoff=0.2) evals=1 samples=1
    @test a == 0
    new_box = CellListMap.update_box(box; unitcell=(2, 2, 2), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with SVector
    box = Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell=SVector(2, 2, 2), cutoff=0.2) evals=1 samples=1
    @test a == 0
    new_box = CellListMap.update_box(box; unitcell=(2, 2, 2), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with Limits
    x = rand(SVector{3,Float64}, 1000)
    box = Box(limits(x), 0.1)
    new_x = rand(SVector{3,Float64}, 1500)
    a = @ballocated CellListMap.update_box($box; unitcell=$(limits(new_x)), cutoff=0.2) evals=1 samples=1
    @test a == 0
    new_box = CellListMap.update_box(box; unitcell=limits(new_x), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test diag(new_box.unit_cell.matrix) == limits(new_x).limits .+ 0.2

    # Update with SMatrix
    box = Box([1 0 0; 0 1 0; 0 0 1], 0.1)
    new_matrix = SMatrix{3,3,Float64,9}(2, 0, 0, 0, 2, 0, 0, 0, 2)
    a = @ballocated CellListMap.update_box($box; unitcell=$new_matrix, cutoff=0.2) evals=1 samples=1
    @test a == 0
    new_box = CellListMap.update_box(box; unitcell=new_matrix, cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

end