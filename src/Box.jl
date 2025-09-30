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

#=

$(TYPEDEF)

# Extended help

$(TYPEDFIELDS)

Structure that contains the maximum lengths on each direction,
to dispatch on the construction of boxes without periodic boundary
conditions.

=#
struct Limits{N,T}
    limits::SVector{N,T}
end

# Set cell size from Limits, when no periodic boundary conditions are 
# used. Adding a fraction of the cutoff to the result avoids 
# the condition of 2*cutoff and the unit cell check check fail
_sides_from_limits(unitcell, cutoff::T) where {T} = unitcell.limits .+ (T(2.1) * cutoff)

#=

$(TYPEDEF)

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

=#
Base.@kwdef struct Box{UnitCellType,N,T,TSQ,M,TR}
    input_unit_cell::UnitCell{UnitCellType,N,T,M}
    aligned_unit_cell::UnitCell{UnitCellType,N,T,M}
    rotation::SMatrix{N,N,TR,M}
    inv_rotation::SMatrix{N,N,TR,M}
    lcell::Int
    nc::SVector{N,Int}
    cutoff::T
    cutoff_sqr::TSQ
    computing_box::Tuple{SVector{N,T},SVector{N,T}}
    cell_size::SVector{N,T}
end

"""
    unitcelltype(::Box{T}) where T = T

Returns the type of a unitcell from the `Box` structure.

## Example

```jldoctest
julia> using CellListMap

julia> box = Box([1,1,1], 0.1);

julia> unitcelltype(box)
OrthorhombicCell

julia> box = Box([1 0 0; 0 1 0; 0 0 1], 0.1);

julia> unitcelltype(box)
TriclinicCell
```

"""
unitcelltype(::Box{T}) where {T} = T

@testitem "unitcelltype" begin
    using CellListMap
    @test unitcelltype(Box([1, 1, 1], 0.1)) == OrthorhombicCell
    @test unitcelltype(Box([1 0 0; 0 1 0; 0 0 1], 0.1)) == TriclinicCell
    x = rand(3, 100)
    @test unitcelltype(Box(limits(x), 0.1)) == NonPeriodicCell
end

#=
    _promote_types(cell,cutoff)

# Extended help

Promotes the types of the unit cell matrix (or sides) and cutoff to floats if one or both were input as integers. 

=#
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
    Box(unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int=1, UnitCellType=TriclinicCell)

Construct box structure given the cell matrix of lattice vectors. This 
constructor will always return a `TriclinicCell` box type, unless the
`UnitCellType` parameter is set manually to `OrthorhombicCell`

## Examples

Building a box with a triclinic unit cell matrix:

```jldoctest
julia> using CellListMap 

julia> unit_cell = [ 100   50    0 
                       0  120    0
                       0    0  130 ];

julia> box = Box(unit_cell, 10.0)
Box{TriclinicCell, 3}
  unit cell matrix = [ 100.0 0.0 0.0; 50.0 120.0 0.0; 0.0 0.0 130.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [20, 13, 16]
  computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells = 4160

```

Building a box with a orthorhombic unit cell matrix, from a square matrix:

```jldoctest
julia> using CellListMap

julia> unit_cell = [ 100 0 0; 0 120 0; 0 0 150 ]; # cell is orthorhombic

julia> box = Box(unit_cell, 10.0, UnitCellType=OrthorhombicCell) # forcing OrthorhombicCell
Box{OrthorhombicCell, 3}
  unit cell matrix = [ 100.0 0.0 0.0; 0.0 120.0 0.0; 0.0 0.0 150.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [13, 15, 18]
  computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells = 3510

```

"""
function Box(input_unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int, ::Type{UnitCellType}) where {UnitCellType}
    input_unit_cell_matrix, cutoff = _promote_types(input_unit_cell_matrix, cutoff)
    T = eltype(input_unit_cell_matrix)

    s = size(input_unit_cell_matrix)
    input_unit_cell_matrix = SMatrix{s[1],s[2],T,s[1] * s[2]}(input_unit_cell_matrix)

    lcell >= 1 || throw(ArgumentError("lcell must be greater or equal to 1"))
    N = size(input_unit_cell_matrix)[1]
    N == size(input_unit_cell_matrix)[2] || throw(ArgumentError("Unit cell matrix must be square."))

    input_unit_cell = UnitCell{UnitCellType,N,T,N * N}(input_unit_cell_matrix)

    return _construct_box(input_unit_cell, lcell, cutoff)
end

Box(unit_cell_matrix::AbstractMatrix, cutoff; lcell::Int=1, UnitCellType=TriclinicCell) =
    Box(unit_cell_matrix, cutoff, lcell, UnitCellType)

# Number of computing cells: for Orthorhombic cells we adjust the cell size such that
# the system has dimensions multiple of cell_size, such that we can use the forward-cell
# method of running over computing cells without complications associated to boundaries
# having fractional cells.
function _compute_nc_and_cell_size(::Type{<:OrthorhombicCellType}, xmin, xmax, cutoff, lcell)
    _nc = floor.(Int, (xmax .- xmin) / (cutoff / lcell))
    cell_size = (xmax .- xmin) ./ _nc
    nc = _nc .+ 2 * lcell .+ 1
    return nc, cell_size
end
function _compute_nc_and_cell_size(::Type{TriclinicCell}, xmin::SVector{N,T}, xmax::SVector{N,T}, cutoff, lcell) where {N,T}
    _nc = ceil.(Int, (xmax .- xmin) / (cutoff / lcell))
    cell_size = SVector{N,T}(ntuple(_ -> cutoff/lcell, Val(N)))
    nc = _nc .+ 2 * lcell .+ 1
    return nc, cell_size
end

# Triclinic cells are aligned such that the largest axis is aligned to x
_align_cell(::Type{TriclinicCell}, m) = align_cell(m)
_align_cell(::Type{<:OrthorhombicCellType}, m) = m, _identity_smatrix(m)

# Creates an identity SMatrix without units
function _identity_smatrix(::SMatrix{N,N,T}) where {N,T}
    oneT = one(T) * inv(one(T)) # remove units
    I = zeros(MMatrix{N,N,typeof(oneT),N*N})
    for i in 1:N
        I[i,i] = oneT
    end
    return SMatrix(I)
end

#
# This function construct the box once the concrete type of the unit cell was obtained.
# This function is useful to perform non-allocating box updates.
#
function _construct_box(input_unit_cell::UnitCell{UnitCellType,N,T}, lcell, cutoff) where {UnitCellType,N,T}

    aligned_unit_cell_matrix, rotation = _align_cell(UnitCellType, input_unit_cell.matrix)
    check_unit_cell(aligned_unit_cell_matrix, cutoff) || throw(ArgumentError(" Unit cell matrix does not satisfy required conditions."))

    aligned_unit_cell = UnitCell{UnitCellType,N,T,N * N}(aligned_unit_cell_matrix)

    # The computing limits are the minimum and maximum coordinates where all particles must be found,
    # including images away from the primitive cell but within the cutoff
    xmin, xmax = cell_limits(aligned_unit_cell.matrix)

    nc, cell_size = _compute_nc_and_cell_size(UnitCellType, xmin, xmax, cutoff, lcell) 
    computing_box = (xmin .- lcell * cell_size, xmax .+ lcell * cell_size)

    # Carry on the squared cutoff, to avoid repeated computation at hot inner loop
    cutoff_sqr = cutoff^2

    return Box{UnitCellType,N,T,typeof(cutoff_sqr),N * N,eltype(rotation)}(
        input_unit_cell=input_unit_cell,
        aligned_unit_cell=aligned_unit_cell,
        rotation=rotation,
        inv_rotation=inv(rotation),
        lcell=lcell,
        nc=nc,
        cutoff=cutoff,
        cutoff_sqr=cutoff_sqr,
        computing_box=computing_box,
        cell_size=cell_size
    )
end

function Base.show(io::IO, ::MIME"text/plain", box::Box{UnitCellType,N}) where {UnitCellType,N}
    _println(io, "Box{$UnitCellType, $N}")
    _print(io, "  unit cell matrix = [ ")
    print(io, join(_uround.(box.input_unit_cell.matrix[1:N, 1]), " "))
    for i in 2:N
        print(io, "; ", join(_uround.(box.input_unit_cell.matrix[1:N, i]), " "))
    end
    println(io, " ]")
    _println(io, "  cutoff = ", box.cutoff)
    _println(io, "  number of computing cells on each dimension = ", box.nc)
    _println(io, "  computing cell sizes = [",
        join(_uround.(box.cell_size), ", "), "] (lcell: ", box.lcell, ")"
    )
    _print(io, "  Total number of cells = ", prod(box.nc))
end

#=
    cell_matrix_from_sides(sides::AbstractVector)

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

=#
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
    Box(sides::AbstractVector, cutoff, lcell::Int=1, UnitCellType=OrthorhombicCell)

For orthorhombic unit cells, `Box` can be initialized with a vector of the length of each side. 

## Example
```jldoctest
julia> using CellListMap

julia> box = Box([120,150,100],10)
Box{OrthorhombicCell, 3}
  unit cell matrix = [ 120.0 0.0 0.0; 0.0 150.0 0.0; 0.0 0.0 100.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [15, 18, 13]
  computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells = 3510
```

"""
function Box(sides::AbstractVector, cutoff, lcell::Int, ::Type{UnitCellType}) where {UnitCellType}
    sides, cutoff = _promote_types(sides, cutoff)
    unit_cell_matrix = cell_matrix_from_sides(sides)
    return Box(unit_cell_matrix, cutoff, lcell, UnitCellType)
end
Box(sides::AbstractVector, cutoff; lcell::Int=1, UnitCellType=OrthorhombicCell) = Box(sides, cutoff, lcell, UnitCellType)

"""
    Box(unitcell::Limits, cutoff; lcell::Int=1)

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus twice the cutoff. This is used to 
emulate pairwise interactions in non-periodic boxes. The output box is an `NonPeriodicCell`
box type, which internally is treated as Orthorhombic with boundaries that guarantee that
particles do not see images of each other. 

## Examples

```jldoctest ;filter = r"\\d+" => ""
julia> using CellListMap, PDBTools

julia> x = coor(readPDB(CellListMap.argon_pdb_file));

julia> box = Box(limits(x), 10.0)
Box{NonPeriodicCell, 3}
  unit cell matrix = [ 39.83 0.0 0.0; 0.0 39.96 0.0; 0.0 0.0 39.99 ]
  cutoff = 10.0
  number of computing cells on each dimension = [6, 6, 6]
  computing cell sizes = [13.28, 13.32, 13.33] (lcell: 1)
  Total number of cells = 216

julia> y = 1.2 .* x;

julia> box = Box(limits(x,y),10)
Box{NonPeriodicCell, 3}
  unit cell matrix = [ 43.6 0.0 0.0; 0.0 43.76 0.0; 0.0 0.0 43.79 ]
  cutoff = 10.0
  number of computing cells on each dimension = [7, 7, 7]
  computing cell sizes = [10.9, 10.94, 10.95] (lcell: 1)
  Total number of cells = 343
```

"""
function Box(unitcell::Limits, cutoff::T; lcell::Int=1) where {T}
    sides = _sides_from_limits(unitcell, cutoff)
    return Box(sides, cutoff, lcell, NonPeriodicCell)
end

# Types of input variables that are acceptable as a unitcell, to construct
# the unitcell matrix. 
const InputUnitCellTypes = Union{Nothing,AbstractVector,AbstractMatrix,Limits,Tuple}

#=
    update_box(
        box::Box{UnitCellType,N,T,TSQ,M};
        unitcell::Union{Nothing,AbstractVector{T},AbstractMatrix{T},Limits,Tuple}=nothing,
        cutoff::Union{Nothing,T}=nothing,
        lcell::Union{Nothing,Int}=nothing
    )

Function that returns an updated system box in a type-stable manner, given possible 
variations in the `unitcell`, `cutoff`, or `lcell` parameters. 

=#
function update_box(
    box::Box{UnitCellType,N,T,TSQ,M};
    unitcell::InputUnitCellTypes=nothing,
    cutoff::Union{Nothing,T}=nothing,
    lcell::Union{Nothing,Int}=nothing
) where {UnitCellType,N,T,TSQ,M}
    _lcell = isnothing(lcell) ? box.lcell : lcell
    _cutoff = isnothing(cutoff) ? box.cutoff : cutoff
    _unitcell = if isnothing(unitcell)
        box.input_unit_cell
    elseif unitcell isa AbstractVector
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa Tuple
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa AbstractMatrix
        UnitCell{UnitCellType,N,T,M}(SMatrix{N,N,T,M}(unitcell))
    elseif unitcell isa Limits
        sides = _sides_from_limits(unitcell, _cutoff)
        UnitCell{UnitCellType,N,T,M}(cell_matrix_from_sides(sides))
    end
    return _construct_box(_unitcell, _lcell, _cutoff)
end

@testitem "Update box with Limits" begin
    using CellListMap
    r = [[1.0, 1.0, 1.0]]
    system = InPlaceNeighborList(x=r, cutoff=3.0, parallel=false)
    list = neighborlist!(system)
    update!(system, r)
    @test list == Tuple{Int64,Int64,Float64}[]
    r = [[1.0, 1.0, 1.0], [10.0, 1.0, 1.0], [3.0, 1.0, 1.0]]
    update!(system, r)
    list = neighborlist!(system)
    @test list == [(1, 3, 2.0)]
    r = [[7,10,10], [18,10,10]]
    system = InPlaceNeighborList(x=r, cutoff=3.0, parallel=false)
    list = neighborlist!(system)
    @test list == Tuple{Int64,Int64,Float64}[]
    r = [[10.89658911843461, 3.709237933444153, 10.0], [13.894156281793144, 11.054172259416013, 10.0]]
    update!(system, r)
    list = neighborlist!(system)
    @test list == Tuple{Int64,Int64,Float64}[]
end

@testitem "Stable Box update" setup=[AllocTest] begin
    using CellListMap
    using StaticArrays
    using BenchmarkTools
    using LinearAlgebra: diag
    using .AllocTest: Allocs

    # update with tuples
    box = Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell=(2, 2, 2), cutoff=0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell=(2, 2, 2), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with SVector
    box = Box([1, 1, 1], 0.1)
    a = @ballocated CellListMap.update_box($box; unitcell=SVector(2, 2, 2), cutoff=0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell=(2, 2, 2), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

    # update with Limits
    x = rand(SVector{3,Float64}, 1000)
    box = Box(limits(x), 0.1)
    new_x = rand(SVector{3,Float64}, 1500)
    a = @ballocated CellListMap.update_box($box; unitcell=$(limits(new_x)), cutoff=0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell=limits(new_x), cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test diag(new_box.input_unit_cell.matrix) == limits(new_x).limits .+ 2.1 * 0.2

    # Update with SMatrix
    box = Box([1 0 0; 0 1 0; 0 0 1], 0.1)
    new_matrix = SMatrix{3,3,Float64,9}(2, 0, 0, 0, 2, 0, 0, 0, 2)
    a = @ballocated CellListMap.update_box($box; unitcell=$new_matrix, cutoff=0.2) evals = 1 samples = 1
    @test a == Allocs(0)
    new_box = CellListMap.update_box(box; unitcell=new_matrix, cutoff=0.2)
    @test new_box.cutoff == 0.2
    @test new_box.input_unit_cell.matrix == [2 0 0; 0 2 0; 0 0 2]

end

#=
    neighbor_cells_forward(box::Box{UnitCellType,N}) where UnitCellType 

# Extended help

Function that returns the iterator of the cartesian indices of the cells that must be 
evaluated (forward, i. e. to avoid repeated interactions) 
if the cells have sides of length `box.cell_size`. `N` can be
`2` or `3`, for two- or three-dimensional systems.

=#
function neighbor_cells_forward(box::Box{UnitCellType,3}) where {UnitCellType}
    @unpack lcell = box
    nb = Iterators.flatten((
        CartesianIndices((1:lcell, -lcell:lcell, -lcell:lcell)),
        CartesianIndices((0:0, 1:lcell, -lcell:lcell)),
        CartesianIndices((0:0, 0:0, 1:lcell))
    ))
    return nb
end

function neighbor_cells_forward(box::Box{UnitCellType,2}) where {UnitCellType}
    @unpack lcell = box
    nb = Iterators.flatten((
        CartesianIndices((1:lcell, -lcell:lcell)),
        CartesianIndices((0:0, 1:lcell))
    ))
    return nb
end

#=
    neighbor_cells(box::Box{UnitCellType,N}) where {UnitCellType,N}

# Extended help

Function that returns the iterator of the cartesian indices of all neighboring
cells of a cell where the computing cell index is `box.lcell`.

=#
function neighbor_cells(box::Box{UnitCellType,N}) where {UnitCellType,N}
    @unpack lcell = box
    return Iterators.filter(
        !isequal(CartesianIndex(ntuple(i -> 0, Val(N)))),
        CartesianIndices(ntuple(i -> -lcell:lcell, Val(N)))
    )
end

#=
    current_and_neighbor_cells(box::Box{UnitCellType,N}) where {UnitCellType,N}

# Extended help

Returns an iterator over all neighbor cells, including the center one.

=#
function current_and_neighbor_cells(box::Box{UnitCellType,N}) where {UnitCellType,N}
    @unpack lcell = box
    return CartesianIndices(ntuple(i -> -lcell:lcell, Val(N)))
end

#=
    particle_cell(x::SVector{N,T}, box::Box) where {N,T}

# Extended help

Returns the coordinates of the *computing cell* to which a particle belongs, given its coordinates
and the `cell_size` vector. The computing box is always Orthorhombic, and the first
computing box with positive coordinates has indexes `Box.lcell + 1`.

=#
@inline function particle_cell(x::SVector{N}, box::Box) where {N}
    CartesianIndex(
        ntuple(Val(N)) do i
            xmin = box.computing_box[1][i]
            xi = (x[i] - xmin) / box.cell_size[i]
            index = floor(Int, xi) + 1
            return index
        end
    )
end

#=
    cell_center(c::CartesianIndex{N},box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}

Computes the geometric center of a computing cell, to be used in the projection
of points. Returns a `SVector{N,T}`

=#
@inline function cell_center(c::CartesianIndex{N}, box::Box{UnitCellType,N,T}) where {UnitCellType,N,T}
    SVector{N,T}(
        ntuple(Val(N)) do i
            xmin = box.computing_box[1][i]
            ci = xmin + box.cell_size[i] * c[i] - box.cell_size[i] / 2
            return ci
        end
    )
end

#=
    in_computing_box(x::SVector{N},box::Box) where N

# Extended help

Function that evaluates if a particle is inside the computing bounding box,
defined by the maximum and minimum unit aligned cell coordinates.

=#
function in_computing_box(x::SVector{N}, box::Box) where {N}
    min = box.computing_box[1]
    max = box.computing_box[2]
    inbox = true
    for i in 1:N
        if !(min[i] <= x[i] < max[i])
            inbox = false
        end
    end
    return inbox
end

#=
    replicate_particle!(ip,p::SVector{N},box,cl) where N

# Extended help

Replicates the particle as many times as necessary to fill the computing box.

=#
function replicate_particle!(ip, p::SVector{N}, box, cl) where {N}
    itr = Iterators.product(ntuple(i -> -1:1, Val(N))...)
    for indices in itr
        (count(isequal(0), indices) == N) && continue
        x = translation_image(p, box.aligned_unit_cell.matrix, indices)
        if in_computing_box(x, box)
            cl = add_particle_to_celllist!(ip, x, box, cl; real_particle=false)
        end
    end
    return cl
end

#=
    check_unit_cell(box::Box)

# Extended help

Checks if the unit cell satisfies the conditions for using the minimum-image
convention. 

=#
check_unit_cell(box::Box) = check_unit_cell(box.aligned_unit_cell.matrix, box.cutoff)

function check_unit_cell(unit_cell_matrix::SMatrix{3}, cutoff; printerr=true)
    a = @view(unit_cell_matrix[:, 1])
    b = @view(unit_cell_matrix[:, 2])
    c = @view(unit_cell_matrix[:, 3])
    check = true

    if size(unit_cell_matrix) != (3, 3)
        printerr && println("UNIT CELL CHECK FAILED: unit cell matrix must have dimensions (3,3).")
        check = false
    end

    bc = cross(b, c)
    bc = bc / norm(bc)
    aproj = dot(a, bc)

    ab = cross(a, b)
    ab = ab / norm(ab)
    cproj = dot(c, ab)

    ca = cross(c, a)
    ca = ca / norm(ca)
    bproj = dot(b, ca)

    if (aproj <= 2 * cutoff) || (bproj <= 2 * cutoff) || (cproj <= 2 * cutoff)
        printerr && println("UNIT CELL CHECK FAILED: distance between cell planes too small relative to cutoff.")
        check = false
    end

    return check
end

function check_unit_cell(unit_cell_matrix::SMatrix{2}, cutoff; printerr=true)
    a = @view(unit_cell_matrix[:, 1])
    b = @view(unit_cell_matrix[:, 2])
    check = true

    if size(unit_cell_matrix) != (2, 2)
        printerr && println("UNIT CELL CHECK FAILED: unit cell matrix must have dimensions (2,2).")
        check = false
    end

    i = a / norm(a)
    bproj = sqrt(norm_sqr(b) - dot(b, i)^2)

    j = b / norm(b)
    aproj = sqrt(norm_sqr(a) - dot(a, j)^2)

    if (aproj <= 2 * cutoff) || (bproj <= 2 * cutoff)
        printerr && println("UNIT CELL CHECK FAILED: distance between cell planes too small relative to cutoff.")
        check = false
    end

    return check
end


