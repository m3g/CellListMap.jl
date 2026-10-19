#
# Types of unit cells. Orthorhombic and triclinic cells are treated identically by the
# cell lists (the grid is defined in fractional coordinates), but the type is kept for
# dispatch and for the interpretation of the `unitcell` input (sides vs. matrix).
# Non-periodic systems are represented by an orthorhombic box enclosing all particles,
# in which neighboring cells are never wrapped.
#
struct TriclinicCell end
struct OrthorhombicCell end
struct NonPeriodicCell end
const OrthorhombicCellType = Union{OrthorhombicCell, NonPeriodicCell}
const PeriodicCellType = Union{OrthorhombicCell, TriclinicCell}

# Wrapper for the unitcell matrix, to be able to dispatch on the unitcell type that contains it
struct UnitCell{UnitCellType, N, T, M}
    matrix::SMatrix{N, N, T, M}
end

#=
Structure that contains the maximum lengths on each direction,
to dispatch on the construction of boxes without periodic boundary
conditions.

=#
struct Limits{N, T}
    limits::SVector{N, T}
    origin::SVector{N, T}
end
Limits(limits::SVector{N, T}) where {N, T} = Limits(limits, zeros(SVector{N, T}))

# Set cell size from Limits, when no periodic boundary conditions are
# used. A margin larger than the cutoff is added, such that particles
# are always far from the boundaries of the box.
_sides_from_limits(unitcell, cutoff::T) where {T} = unitcell.limits .+ (210 * cutoff / 100)

#=
Structure that contains the data required to compute the cell lists: the unit cell,
the cutoff, and the grid of computing cells.

The computing cells are defined in fractional coordinates: the unit cell is divided into
`nc[i]` slices along each lattice vector `i`, such that each computing cell is a small
parallelepiped similar to the unit cell. The number of cells is chosen such that the distance
between opposite faces of the computing cells is at least `cutoff / lcell`, and the number
of neighboring cells that must be searched in each direction is `n_search[i]`. There is no
restriction on the size of the unit cell relative to the cutoff.

## Examples

```julia-repl
julia> using CellListMap

julia> box = Box([250, 250, 250], 10)
Box{OrthorhombicCell, 3}
  unit cell matrix = [ 250.0 0.0 0.0; 0.0 250.0 0.0; 0.0 0.0 250.0 ]
  cutoff = 10.0
  number of computing cells on each dimension = [25, 25, 25]
  computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
  Total number of cells = 15625
```

=#
Base.@kwdef struct Box{UnitCellType, N, T, TSQ, M, TI}
    input_unit_cell::UnitCell{UnitCellType, N, T, M}
    # inverse of the unit cell matrix (without units), to compute fractional coordinates
    inv_unit_cell::SMatrix{N, N, TI, M}
    lcell::Int
    # number of computing cells along each lattice vector
    nc::SVector{N, Int}
    # number of neighboring cells to search along each lattice vector
    n_search::SVector{N, Int}
    cutoff::T
    cutoff_sqr::TSQ
    computing_box::Tuple{SVector{N, T}, SVector{N, T}}
    origin::SVector{N, T}
end

#=
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

=#
unitcelltype(::Box{T}) where {T} = T

#=
    _promote_types(cell,cutoff)

# Extended help

Promotes the types of the unit cell matrix (or sides) and cutoff to floats if one or both were input as integers,
returns the minimum float type necessary.

=#
function _promote_types(unitcell, cutoff)
    input_type = if isnothing(unitcell)
        typeof(cutoff)
    elseif isnothing(cutoff)
        eltype(unitcell)
    else
        promote_type(eltype(unitcell), typeof(cutoff))
    end
    float_type = input_type == Int ? Float64 : input_type
    return float_type
end

#=
    Box(unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int=1, UnitCellType=TriclinicCell)

Construct box structure given the cell matrix, where
columns correspond to the lattice vectors. This constructor will always return
a `TriclinicCell` box type, unless the `UnitCellType` parameter is set
manually to `OrthorhombicCell`

=#
function Box(input_unit_cell_matrix::AbstractMatrix, cutoff, lcell::Int, ::Type{UnitCellType}; origin = nothing) where {UnitCellType}
    lcell >= 1 || throw(ArgumentError("lcell must be greater or equal to 1"))
    s = size(input_unit_cell_matrix)
    s[1] == s[2] || throw(ArgumentError("Unit cell matrix must be square."))
    T = _promote_types(input_unit_cell_matrix, cutoff)
    unit_cell_matrix = SMatrix{s[1], s[2], T, s[1] * s[2]}(input_unit_cell_matrix)
    unit_cell = UnitCell{UnitCellType, s[1], T, s[1] * s[2]}(unit_cell_matrix)
    _origin = isnothing(origin) ? zeros(SVector{s[1], T}) : SVector{s[1], T}(origin)
    return _construct_box(unit_cell, lcell, T(cutoff); origin = _origin)
end

Box(unit_cell_matrix::AbstractMatrix, cutoff; lcell::Int = 1, UnitCellType = TriclinicCell) =
    Box(unit_cell_matrix, cutoff, lcell, UnitCellType)

#
# This function construct the box once the concrete type of the unit cell was obtained.
# This function is useful to perform non-allocating box updates.
#
function _construct_box(input_unit_cell::UnitCell{UnitCellType, N, T}, lcell, cutoff; origin::SVector{N, T} = zeros(SVector{N, T})) where {UnitCellType, N, T}
    m = input_unit_cell.matrix
    face_distances = _face_distances(m)
    if !all(>(zero(T)), face_distances)
        throw(ArgumentError("Unit cell matrix is singular: the unit cell has zero volume."))
    end
    # Number of computing cells along each lattice vector, such that the distance between
    # the faces of the cells is at least cutoff / lcell, and number of cells to be searched
    # in each direction to cover the cutoff.
    nc = SVector{N, Int}(ntuple(i -> max(1, floor(Int, lcell * face_distances[i] / cutoff)), Val(N)))
    n_search = SVector{N, Int}(ntuple(i -> max(1, ceil(Int, cutoff * nc[i] / face_distances[i])), Val(N)))
    inv_unit_cell = inv(m ./ oneunit(T))
    # The computing box contains all particles (wrapped to the unit cell), plus the cutoff
    xmin, xmax = cell_limits(m)
    computing_box = (xmin .+ origin .- cutoff, xmax .+ origin .+ cutoff)
    cutoff_sqr = cutoff^2
    return Box{UnitCellType, N, T, typeof(cutoff_sqr), N * N, eltype(inv_unit_cell)}(
        input_unit_cell = input_unit_cell,
        inv_unit_cell = inv_unit_cell,
        lcell = lcell,
        nc = nc,
        n_search = n_search,
        cutoff = cutoff,
        cutoff_sqr = cutoff_sqr,
        computing_box = computing_box,
        origin = origin,
    )
end

# Size of the computing cells: distances between their opposite faces
cell_sizes(box::Box) = _face_distances(box.input_unit_cell.matrix) ./ box.nc

function Base.show(io::IO, ::MIME"text/plain", box::Box{UnitCellType, N}) where {UnitCellType, N}
    _println(io, "Box{$UnitCellType, $N}")
    _print(io, "  unit cell matrix = [ ")
    print(io, join(_uround.(box.input_unit_cell.matrix[1:N, 1]), " "))
    for i in 2:N
        print(io, "; ", join(_uround.(box.input_unit_cell.matrix[1:N, i]), " "))
    end
    println(io, " ]")
    _println(io, "  cutoff = ", box.cutoff)
    _println(io, "  number of computing cells on each dimension = ", box.nc)
    _println(
        io, "  computing cell sizes = [",
        join(_uround.(cell_sizes(box)), ", "), "] (lcell: ", box.lcell, ")"
    )
    return _print(io, "  Total number of cells = ", prod(box.nc))
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

#=
    Box(sides::AbstractVector, cutoff, lcell::Int=1, UnitCellType=OrthorhombicCell)

For orthorhombic unit cells, `Box` can be initialized with a vector of the length of each side.

=#
function Box(sides::AbstractVector, cutoff, lcell::Int, ::Type{UnitCellType}; origin = nothing) where {UnitCellType}
    T = _promote_types(sides, cutoff)
    unit_cell_matrix = cell_matrix_from_sides(T.(sides))
    return Box(unit_cell_matrix, T(cutoff), lcell, UnitCellType; origin)
end
Box(sides::AbstractVector, cutoff; lcell::Int = 1, UnitCellType = OrthorhombicCell) = Box(sides, cutoff, lcell, UnitCellType)

#=
    Box(unitcell::Limits, cutoff; lcell::Int=1)

This constructor receives the output of `limits(x)` or `limits(x,y)` where `x` and `y` are
the coordinates of the particles involved, and constructs a `Box` with size larger than
the maximum coordinates ranges of all particles plus twice the cutoff. This is used to
compute pairwise interactions in non-periodic systems. The output box is an `NonPeriodicCell`
box type, for which neighboring cells are never wrapped.

=#
function Box(unitcell::Limits, cutoff::T; lcell::Int = 1) where {T}
    sides = _sides_from_limits(unitcell, cutoff)
    return Box(sides, cutoff, lcell, NonPeriodicCell; origin = unitcell.origin)
end

# Types of input variables that are acceptable as a unitcell, to construct
# the unitcell matrix.
const InputUnitCellTypes = Union{Nothing, AbstractVector, AbstractMatrix, Limits, Tuple}

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
        box::Box{UnitCellType, N, T, TSQ, M};
        unitcell::InputUnitCellTypes = nothing,
        cutoff::Union{Nothing, T} = nothing,
        lcell::Union{Nothing, Int} = nothing
    ) where {UnitCellType, N, T, TSQ, M}
    _lcell = isnothing(lcell) ? box.lcell : lcell
    _cutoff = isnothing(cutoff) ? box.cutoff : cutoff
    _unitcell = if isnothing(unitcell)
        box.input_unit_cell
    elseif unitcell isa AbstractVector
        UnitCell{UnitCellType, N, T, M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa Tuple
        UnitCell{UnitCellType, N, T, M}(cell_matrix_from_sides(unitcell))
    elseif unitcell isa AbstractMatrix
        UnitCell{UnitCellType, N, T, M}(SMatrix{N, N, T, M}(unitcell))
    elseif unitcell isa Limits
        sides = _sides_from_limits(unitcell, _cutoff)
        UnitCell{UnitCellType, N, T, M}(cell_matrix_from_sides(sides))
    end
    _origin = if unitcell isa Limits
        SVector{N, T}(unitcell.origin)
    elseif isnothing(unitcell)
        box.origin
    else
        zeros(SVector{N, T})
    end
    return _construct_box(_unitcell, _lcell, _cutoff; origin = _origin)
end

#=
    _face_distances(unit_cell_matrix::SMatrix)

# Extended help

Returns the distances between opposite faces of the unit cell, associated with each
lattice vector (the height of the cell in the direction normal to the face spanned by
the other lattice vectors).

=#
function _face_distances(unit_cell_matrix::SMatrix{3, 3, T, 9}) where {T}
    a = unit_cell_matrix[:, 1]
    b = unit_cell_matrix[:, 2]
    c = unit_cell_matrix[:, 3]
    bc = cross(b, c)
    ca = cross(c, a)
    ab = cross(a, b)
    return SVector(abs(dot(a, bc)) / norm(bc), abs(dot(b, ca)) / norm(ca), abs(dot(c, ab)) / norm(ab))
end

function _face_distances(unit_cell_matrix::SMatrix{2, 2, T, 4}) where {T}
    a = unit_cell_matrix[:, 1]
    b = unit_cell_matrix[:, 2]
    i = a / norm(a)
    j = b / norm(b)
    return SVector(sqrt(max(zero(T)^2, sum(abs2, a) - dot(a, j)^2)), sqrt(max(zero(T)^2, sum(abs2, b) - dot(b, i)^2)))
end

#=
    check_unit_cell(unit_cell_matrix, cutoff)

# Extended help

Checks if the unit cell is large enough relative to the cutoff for the minimum-image
convention to apply, that is, if at most one image of each particle can be within the
cutoff of any other particle. This is not required (smaller cells are supported), but
in that case pairs are reported once for each image within the cutoff.

=#
check_unit_cell(box::Box) = check_unit_cell(box.input_unit_cell.matrix, box.cutoff)
check_unit_cell(unit_cell_matrix::SMatrix, cutoff) = all(>(2 * cutoff), _face_distances(unit_cell_matrix))

#=
    _cell_offsets!(half, full, box)

# Extended help

Computes the offsets (in units of computing cells) of the neighboring cells that may contain
particles within the cutoff of the particles of a cell. `full` contains all such offsets, including
the zero offset, and `half` contains one of each pair `{v, -v}` of nonzero offsets (the one with
the first nonzero component positive), to visit each pair of cells only once in self-computations.
Offsets for which the minimum distance between the cells is greater than the cutoff are discarded.
The vectors are emptied and filled in place (not allocating if they have enough capacity).

=#
function _cell_offsets!(half::Vector{SVector{N, Int}}, full::Vector{SVector{N, Int}}, box::Box{UnitCellType, N, T}) where {UnitCellType, N, T}
    empty!(half)
    empty!(full)
    m = box.input_unit_cell.matrix ./ oneunit(T)
    # vectors spanning a computing cell
    C = m .* (1 ./ box.nc')
    cutoff = box.cutoff / oneunit(T)
    for I in CartesianIndices(ntuple(i -> -box.n_search[i]:box.n_search[i], Val(N)))
        o = SVector{N, Int}(Tuple(I))
        _min_cell_distance(C, o) > cutoff * (1 + sqrt(eps(typeof(cutoff)))) && continue
        push!(full, o)
        _lexicographically_positive(o) && push!(half, o)
    end
    return half, full
end

@inline function _lexicographically_positive(s::SVector{N}) where {N}
    for k in 1:N
        s[k] > 0 && return true
        s[k] < 0 && return false
    end
    return false
end

#=
    _min_cell_distance(C, o)

# Extended help

Minimum distance between two parallelepiped cells spanned by the columns of `C`, separated
by the integer offset `o`, that is, the minimum of `|C t|` for `t` in `o .+ [-1, 1]^N`.
This is a box-constrained convex quadratic problem, solved by coordinate descent.

=#
function _min_cell_distance(C::SMatrix{N, N, T}, o::SVector{N, Int}) where {N, T}
    G = C' * C
    lo = SVector{N, T}(o) .- 1
    hi = SVector{N, T}(o) .+ 1
    t = clamp.(zero(SVector{N, T}), lo, hi)
    for _ in 1:100
        for d in 1:N
            s = zero(T)
            for e in 1:N
                e == d && continue
                s += G[d, e] * t[e]
            end
            t = setindex(t, clamp(-s / G[d, d], lo[d], hi[d]), d)
        end
    end
    return sqrt(max(zero(T), dot(t, G * t)))
end
