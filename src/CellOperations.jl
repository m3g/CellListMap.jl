
#= 

Function that validates the coordinates. 

=#
function _validate_coordinates(x)
    for (i, v) in enumerate(x)
        if any(isnan, v) || any(ismissing, v)
            throw(ArgumentError("""\n

                Invalid coordinates found: $v for particle of index $i.

            """))
        end
    end
    return nothing
end



#=
    fastmod1(x)

Computes `mod(x,1)`, quickly, using `x - floor(x)`. Maybe irrelevant.

=#
@inline fastmod1(x) = x - floor(x)

#=
    wrap_cell_fraction(x,unit_cell_matrix)

# Extended help

Obtaint the coordinates of `x` as a fraction of unit cell vectors, first
positive cell. `x` is a vector of dimension `N` and `cell` a matrix of 
dimension `NxN`

## Example

```julia-repl
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_cell_fraction(x,unit_cell_matrix)
2-element Vector{Float64}:
 0.5
 0.3
```

=#
@inline function wrap_cell_fraction(x::AbstractVector{T}, unit_cell_matrix::AbstractMatrix) where {T}
    # Division by `oneunit` is to support Unitful quantities. 
    # this workaround works here because the units cancel.
    # see: https://github.com/PainterQubits/Unitful.jl/issues/46
    x_stripped = x ./ oneunit(T)
    m_stripped = unit_cell_matrix ./ oneunit(T)
    p = fastmod1.(m_stripped \ x_stripped)
    # Boundary coordinates belong to the lower boundary
    p = ifelse.(p .== one(eltype(x_stripped)), zero(eltype(x_stripped)), p)
    return p
end

#=
    wrap_to_first(x,unit_cell_matrix)

# Extended help

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates. 

## Example

```julia-repl
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_to_first(x,unit_cell_matrix)
2-element Vector{Float64}:
 5.0
 3.0000000000000004
```

=#
@inline function wrap_to_first(x, unit_cell_matrix)
    p = wrap_cell_fraction(x, unit_cell_matrix)
    return unit_cell_matrix * p
end

@testitem "wrap_to_first" begin
    x = [15.0, 13.0]
    unit_cell_matrix = [10.0 0.0; 0.0 10.0]
    @test CellListMap.wrap_to_first(x, unit_cell_matrix) ≈ [5.0, 3.0]
end

#=
    wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}

# Extended help

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

=#
@inline function wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    invu = inv(oneunit(T))
    unit_cell_matrix = invu * unit_cell_matrix
    x_f = wrap_cell_fraction(invu * x, unit_cell_matrix)
    xref_f = wrap_cell_fraction(invu * xref, unit_cell_matrix)
    xw = wrap_relative_to(x_f, xref_f, SVector{N,eltype(x_f)}(ntuple(i -> 1, Val(N))))
    return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref
end

#=
    wrap_relative_to(x,xref,sides::AbstractVector)

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`,
for an Orthorhombic cell of which only the sides are provided.

=#
@inline function wrap_relative_to(x, xref, sides::AbstractVector)
    lengths = abs.(sides)
    xw = mod.(x - xref, lengths)
    # Wrap a single coordinate: there is no need for the x <= -s/2 case because
    # the mod function already takes care of it, with positive lengths.
    xw = @. ifelse(xw >= lengths / 2, xw - lengths, xw)
    return xw + xref
end

@testitem "wrap_relative_to" begin
    using StaticArrays
    using Unitful
    for (x, y, xy, yx) in [
        ([15.0, 13.0], [4.0, 2.0], [5.0, 3.0], [14.0, 12.0]),
        ([-7.0, -6.0], [1.0, 2.0], [3.0, 4.0], [-9.0, -8.0]),
    ]
        # triclinic cells: unitcell is a matrix
        unit_cell_matrix = SMatrix{2,2}(10.0, 0.0, 0.0, 10.0)
        @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, unit_cell_matrix) ≈ yx
        unit_cell_matrix = SMatrix{2,2}(-10.0, 0.0, 0.0, -10.0)
        @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, unit_cell_matrix) ≈ yx
        # orthorhombic cells: sides is a vector
        sides = [10.0, 10.0]
        @test CellListMap.wrap_relative_to(x, y, sides) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, sides) ≈ yx
        sides = [-10.0, -10.0]
        @test CellListMap.wrap_relative_to(x, y, sides) ≈ xy
        @test CellListMap.wrap_relative_to(y, x, sides) ≈ yx
    end
    # test unit propagations
    x = [15.0, 13.0]u"nm"
    y = [4.0, 2.0]u"nm"
    unit_cell_matrix = SMatrix{2,2}(10.0, 0.0, 0.0, 10.0)u"nm"
    @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ [5.0, 3.0]u"nm"
    unit_cell_matrix = [10.0, 10.0]u"nm"
    @test CellListMap.wrap_relative_to(x, y, unit_cell_matrix) ≈ [5.0, 3.0]u"nm"
end

#=
    translation_image(x::SVector{N,T},unit_cell_matrix,indices) where {N,T}

# Extended help

Translate vector `x` according to the `unit_cell_matrix` lattice vectors and the `indices`
provided.

=#
@inline translation_image(x::SVector{N,T}, unit_cell_matrix, indices) where {N,T} =
    x + unit_cell_matrix * SVector{N,Int}(ntuple(i -> indices[i], Val(N)))

#=
    translation_image(x::AbstractVector{<:AbstractVector},unit_cell_matrix,indices)

Translates a complete set of coordinates given a set of indexes of unit-cells. Returns a new
set of coordinates. 

## Example

```julia-repl
julia> x = rand(SVector{2,Float64},100);

julia> box = Box([1,1],0.1);

julia> CellListMap.translation_image(x,box.unit_cell.matrix,(1,1))
100-element Vector{SVector{2, Float64}}:
 [1.847791110439223, 1.5989103939725295]
 [1.3493293666090889, 1.4002971843576644]
 [1.4111736701313218, 1.3471780214994182]
 ⋮
 [1.1548437388991908, 1.7034501001177493]
 [1.4066300885242247, 1.2907398318754952]
```

=#
function translation_image(x::AbstractVector{<:AbstractVector}, unit_cell_matrix, indices)
    x_new = similar(x)
    for i in eachindex(x)
        x_new[i] = translation_image(x[i], unit_cell_matrix, indices)
    end
    return x_new
end

@testitem "translation image" begin
    using StaticArrays
    x = SVector{2}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    unitcell = [10.0 0.0; 0.0 10.0]
    @test CellListMap.translation_image(x, unitcell, [0, 0]) ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    @test CellListMap.translation_image(x, unitcell, [1, 1]) ≈ [[11.0, 11.0], [12.0, 12.0], [13.0, 13.0]]
    @test CellListMap.translation_image(x, unitcell, [1, 2]) ≈ [[11.0, 21.0], [12.0, 22.0], [13.0, 23.0]]
    @test CellListMap.translation_image(x, unitcell, [2, 1]) ≈ [[21.0, 11.0], [22.0, 12.0], [23.0, 13.0]]
    @test CellListMap.translation_image(x, unitcell, [-1, -1]) ≈ [[-9.0, -9.0], [-8.0, -8.0], [-7.0, -7.0]]
end

#=
    replicate_system!(
        x::AbstractVector{SVector{N,T}},
        unit_cell_matrix::AbstractMatrix,
        ranges::Tuple
    ) where {N,T}

# Extended help

Replicate the system (modifying the original array of coordinates) in all
directions defined by the periodic system and by the range of unitary cells 
of interest. `x` can be a `(N,M)` matrix, and the unit cell matrix can be
provided instead of the `box`.

## Example

```julia-repl
julia> x = rand(SVector{2,Float64},100);

julia> box = Box([1,1],0.1);

julia> CellListMap.replicate_system!(x,box,(0:0,-1:1))
300-element Vector{SVector{2, Float64}}:
 [0.7119987163255118, 0.6788616154460262]
 [0.6188407316804118, 0.8497116428720384]
 [0.21328895963244354, 0.48932085643862977]
 ⋮
 [0.4114499470191678, 1.1034376619603892]
 [0.6094126258851252, 1.2328989485215263]
```

=#
function replicate_system!(
    x::AbstractVector{SVector{N,T}},
    unit_cell_matrix::AbstractMatrix,
    ranges::Tuple
) where {N,T}
    length(ranges) == N || throw(DimensionMismatch("Tuple of ranges must have the same dimension as the vectors: $N"))
    i0 = ntuple(i -> 0, Val(N))
    imgs = Iterators.filter(!isequal(i0),
        Iterators.product(ranges...)
    )
    x0 = copy(x)
    for img in imgs
        x_new = translation_image(x0, unit_cell_matrix, img)
        append!(x, x_new)
    end
    return x
end

# Warning: this function is non-mutating, because the matrices cannot be
# resized in-place.
function replicate_system(x::AbstractMatrix{T}, cell, ranges) where {T}
    N = size(x, 1)
    x_re = [SVector{N,T}(ntuple(i -> x[i, j], Val(N))) for j in axes(x, 2)]
    replicate_system!(x_re, cell, ranges)
    x = Matrix(reinterpret(reshape, Float64, x_re))
    return x
end

@testitem "replicate system" begin
    using StaticArrays
    using CellListMap: replicate_system!, replicate_system
    unitcell = [10.0 0.0; 0.0 10.0]
    x = SVector{2,Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    replicate_system!(x, unitcell, (1, 0))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [11.0, 1.0], [12.0, 2.0], [13.0, 3.0]]
    x = SVector{2,Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    replicate_system!(x, unitcell, (1, 1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [11.0, 11.0], [12.0, 12.0], [13.0, 13.0]]
    x = SVector{2,Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    replicate_system!(x, unitcell, (0, 1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [1.0, 11.0], [2.0, 12.0], [3.0, 13.0]]
    x = SVector{2,Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    replicate_system!(x, unitcell, (-1, -1))
    @test x ≈ [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [-9.0, -9.0], [-8.0, -8.0], [-7.0, -7.0]]
    # with a matrix input
    x = [1.0 2.0 3.0; 1.0 2.0 3.0]
    y = replicate_system(x, unitcell, (1, 0))
    @test y ≈ [1.0 2.0 3.0 11.0 12.0 13.0; 1.0 2.0 3.0 1.0 2.0 3.0]
    # throw error if dimensions do not match
    x = SVector{2,Float64}[[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
    @test_throws DimensionMismatch replicate_system!(x, unitcell, (1, 0, 1))
end


#=
    cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}

# Extended help

Given the linear index of the cell in the cell list, returns the cartesian indices 
of the cell (for arbitrary dimension N).

=#
@inline cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} =
    CartesianIndices(ntuple(i -> nc[i], Val(N)))[i1D]

#=
    cell_linear_index(nc::SVector{N,Int}, indices) where N

# Extended help

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

=#
@inline cell_linear_index(nc::SVector{N,Int}, indices) where {N} =
    LinearIndices(ntuple(i -> nc[i], Val(N)))[ntuple(i -> indices[i], Val(N))...]

@testitem "cell cartesian/linear indices" begin
    using StaticArrays
    using CellListMap: cell_cartesian_indices, cell_linear_index
    @test cell_cartesian_indices(SVector(3, 3), 1) == CartesianIndex(1, 1)
    @test cell_cartesian_indices(SVector(3, 3), 2) == CartesianIndex(2, 1)
    @test cell_cartesian_indices(SVector(3, 3), 4) == CartesianIndex(1, 2)
    @test cell_cartesian_indices(SVector(3, 3), 9) == CartesianIndex(3, 3)
    @test cell_linear_index(SVector(3, 3), CartesianIndex(1, 1)) == 1
    @test cell_linear_index(SVector(3, 3), CartesianIndex(2, 1)) == 2
    @test cell_linear_index(SVector(3, 3), CartesianIndex(1, 2)) == 4
    @test cell_linear_index(SVector(3, 3), CartesianIndex(3, 3)) == 9
end

#
# Compute the maximum and minimum coordinates of the vectors composing the particle sets
#
function _minmax(x::AbstractVector{<:AbstractVector})
    length(x) <= 0 && throw(ArgumentError("Cannot set unitcell box from empty coordinates vector."))
    N = size(x[begin], 1)
    T = eltype(x[begin])
    xmin = fill(typemax(T), MVector{N,T})
    xmax = fill(typemin(T), MVector{N,T})
    for v in x
        @. xmin = min(xmin, v)
        @. xmax = max(xmax, v)
    end
    return SVector(xmin), SVector(xmax)
end

@testitem "_minmax" setup=[AllocTest] begin
    using BenchmarkTools
    using StaticArrays
    import CellListMap: _minmax
    using .AllocTest: Allocs
    x = [[0.0, 0.5, 1.0], [0.5, 1.0, 0.0], [1.0, 0.0, 0.5]]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    x = [SVector(0.0, 0.5, 1.0), SVector(0.5, 1.0, 0.0), SVector(1.0, 0.0, 0.5)]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    a = @ballocated _minmax($x) evals = 1 samples = 1
    @test a == Allocs(0)
end

"""
    limits(x; validate_coordinates::Union{Nothing,Function})

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`, 
to be used to set a box without effective periodic boundary conditions.

The `validate_coordinates` function is used to validate the coordinates of the particles.
By default, it will throw an error if any of the coordinates contain `NaN` or `missing` values.
To disable this validation, set `validate_coordinates = nothing`. Custom checks can be implemented
by passing a function that takes the coordinates as input and throws an error if the coordinates
are invalid.

"""
function limits(x::AbstractVector{<:AbstractVector}; validate_coordinates::Union{Nothing,Function}=_validate_coordinates)
    isnothing(validate_coordinates) || validate_coordinates(x)
    xmin, xmax = _minmax(x)
    return Limits(xmax .- xmin)
end

function limits(x::AbstractMatrix; validate_coordinates::Union{Nothing,Function}=_validate_coordinates)
    N = size(x, 1)
    (N == 2 || N == 3) || throw(DimensionMismatch("The first dimension of the matrix must be the dimension (2 or 3)"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return limits(x_re; validate_coordinates)
end

"""
    limits(x,y; validate_coordinates::Union{Nothing, Function})

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`
and `y`, to used to set a box without effective periodic boundary conditions.

The `validate_coordinates` function is used to validate the coordinates of the particles.
By default, it will throw an error if any of the coordinates contain `NaN` or `missing` values.
To disable this validation, set `validate_coordinates = nothing`. Custom checks can be implemented
by passing a function that takes the coordinates as input and throws an error if the coordinates
are invalid.

"""
function limits(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    validate_coordinates::Union{Nothing,Function}=_validate_coordinates
)
    isnothing(validate_coordinates) || validate_coordinates(x)
    isnothing(validate_coordinates) || validate_coordinates(y)
    xmin, xmax = _minmax(x)
    ymin, ymax = _minmax(y)
    xymin = min.(xmin, ymin)
    xymax = max.(xmax, ymax)
    return Limits(xymax .- xymin)
end

function limits(x::AbstractMatrix, y::AbstractMatrix; validate_coordinates::Union{Nothing,Function}=_validate_coordinates)
    N = size(x, 1)
    M = size(y, 1)
    N == M || throw(DimensionMismatch("The first dimension of the input matrices must be equal. "))
    (N == 2 || N == 3) || throw(DimensionMismatch("The first dimension of the matrix must be the dimension (2 or 3)"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return limits(x_re, y_re; validate_coordinates)
end

@testitem "limits - invalid coordinates" begin
    using CellListMap: limits
    using StaticArrays: SVector
    x = rand(SVector{3,Float64}, 100)
    x[50] = SVector(1.0, NaN, 1.0)
    y = rand(SVector{3,Float64}, 100)
    @test_throws ArgumentError limits(x)
    @test_throws ArgumentError limits(x, y)
    @test_throws ArgumentError limits(y, x)
    x = rand(3, 100)
    x[2, 50] = NaN
    y = rand(3, 100)
    @test_throws ArgumentError limits(x)
    @test_throws ArgumentError limits(x, y)
    @test_throws ArgumentError limits(y, x)
end


#=
    align_cell(m::StaticMatrix)
    align_cell!(m::AbstractMatrix)

# Extended help

These functions rotate the unit cell matrix such that the largest lattice vector is oriented
along the x-axis and, for 3D cells, also that the the plane formed by the largest and 
second largest lattice vectors is oriented perpendicular to the z-axis. 

=#
function align_cell end
align_cell(m::AbstractMatrix) = align_cell!(copy(m))
align_cell(m::SMatrix{2}) = _align_cell2D!(m)
align_cell(m::SMatrix{3}) = _align_cell3D!(m)

function align_cell!(m::AbstractMatrix)
    if size(m) == (2, 2)
        m, R = _align_cell2D!(m)
    elseif size(m) == (3, 3)
        m, R = _align_cell3D!(m)
    else
        throw(ArgumentError("align_cell! only supports square matrices in 2 or 3 dimensions."))
    end
    return m, R
end

function _align_cell2D!(m::AbstractMatrix{T}) where {T}
    m = m ./ oneunit(T)
    x, y = 1, 2
    a = @view(m[:, 1])
    b = @view(m[:, 2])
    if norm(b) > norm(a)
        a = b
    end
    # cell is already properly rotated
    if a[y] ≈ zero(T)
        R = one(m)
    else
        # rotate first axis to be parallel to x (clockwise)
        sinθ = -norm(a) / (a[x]^2 / a[y] + a[y])
        cosθ = -a[x] * sinθ / a[y]
        #! format: off
        R = @SMatrix[cosθ -sinθ 
                     sinθ  cosθ]
        #! format: on
        m = R * m
    end
    return oneunit(T) .* m, R
end

function _align_cell3D!(m::AbstractMatrix{T}) where {T}
    m = inv(oneunit(T)) * m
    x, y, z = 1, 2, 3
    # Choose a and b to be the largest lattice vectors, in order
    n = SVector{3}(norm_sqr(v) for v in eachcol(m))
    combinations = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
    local a, b, c, axis_on_x
    for (i, j, k) in combinations
        if n[i] >= n[j] >= n[k]
            a, b, c = (@view(m[:, i]), @view(m[:, j]), @view(m[:, k]))
            axis_on_x = i
            break
        end
    end

    # Find rotation that aligns a with x
    a1 = normalize(a)
    v = SVector(0, a1[z], -a1[y]) # a1 × i 
    if norm_sqr(v) ≈ 0
        R1 = one(m)
    else
        #! format: off
        vₛ = @SMatrix[     0  -v[z]    v[y]
                        v[z]      0   -v[x]
                       -v[y]   v[x]       0  ]
        #! format: on
        R1 = one(m) + vₛ + vₛ^2 * inv(1 + a1[x])
    end
    m = R1 * m
    # Find rotation along x-axis that makes b orthogonal to z
    x, y, z = @view(m[:, 2])
    if (y^2 + z^2) ≈ 0
        R2 = one(m)
    else
        a = x
        b = sqrt(norm_sqr(m[:, 2]) - x^2)
        sinθ = -z * b / (y^2 + z^2)
        cosθ = sqrt(1 - sinθ^2)
        #! format: off
        R2 = @SMatrix[ 1    0     0
                       0 cosθ -sinθ
                       0 sinθ  cosθ ]
        #! format: on
    end
    m = R2 * m
    return oneunit(T) .* m, R2 * R1
end

@testitem "align_cell" begin
    using CellListMap.TestingNeighborLists: random_rotation
    import CellListMap: align_cell
    using StaticArrays
    using LinearAlgebra

    l = sqrt(2) / 2

    m = @SMatrix[1.0 0.0; 0.0 1.0]
    @test align_cell(m) == (one(m), one(m))

    m = @SMatrix[l 0; l 1]
    mt, R = align_cell(m)
    @test mt ≈ [1 l; 0 l]
    @test R ≈ [l l; -l l]

    m = @SMatrix[-l 0; l 1]
    mt, R = align_cell(m)
    @test mt ≈ [1 l; 0 -l]
    @test R ≈ [-l l; -l -l]

    #! format: off
    m = @SMatrix[ 
        1  0  0
        0  1  0 
        0  0  1
    ]
    #! format: on
    @test align_cell(m) == (one(m), one(m))

    #! format: off
    m = @SMatrix[ 
        3  0  0
        0  2  0 
        0  0  1
    ]
    #! format: on
    for _ in 1:5
        local R
        R = random_rotation()
        mr = R * m
        ma, Ra = align_cell(mr)
        @test ma[:, 1] ≈ m[:, 1]
        @test cross([1, 0, 0], cross(ma[:, 2], ma[:, 3])) ≈ zeros(3) atol = 1e-10
    end

    # throw error if not 2D or 3D
    @test_throws ArgumentError align_cell(rand(4, 4))

end

#=
    cell_limits(m::AbstractMatrix)

For 2D and 3D matrices, returns the maximum and minimum coordinates of all vertices. 

=#
function cell_limits(m::AbstractMatrix{T}) where {T}
    vertices = cell_vertices(m)
    xmin = MVector(vertices[begin])
    xmax = MVector(vertices[begin])
    for v in @view(vertices[begin+1:end])
        for j in eachindex(v)
            xmin[j] = min(xmin[j], v[j])
            xmax[j] = max(xmax[j], v[j])
        end
    end
    return SVector(xmin), SVector(xmax)
end

@testitem "cell_limts" begin
    import CellListMap: cell_limits, align_cell
    m = [1 0; 0 1]
    @test cell_limits(m) == ([0.0, 0.0], [1.0, 1.0])

    m = [10 5; 5 10]
    @test cell_limits(m) == ([0.0, 0.0], [15.0, 15.0])

    mr, _ = align_cell(m)
    cl = cell_limits(mr)
    @test cl[1] ≈ [0.0, 0.0] atol = 1e-10
    @test cl[2] ≈ [20.12461179749811, 6.708203932499369]

    m = [10 5; 0 10]
    @test cell_limits(m) == ([0.0, 0.0], [15.0, 10.0])

    mr, _ = align_cell(m)
    cl = cell_limits(mr)
    @test cl[1] ≈ [0.0, -8.94427190999916]
    @test cl[2] ≈ [15.652475842498529, 0.0]

    m = [1 0 0; 0 1 0; 0 0 1]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])

    m = [1 0 0; 0 2 0; 0 0 1]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 2.0, 1.0])

    mr, _ = align_cell(m)
    @test cell_limits(mr) == ([0.0, -1.0, 0.0], [2.0, 0.0, 1.0])

    m = [1 0 0; 0 2 0; 0 0 3]
    @test cell_limits(m) == ([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])

    mr, _ = align_cell(m)
    @test cell_limits(mr) == ([0.0, 0.0, -1.0], [3.0, 2.0, 0.0])
end

#=
    cell_vertices(m::AbstractMatrix)

Function that returns the vertices of a unit cell in 2D or 3D, given the unit cell matrix.

=#
function cell_vertices(m::AbstractMatrix)
    if size(m) == (2, 2)
        x = _cell_vertices2D(m)
    elseif size(m) == (3, 3)
        x = _cell_vertices3D(m)
    else
        throw(ArgumentError("cell_vertices only supports square matrices in 2 or 3 dimensions."))
    end
end

@views function _cell_vertices2D(m::AbstractMatrix{T}) where {T}
    S = SVector{2,T}
    x = SVector{4,S}(S(zero(T), zero(T)), S(m[:, 1]), S(m[:, 1]) + S(m[:, 2]), S(m[:, 2]))
    return x
end

@views function _cell_vertices3D(m::AbstractMatrix{T}) where {T}
    S = SVector{3,T}
    x = SVector{8,S}(
        S(zero(T), zero(T), zero(T)),
        S(m[:, 1]),
        S(m[:, 1]) + S(m[:, 2]),
        S(m[:, 2]),
        S(m[:, 1]) + S(m[:, 3]),
        S(m[:, 3]),
        S(m[:, 2]) + S(m[:, 3]),
        S(m[:, 1]) + S(m[:, 2]) + S(m[:, 3]),
    )
    return x
end

@testitem "cell vertices" begin
    import CellListMap: cell_vertices
    m = [1 0; 0 1]
    @test all(
        isapprox.(cell_vertices(m), [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], atol=1e-6)
    )
    m = [1 0.5; 0 1]
    @test all(
        isapprox.(cell_vertices(m), [[0.0, 0.0], [1.0, 0.0], [1.5, 1.0], [0.5, 1.0]], atol=1e-6)
    )
    m = [1 0 0; 0 1 0; 0 0 1]
    @test all(
        isapprox.(cell_vertices(m), 
        [ [0, 0, 0],
          [1, 0, 0],
          [1, 1, 0],
          [0, 1, 0],
          [1, 0, 1],
          [0, 0, 1],
          [0, 1, 1],
          [1, 1, 1] ],
        atol = 1e-6)
    )
    m = [1 0.5 0; 0 1 0; 0 0 1]
    @test all(
        isapprox.(cell_vertices(m), 
        [ [0.0, 0.0, 0.0],
          [1.0, 0.0, 0.0],
          [1.5, 1.0, 0.0],
          [0.5, 1.0, 0.0],
          [1.0, 0.0, 1.0],
          [0.0, 0.0, 1.0],
          [0.5, 1.0, 1.0],
          [1.5, 1.0, 1.0] ],
        atol=1e-6)
    )
    # throw error if not 2D or 3D
    @test_throws ArgumentError cell_vertices(rand(4, 4))
end

# Obs: auxiliary functions to draw unit cells are available in the 
# testing.jl file.
