"""
    fastmod1(x)

$(INTERNAL)

Computes `mod(x,1)`, quickly, using `x - floor(x)`. Maybe irrelevant.

"""
@inline fastmod1(x) = x - floor(x)

"""
    wrap_cell_fraction(x,unit_cell_matrix)

$(INTERNAL)

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

"""
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

"""
    wrap_to_first(x,unit_cell_matrix)

$(INTERNAL)

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

"""
@inline function wrap_to_first(x, unit_cell_matrix)
    p = wrap_cell_fraction(x, unit_cell_matrix)
    return unit_cell_matrix * p
end

"""
    wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}

$(INTERNAL)

# Extended help

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. 

"""
@inline function wrap_relative_to(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    invu = inv(oneunit(T))
    unit_cell_matrix = invu * unit_cell_matrix
    x_f = wrap_cell_fraction(invu*x, unit_cell_matrix)
    xref_f = wrap_cell_fraction(invu*xref, unit_cell_matrix)
    xw = wrap_relative_to(x_f, xref_f, SVector{N,eltype(x_f)}(ntuple(i -> 1, N)))
    return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref
end

"""
    wrap_relative_to(x,xref,sides::AbstractVector)

$(INTERNAL)

# Extended help

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`,
for an Orthorhombic cell of which only the side lengths are provided.

"""
@inline function wrap_relative_to(x, xref, sides::AbstractVector)
    xw = mod.(x - xref, sides)
    xw = _wrap_single_coordinate.(xw, sides)
    return xw + xref
end

#
# Wrap a single coordinate
#
@inline function _wrap_single_coordinate(x, s)
    if x >= s / 2
        x = x - s
    elseif x < -s / 2
        x = x + s
    end
    return x
end

"""
    translation_image(x::SVector{N,T},unit_cell_matrix,indices) where {N,T}

$(INTERNAL)

# Extended help

Translate vector `x` according to the `unit_cell_matrix` lattice vectors and the `indices`
provided.

"""
@inline translation_image(x::SVector{N,T}, unit_cell_matrix, indices) where {N,T} =
    x + unit_cell_matrix * SVector{N,Int}(ntuple(i -> indices[i], N))

"""
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

"""
function translation_image(x::AbstractVector{<:AbstractVector}, unit_cell_matrix, indices)
    x_new = similar(x)
    for i in eachindex(x)
        x_new[i] = translation_image(x[i], unit_cell_matrix, indices)
    end
    return x_new
end

"""
    replicate_system!(
        x::AbstractVector{SVector{N,T}},
        unit_cell_matrix::AbstractMatrix,
        ranges::Tuple
    ) where {N,T}

$(INTERNAL)

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

"""
function replicate_system!(
    x::AbstractVector{SVector{N,T}},
    unit_cell_matrix::AbstractMatrix,
    ranges::Tuple
) where {N,T}
    length(ranges) == N || throw(DimensionMismatch("Tuple of ranges must have the same dimension as the vectors: $N"))
    i0 = ntuple(i -> 0, N)
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

function replicate_system!(x::AbstractMatrix{T}, cell, ranges) where {T}
    N = size(x, 1)
    x_re = [SVector{N,T}(ntuple(i -> x[i, j], N)) for j in axes(x, 2)]
    replicate_system!(x_re, cell, ranges)
    x = Matrix(reinterpret(reshape, Float64, x_re))
    return x
end


"""
    cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}

$(INTERNAL)

# Extended help

Given the linear index of the cell in the cell list, returns the cartesian indices 
of the cell (for arbitrary dimension N).

"""
@inline cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} =
    CartesianIndices(ntuple(i -> nc[i], N))[i1D]

"""
    cell_linear_index(nc::SVector{N,Int}, indices) where N

$(INTERNAL)

# Extended help

Returns the index of the cell, in the 1D representation, from its cartesian coordinates. 

"""
@inline cell_linear_index(nc::SVector{N,Int}, indices) where {N} =
    LinearIndices(ntuple(i -> nc[i], N))[ntuple(i -> indices[i], N)...]

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

@testitem "_minmax" begin
    using BenchmarkTools
    using StaticArrays
    import CellListMap: _minmax
    x = [[0.0, 0.5, 1.0], [0.5, 1.0, 0.0], [1.0, 0.0, 0.5]]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    x = [SVector(0.0, 0.5, 1.0), SVector(0.5, 1.0, 0.0), SVector(1.0, 0.0, 0.5)]
    @test _minmax(x) === (SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
    a = @ballocated _minmax($x) evals = 1 samples = 1
    @test a == 0
end

"""
    limits(x)

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`, 
to be used to set a box without effective periodic boundary conditions.

"""
function limits(x::AbstractVector{<:AbstractVector})
    xmin, xmax = _minmax(x)
    return Limits(xmax .- xmin)
end

function limits(x::AbstractMatrix)
    N = size(x, 1)
    (N == 2 || N == 3) || throw(DimensionMismatch("The first dimension of the matrix must be the dimension (2 or 3)"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    return limits(x_re)
end

"""
    limits(x,y)

Returns the lengths of a orthorhombic box that encompasses all the particles defined in `x`
and `y`, to used to set a box without effective periodic boundary conditions.

"""
function limits(x::T, y::T) where {T<:AbstractVector{<:AbstractVector}}
    xmin, xmax = _minmax(x)
    ymin, ymax = _minmax(y)
    xymin = min.(xmin, ymin)
    xymax = max.(xmax, ymax)
    return Limits(xymax .- xymin)
end

function limits(x::T, y::T) where {T<:AbstractMatrix}
    N = size(x, 1)
    M = size(y, 1)
    N == M || throw(DimensionMismatch("The first dimension of the input matrices must be equal. "))
    (N == 2 || N == 3) || throw(DimensionMismatch("The first dimension of the matrix must be the dimension (2 or 3)"))
    x_re = reinterpret(reshape, SVector{N,eltype(x)}, x)
    y_re = reinterpret(reshape, SVector{N,eltype(y)}, y)
    return limits(x_re, y_re)
end


"""
    align_cell(m::StaticMatrix)
    align_cell!(m::AbstractMatrix)

$(INTERNAL)

# Extended help

These functions rotate the unit cell matrix such that the largest lattice vector is oriented
along the x-axis and, for 3D cells, also that the the plane formed by the largest and 
second largest lattice vectors is oriented perpendicular to the z-axis. 

"""
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
        R = random_rotation()
        mr = R * m
        ma, Ra = align_cell(mr)
        @test ma[:,1] ≈ m[:,1]
        @test cross([1,0,0], cross(ma[:,2],ma[:,3])) ≈ zeros(3) atol=1e-10
    end

end

"""
    cell_vertices(m::AbstractMatrix)

$(INTERNAL)

Function that returns the vertices of a unit cell in 2D or 3D, given the unit cell matrix.

"""
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

"""
    draw_cell_vertices(m::AbstractMatrix)

$(INTERNAL)

Function that returns the vertices of a unit cell matrix in 2D or 3D, as a vector
of static vectors, in a proper order for ploting the cell (the first vertex, in the
origin, is repeated at the end of the list, to close the figure)

"""
function draw_cell_vertices(m::AbstractMatrix)
    if size(m) == (2, 2)
        x = _draw_cell_vertices2D(m)
    elseif size(m) == (3, 3)
        x = _draw_cell_vertices3D(m)
    else
        throw(ArgumentError("draw_cell_vertices only supports square matrices in 2 or 3 dimensions."))
    end
end

function _draw_cell_vertices2D(m::AbstractMatrix)
    S = SVector{2,Float64}
    x = [
        S(0.0, 0.0),
        S(m[:, 1]),
        S(m[:, 1] + m[:, 2]),
        S(m[:, 2]),
        S(0.0, 0.0)
    ]
    return x
end

function _draw_cell_vertices3D(m::AbstractMatrix)
    S = SVector{3,Float64}
    x = [
        S(0.0, 0.0, 0.0),
        S(m[:, 1]),
        S(m[:, 1] + m[:, 2]),
        S(m[:, 2]),
        S(0.0, 0.0, 0.0),
        S(m[:, 1]),
        S(m[:, 1] + m[:, 3]),
        S(m[:, 3]),
        S(0.0, 0.0, 0.0),
        S(m[:, 2]),
        S(m[:, 2] + m[:, 3]),
        S(m[:, 2]),
        S(0.0, 0.0, 0.0),
        S(m[:, 1]),
        S(m[:, 1] + m[:, 2]),
        S(m[:, 1] + m[:, 2]) + m[:, 3],
        S(m[:, 1] + m[:, 3]),
        S(m[:, 3]),
        S(m[:, 3]) + m[:, 2],
        S(m[:, 1] + m[:, 2]) + m[:, 3]
    ]
    return x
end

@testitem "draw_cell_vertices" begin
    using StaticArrays
    import CellListMap: draw_cell_vertices
    m = [1 0; 0 1]
    @test draw_cell_vertices(m) == SVector{2,Float64}[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]
    m = [1 0 0; 0 1 0; 0 0 1]
    @test draw_cell_vertices(m) == SVector{3,Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
end

"""
    cell_limits(m::AbstractMatrix)

$(INTERNAL)

For 2D and 3D matrices, returns the maximum and minimum coordinates of all vertices. 

"""
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
    @test cell_limits(mr) == ([0.0, 0.0], [20.12461179749811, 6.708203932499369])

    m = [10 5; 0 10]
    @test cell_limits(m) == ([0.0, 0.0], [15.0, 10.0])

    mr, _ = align_cell(m)
    @test cell_limits(mr) == ([0.0, -8.94427190999916], [15.652475842498529, 0.0])

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

"""
    draw_cell(m::AbstractMatrix; aspect_ratio=:auto)

$(INTERNAL)

Draw the unit cell in a 2D or 3D plot. Requires `using Plots`.

"""
function draw_cell(m::AbstractMatrix; aspect_ratio=:auto)
    plot = Main.plot
    plot! = Main.plot!
    vertices = draw_cell_vertices(m)
    plt = plot(Tuple.(vertices), label=nothing)
    lims = cell_limits(m)
    dx = (lims[2][1] - lims[1][1]) / 10
    dy = (lims[2][2] - lims[1][2]) / 10
    plot!(plt,
        xlims=(lims[1][1] - dx, lims[2][1] + dx),
        ylims=(lims[1][2] - dy, lims[2][2] + dy),
        aspect_ratio=aspect_ratio,
        xlabel="x",
        ylabel="y",
        zlabel="z",
    )
    return plt
end



