#
# Prototype: cell lists on a grid in fractional coordinates, without replication
# of periodic images. Periodic images are reached in the pair loop: the neighbor
# cell index is wrapped with `divmod`, and the quotient is the lattice translation
# of the neighbor cell. Self-pairs, fully periodic systems only.
#
module FracCellList

using StaticArrays
using LinearAlgebra: cross, dot, norm, det
using Base.Threads: @threads, @spawn, nthreads, Atomic, atomic_add!

struct Particle{N, T}
    index::Int
    x::SVector{N, T}
end

# Distances between opposite faces of the unit cell (columns are lattice vectors)
function face_distances(M::SMatrix{3, 3})
    a, b, c = M[:, 1], M[:, 2], M[:, 3]
    bc, ca, ab = cross(b, c), cross(c, a), cross(a, b)
    return SVector(abs(dot(a, bc)) / norm(bc), abs(dot(b, ca)) / norm(ca), abs(dot(c, ab)) / norm(ab))
end
function face_distances(M::SMatrix{2, 2, T}) where {T}
    a, b = M[:, 1], M[:, 2]
    return SVector(
        sqrt(max(zero(T), sum(abs2, a) - dot(a, b / norm(b))^2)),
        sqrt(max(zero(T), sum(abs2, b) - dot(b, a / norm(a))^2)),
    )
end

# Offsets with the first nonzero component positive: exactly one of {v, -v}
function lexicographically_positive(o)
    for k in eachindex(o)
        o[k] > 0 && return true
        o[k] < 0 && return false
    end
    return false
end

# Minimum distance between two grid cells (parallelepipeds spanned by the columns
# of `C`) separated by the integer offset `o`: min |C t| for t in o .+ [-1, 1]^N.
# Box-constrained convex quadratic, solved by coordinate descent.
function min_cell_distance(C::SMatrix{N, N, T}, o) where {N, T}
    G = C' * C
    lo = SVector{N, T}(o) .- 1
    hi = SVector{N, T}(o) .+ 1
    t = clamp.(zero(SVector{N, T}), lo, hi)
    for _ in 1:200
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

struct Grid{N, T, NN}
    M::SMatrix{N, N, T, NN}
    invM::SMatrix{N, N, T, NN}
    cutoff::T
    cutoff_sqr::T
    n_cells::SVector{N, Int}
    n_search::SVector{N, Int}
    offsets::Vector{SVector{N, Int}}
end

function Grid(M::SMatrix{N, N, T}, cutoff::T; cells_per_cutoff::Int = 1) where {N, T}
    h = face_distances(M)
    n_cells = SVector{N, Int}(max.(1, floor.(Int, cells_per_cutoff .* h ./ cutoff)))
    n_search = SVector{N, Int}(max.(1, ceil.(Int, cutoff .* n_cells ./ h)))
    # cell vectors: lattice vectors divided by the number of cells in each direction
    C = M .* (one(T) ./ n_cells')
    offsets = SVector{N, Int}[]
    for I in CartesianIndices(ntuple(d -> -n_search[d]:n_search[d], Val(N)))
        o = SVector{N, Int}(Tuple(I))
        lexicographically_positive(o) || continue
        min_cell_distance(C, o) > cutoff * (1 + sqrt(eps(T))) && continue
        push!(offsets, o)
    end
    return Grid{N, T, N * N}(M, inv(M), cutoff, cutoff^2, n_cells, n_search, offsets)
end

mutable struct CellList{N, T, NN}
    grid::Grid{N, T, NN}
    particles::Vector{Particle{N, T}}     # sorted by cell
    shifts::Vector{SVector{N, Int32}}     # wrap shifts, sorted as particles
    cell_start::Vector{Int}               # particles of cell c: cell_start[c]:cell_start[c+1]-1
    # scratch
    cell_of::Vector{Int}
    wrapped::Vector{SVector{N, T}}
    wshift::Vector{SVector{N, Int32}}
    pos::Vector{Int}
end

function CellList(x::AbstractVector{SVector{N, T}}, grid::Grid{N, T}; parallel::Bool = true) where {N, T}
    n = length(x)
    ncell = prod(grid.n_cells)
    cl = CellList(
        grid, Vector{Particle{N, T}}(undef, n), Vector{SVector{N, Int32}}(undef, n),
        zeros(Int, ncell + 1), Vector{Int}(undef, n), Vector{SVector{N, T}}(undef, n),
        Vector{SVector{N, Int32}}(undef, n), Vector{Int}(undef, ncell + 1)
    )
    return build!(cl, x; parallel)
end

@inline function cell_linear_index(c::SVector{N, Int}, n_cells::SVector{N, Int}) where {N}
    li = c[N]
    for d in (N - 1):-1:1
        li = li * n_cells[d] + c[d]
    end
    return li + 1
end

@inline function cell_cartesian_index(li::Int, n_cells::SVector{N, Int}) where {N}
    r = li - 1
    return SVector{N, Int}(
        ntuple(Val(N)) do d
            q, m = divrem(r, n_cells[d])
            r = q
            m
        end
    )
end

@inline function locate(xi::SVector{N, T}, grid::Grid{N, T}) where {N, T}
    f = grid.invM * xi
    fl = floor.(f)
    fw = f - fl
    at_upper = fw .== one(T)
    fw = ifelse.(at_upper, zero(T), fw)
    fl = ifelse.(at_upper, fl .+ one(T), fl)
    c = min.(floor.(Int, fw .* grid.n_cells), grid.n_cells .- 1)
    return cell_linear_index(c, grid.n_cells), grid.M * fw, SVector{N, Int32}(Int32.(-fl))
end

function build!(cl::CellList{N, T}, x; parallel::Bool = true) where {N, T}
    grid = cl.grid
    n = length(x)
    # 1. wrap and locate (embarrassingly parallel)
    if parallel && nthreads() > 1
        @threads :static for i in 1:n
            cl.cell_of[i], cl.wrapped[i], cl.wshift[i] = locate(x[i], grid)
        end
    else
        for i in 1:n
            cl.cell_of[i], cl.wrapped[i], cl.wshift[i] = locate(x[i], grid)
        end
    end
    # 2. counting sort by cell
    cs = cl.cell_start
    fill!(cs, 0)
    for i in 1:n
        cs[cl.cell_of[i] + 1] += 1
    end
    cs[1] = 1
    for c in 2:length(cs)
        cs[c] += cs[c - 1]
    end
    copyto!(cl.pos, cs)
    for i in 1:n
        c = cl.cell_of[i]
        p = cl.pos[c]
        cl.pos[c] = p + 1
        cl.particles[p] = Particle(i, cl.wrapped[i])
        cl.shifts[p] = cl.wshift[i]
    end
    return cl
end

#
# Pair loop for the particles of cell `c`. The callback receives
# (i, j, xi, xj, d2, shift, output), with xj the image of j within the cutoff of xi,
# and shift the lattice translation such that x[j] + M*shift - x[i] == xj - xi.
#
@inline function cell_pairs!(f::F, output, cl::CellList{N, T}, c::Int) where {F, N, T}
    grid = cl.grid
    (; M, cutoff_sqr, n_cells, offsets) = grid
    particles = cl.particles
    shifts = cl.shifts
    cs = cl.cell_start
    first_i, last_i = cs[c], cs[c + 1] - 1
    first_i > last_i && return output
    # same cell
    @inbounds for a in first_i:last_i
        pa = particles[a]
        for b in (a + 1):last_i
            pb = particles[b]
            d2 = sum(abs2, pb.x - pa.x)
            if d2 <= cutoff_sqr
                S = SVector{N, Int}(shifts[b]) - SVector{N, Int}(shifts[a])
                output = f(pa.index, pb.index, pa.x, pb.x, d2, S, output)
            end
        end
    end
    # neighbor cells, through the half set of offsets
    ci = cell_cartesian_index(c, n_cells)
    @inbounds for o in offsets
        t = ci + o
        q = fld.(t, n_cells)
        cj = cell_linear_index(t - q .* n_cells, n_cells)
        first_j, last_j = cs[cj], cs[cj + 1] - 1
        first_j > last_j && continue
        cshift = M * q
        for a in first_i:last_i
            pa = particles[a]
            xa = pa.x - cshift # so that xb - xa is the vector to the image of b
            for b in first_j:last_j
                pb = particles[b]
                d2 = sum(abs2, pb.x - xa)
                if d2 <= cutoff_sqr
                    S = q + SVector{N, Int}(shifts[b]) - SVector{N, Int}(shifts[a])
                    output = f(pa.index, pb.index, pa.x, pb.x + cshift, d2, S, output)
                end
            end
        end
    end
    return output
end

function map_pairs(f::F, output, cl::CellList; parallel::Bool = true, reducer = +, copy_output = deepcopy) where {F}
    ncell = length(cl.cell_start) - 1
    if !parallel || nthreads() == 1
        for c in 1:ncell
            output = cell_pairs!(f, output, cl, c)
        end
        return output
    end
    nt = nthreads()
    nchunks = min(ncell, 16 * nt)
    chunks = [round(Int, (k - 1) * ncell / nchunks) + 1:round(Int, k * ncell / nchunks) for k in 1:nchunks]
    next = Atomic{Int}(1)
    outputs = [copy_output(output) for _ in 1:nt]
    tasks = Vector{Task}(undef, nt)
    for it in 1:nt
        tasks[it] = @spawn begin
            out = outputs[$it]
            while true
                k = atomic_add!(next, 1)
                k > nchunks && break
                for c in chunks[k]
                    out = cell_pairs!(f, out, cl, c)
                end
            end
            outputs[$it] = out
        end
    end
    foreach(wait, tasks)
    for it in 1:nt
        output = reducer(output, outputs[it])
    end
    return output
end

# Neighbor list (i, j, d), as in CellListMap.neighborlist. Per-thread lists are
# pre-sized from the expected number of pairs (uniform density), and concatenated
# into a single list sized once.
mutable struct PairList{T}
    n::Int
    list::Vector{Tuple{Int, Int, T}}
end
@inline function push_pair(i, j, xi, xj, d2, S, pl::PairList)
    pl.n += 1
    if pl.n > length(pl.list)
        push!(pl.list, (i, j, sqrt(d2)))
    else
        @inbounds pl.list[pl.n] = (i, j, sqrt(d2))
    end
    return pl
end

function neighborlist(cl::CellList{N, T}; parallel::Bool = true) where {N, T}
    n = length(cl.particles)
    vol = abs(det(cl.grid.M))
    sphere = (N == 2 ? π : 4π / 3) * cl.grid.cutoff^N
    expected = round(Int, 1.1 * n * n / 2 * sphere / vol)
    nt = (parallel && nthreads() > 1) ? nthreads() : 1
    if nt == 1
        pl = map_pairs(push_pair, PairList(0, Vector{Tuple{Int, Int, T}}(undef, expected)), cl; parallel = false)
        resize!(pl.list, pl.n)
        return pl.list
    end
    lists = PairList{T}[]
    reducer(a, b) = (push!(lists, b); a)
    map_pairs(
        push_pair, PairList(0, Tuple{Int, Int, T}[]), cl; parallel = true, reducer,
        copy_output = _ -> PairList(0, Vector{Tuple{Int, Int, T}}(undef, cld(expected, nt)))
    )
    total = sum(l -> l.n, lists)
    out = Vector{Tuple{Int, Int, T}}(undef, total)
    starts = cumsum([1; [l.n for l in lists[1:(end - 1)]]])
    @threads for k in eachindex(lists)
        copyto!(out, starts[k], lists[k].list, 1, lists[k].n)
    end
    return out
end

end
